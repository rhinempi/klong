"""Reduction against an exact sequence oracle, including packed-word boundaries."""
import collections
import os
from pathlib import Path
import random
import re
import subprocess
import tempfile
import unittest

from test_query_merge import canonical, pack, reverse_complement

CORE = os.environ.get('MEGAHIT_TEST_CORE')
BASELINE = os.environ.get('MEGAHIT_TEST_BASELINE_CORE')


@unittest.skipUnless(CORE, 'set MEGAHIT_TEST_CORE to the newly built core')
class ReductionTest(unittest.TestCase):
    def setUp(self):
        self.workspace = tempfile.TemporaryDirectory(prefix='packed reduction ')
        self.addCleanup(self.workspace.cleanup)
        self.root = Path(self.workspace.name)

    def write_edges(self, name, k, shards, sorted_input=False):
        prefix = self.root / name
        words = (2 * (k + 1) + 16 + 31) // 32
        buckets = []
        for file_id, shard in enumerate(shards):
            Path(str(prefix) + '.edges.' + str(file_id)).write_bytes(
                b''.join(pack(sequence, multiplicity) for sequence, multiplicity in shard))
            if sorted_input:
                buckets.append('{} {} 0 {}\n'.format(file_id, file_id, len(shard)))
        Path(str(prefix) + '.edges.info').write_text(
            'kmer_size {}\nwords_per_edge {}\nnum_files {}\nnum_buckets {}\n'
            'num_edges {}\nis_sorted {}\n'.format(k, words, len(shards), len(buckets),
                                                  sum(map(len, shards)), int(sorted_input))
            + ''.join(buckets))
        return prefix, words

    def records(self, prefix, words):
        result = collections.Counter()
        for path in prefix.parent.glob(prefix.name + '.edges.[0-9]*'):
            if not path.name.rsplit('.', 1)[-1].isdigit():
                continue
            data = path.read_bytes()
            stride = words * 4
            self.assertEqual(len(data) % stride, 0)
            result.update(data[i:i + stride] for i in range(0, len(data), stride))
        return result

    def reduce(self, label, short_k, long_k, short, longer, threads=3, extra=(), binary=None):
        output = self.root / label
        command = [binary or CORE, 'prefixbf', '--mode', 'reduce',
                   '--short_prefix', str(short), '--long_prefix', str(longer),
                   '--short_k', str(short_k), '--long_k', str(long_k),
                   '--output_prefix', str(output), '-t', str(threads),
                   '--max_bloom_gb', '0.02'] + list(extra)
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        return output, result

    def test_word_boundaries_orientations_and_thread_counts(self):
        for short_k, long_k in [(1, 3), (3, 5), (14, 29), (15, 31), (16, 32),
                                (21, 29), (29, 39), (31, 59), (32, 59), (39, 59),
                                (63, 79), (64, 99), (127, 141), (254, 255), (255, 255)]:
            rng = random.Random(short_k * 1000 + long_k)
            dna = lambda length: ''.join(rng.choice('ACGT') for _ in range(length))
            queries = [dna(short_k + 1) for _ in range(30)]
            queries += ['A' * (short_k + 1), 'T' * (short_k + 1), queries[0]]
            edges = []
            for i, sequence in enumerate(queries):
                if i % 3 != 0:
                    edges.append(sequence + dna(long_k - short_k))
                if i % 3 == 2:
                    edges.append(reverse_complement(sequence) + dna(long_k - short_k))
            represented = {sequence[:short_k + 1] for sequence in edges}
            represented.update(reverse_complement(sequence[-short_k - 1:]) for sequence in edges)
            short_shards = [[(sequence, i + 1) for i, sequence in enumerate(queries)], []]
            long_shards = [[(sequence, i + 9) for i, sequence in enumerate(edges)], []]
            short, words = self.write_edges('short', short_k, short_shards)
            longer, _ = self.write_edges('long', long_k, long_shards, sorted_input=True)
            expected = collections.Counter(pack(sequence, count) for sequence, count in short_shards[0]
                                           if sequence not in represented or
                                           reverse_complement(sequence) not in represented)
            for threads in (1, 4):
                with self.subTest(short_k=short_k, long_k=long_k, threads=threads):
                    label = '{}-{}-{}'.format(short_k, long_k, threads)
                    output, _ = self.reduce('new-' + label, short_k, long_k, short, longer, threads,
                                            ['--human_readable', '1'])
                    self.assertEqual(self.records(output, words), expected)
                    if BASELINE:
                        before, _ = self.reduce('old-' + label, short_k, long_k, short, longer, threads,
                                                binary=BASELINE)
                        self.assertEqual(self.records(output, words), self.records(before, words))

    def test_large_shard_crosses_input_and_output_buffers(self):
        sequence = 'ACGTTGCAACCTGATCGTACGA'
        queries = [(sequence, 17)] * 150003
        short, words = self.write_edges('short', 21, [queries])
        longer, _ = self.write_edges('long', 29, [[('A' * 30, 2)]])
        output, result = self.reduce('large', 21, 29, short, longer, threads=4)
        self.assertEqual(self.records(output, words), collections.Counter({pack(sequence, 17): len(queries)}))
        self.assertIn('2 filtering chunks, 4 workers', result.stderr)
        output, _ = self.reduce('large-unique', 21, 29, short, longer, threads=12,
                                extra=['--reduce_dedup', '--sort_memory_mb', '0.0625'])
        self.assertEqual(self.records(output, words), collections.Counter({pack(canonical(sequence), 17): 1}))

    def test_post_reduction_dedup_orientations_counts_and_multipass(self):
        for k, long_k in ((21, 29), (39, 59), (127, 141), (254, 255)):
            rng = random.Random(k)
            dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
            sequences = [canonical(dna(k + 1)) for _ in range(4500)]
            rows = [(seq, 3) for seq in sequences]
            rows += [(reverse_complement(seq), 50) for seq in sequences]
            rows += [(seq, 65535) for seq in sequences[:10]]
            short, words = self.write_edges('duplicates', k, [rows])
            # Remove every orientation of the first sequence, not just its duplicates.
            seq = sequences[0]
            longer, _ = self.write_edges('representatives', long_k,
                                        [[(seq + dna(long_k - k), 7),
                                          (reverse_complement(seq) + dna(long_k - k), 7)]])
            expected = collections.Counter(pack(seq, 65535 if i < 10 else 50)
                                           for i, seq in enumerate(sequences) if i)
            for threads in (1, 4, 12):
                with self.subTest(k=k, threads=threads):
                    output, result = self.reduce('unique-{}-{}'.format(k, threads), k, long_k,
                                                short, longer, threads,
                                                ['--reduce_dedup', '--sort_memory_mb', '0.0625'])
                    self.assertEqual(self.records(output, words), expected)
                    self.assertIn('unique surviving edges', result.stderr)
                    self.assertFalse(list(self.root.glob(output.name + '.sort-*')))

    def test_post_reduction_dedup_empty_survivors(self):
        short, words = self.write_edges('all-short', 21, [[('A' * 22, 3), ('T' * 22, 50)]])
        longer, _ = self.write_edges('all-long', 29, [[('A' * 30, 2)]])
        output, _ = self.reduce('empty-unique', 21, 29, short, longer, 4, ['--reduce_dedup'])
        self.assertFalse(self.records(output, words))

    def test_large_long_shard_uses_multiple_insertion_chunks(self):
        short, words = self.write_edges('short', 21, [[('A' * 22, 17), ('T' * 22, 9)]])
        longer, _ = self.write_edges('long', 29, [[('A' * 30, 2)] * 200003])
        output, result = self.reduce('long-chunks', 21, 29, short, longer, threads=4)
        self.assertIn('3 insertion chunks', result.stderr)
        self.assertFalse(self.records(output, words))
        output, capped = self.reduce('non-power-cap', 21, 29, short, longer, threads=4,
                                     extra=['--max_bloom_gb', '0.012'])
        self.assertFalse(self.records(output, words))
        match = re.search(r'Reduction Bloom: (\d+) bytes.*cap (\d+) bytes', capped.stderr)
        self.assertIsNotNone(match, capped.stderr)
        size, cap = map(int, match.groups())
        self.assertGreater(size, 8 << 20)
        self.assertLessEqual(size, cap)

    def test_memory_budget_and_fpr_are_reported(self):
        short, _ = self.write_edges('short', 21, [[('A' * 22, 1)]])
        longer, _ = self.write_edges('long', 29, [[('A' * 30, 2)] * 6000])
        _, result = self.reduce('budgeted', 21, 29, short, longer, threads=2,
                                extra=['--host_mem', str(6 << 20)])
        match = re.search(r'Reduction Bloom: (\d+) bytes, (\d+) probes, (\d+) active words, '
                          r'estimated FPR ([\deE.+-]+), cap (\d+) bytes', result.stderr)
        self.assertIsNotNone(match, result.stderr)
        size, probes, words, fpr, cap = match.groups()
        self.assertLessEqual(int(size), int(cap))
        self.assertLessEqual(int(cap), (6 << 20) - (4 << 20))
        self.assertEqual(int(size) & (int(size) - 1), 0)
        self.assertEqual(int(words), 1)
        self.assertLessEqual(float(fpr), 1e-15)
        self.assertLess(int(probes), 16)

        _, capped = self.reduce('capped', 21, 29, short, longer, threads=1,
                                extra=['--max_bloom_gb', '0.000001'])
        self.assertIn('WARNING: Bloom memory/hash cap prevents', capped.stderr)

    def test_invalid_budget_or_fpr_fails(self):
        short, _ = self.write_edges('short', 21, [[('A' * 22, 1)]])
        longer, _ = self.write_edges('long', 29, [[('A' * 30, 2)]])
        for option, value in [('--host_mem', '-1'), ('--reduce_fpr', '0'),
                              ('--reduce_fpr', '1'), ('--max_bloom_gb', '-1'),
                              ('--reserve_gb', '-1'), ('--host_mem', '100')]:
            result = subprocess.run([CORE, 'prefixbf', '--short_prefix', str(short),
                                     '--long_prefix', str(longer), '--short_k', '21', '--long_k', '29',
                                     '--output_prefix', str(self.root / 'bad-options'), option, value],
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.assertNotEqual(result.returncode, 0, (option, value))

    def test_empty_inputs(self):
        short, words = self.write_edges('short', 21, [[], []])
        longer, _ = self.write_edges('long', 29, [[]])
        output, _ = self.reduce('empty', 21, 29, short, longer)
        self.assertFalse(self.records(output, words))

    def test_invalid_or_missing_shard_fails(self):
        short, _ = self.write_edges('short', 21, [[('A' * 22, 1)]])
        longer, _ = self.write_edges('long', 29, [[('A' * 30, 1)]])
        for content in (b'broken', None):
            path = Path(str(short) + '.edges.0')
            if content is None:
                path.unlink()
            else:
                path.write_bytes(content)
            result = subprocess.run([CORE, 'prefixbf', '--short_prefix', str(short),
                                     '--long_prefix', str(longer), '--short_k', '21', '--long_k', '29',
                                     '--output_prefix', str(self.root / 'invalid'), '--max_bloom_gb', '0.01'],
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.assertNotEqual(result.returncode, 0)
