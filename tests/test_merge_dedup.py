"""Exact external-sort merge against the sequence oracle and Bloom backend."""
import os
from pathlib import Path
import random
import re
import resource
import signal
import subprocess
import tempfile
import unittest

from test_query_merge import canonical, oracle, pack

CORE = os.environ.get('MEGAHIT_TEST_CORE')


@unittest.skipUnless(CORE, 'set MEGAHIT_TEST_CORE to the built core')
class ExactMergeDedupTest(unittest.TestCase):
    def setUp(self):
        self.workspace = tempfile.TemporaryDirectory(prefix='exact merge ')
        self.addCleanup(self.workspace.cleanup)
        self.root = Path(self.workspace.name)

    def fixture(self, overlap, k, shards, contigs, long_contigs):
        self.shards, self.contigs, self.long_contigs = shards, contigs, long_contigs
        prefix = self.root / 'reads'
        for i, shard in enumerate(shards):
            Path(str(prefix) + '.edges.' + str(i)).write_bytes(
                b''.join(pack(sequence, count) for sequence, count in shard))
        words = (2 * (k + 1) + 16 + 31) // 32
        Path(str(prefix) + '.edges.info').write_text(
            'kmer_size {}\nwords_per_edge {}\nnum_files {}\nnum_buckets 0\n'
            'num_edges {}\nis_sorted 0\n'.format(k, words, len(shards), sum(map(len, shards))))
        for name, sequences in [('short.fa', contigs), ('long.fa', long_contigs)]:
            (self.root / name).write_text(''.join('>c{}\n{}\n'.format(i, seq) for i, seq in enumerate(sequences)))
        self.overlap, self.k = overlap, k
        self.command = [CORE, 'prefixbf', '--mode', 'merge', '--long_prefix', str(prefix),
                        '--short_k', str(overlap), '--long_k', str(k),
                        '--contig', str(self.root / 'short.fa'), '--long_contig', str(self.root / 'long.fa'),
                        '--max_bloom_gb', '0.02', '--human_readable', '0']

    def run_merge(self, label, dedup='sort', threads=4, index='hash', memory='0.125', extra=(), **kwargs):
        output = self.root / label
        result = subprocess.run(self.command + ['--output_prefix', str(output), '-t', str(threads),
                                '--merge_index', index, '--merge_dedup', dedup,
                                '--sort_memory_mb', memory] + list(extra),
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, **kwargs)
        return output, result

    def check_output(self, output, threads):
        expected = oracle(self.shards, self.contigs, self.long_contigs,
                          self.overlap, self.k + 1, threads)
        for i, data in enumerate(expected):
            self.assertEqual(Path(str(output) + '.edges.' + str(i)).read_bytes(), data)
        self.assertFalse(list(self.root.glob(output.name + '.sort-*')))

    def test_core_defaults_to_exact_sort(self):
        self.fixture(3, 5, [[('AACCGG', 7)]], ['CAATTT'], ['AACCGGTTAACC'])
        output = self.root / 'default'
        result = subprocess.run(self.command + ['--output_prefix', str(output), '-t', '2'],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        self.assertIn('no merge Bloom filters', result.stderr)
        self.check_output(output, 2)

    def test_long_only_merge_without_short_contig_input(self):
        self.fixture(3, 5, [[('AACCGG', 7), ('AACCGG', 19)], []], [],
                     ['AACCGGTTAACC', 'AACCGGTTAACC', 'GGTTAACCGGTT'])
        position = self.command.index('--contig')
        del self.command[position:position + 2]
        output, result = self.run_merge('long-only', threads=4)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('Merge index: disabled', result.stderr)
        self.check_output(output, 4)
        # Both original base records survive unchanged; generated copies do not.
        data = Path(str(output) + '.edges.0').read_bytes()
        self.assertTrue(data.startswith(pack('AACCGG', 7) + pack('AACCGG', 19)))

    def test_terminal_multiplicity_changes_only_novel_edges(self):
        self.fixture(3, 5, [[('AACCGG', 7)]], [], ['AACCGGTTAACC'] * 2)
        for dedup in ('sort', 'bloom'):
            output, result = self.run_merge('multi2-' + dedup, dedup=dedup,
                                            extra=['--contig_edge_multiplicity', '2'])
            self.assertEqual(result.returncode, 0, result.stderr)
            expected = oracle(self.shards, self.contigs, self.long_contigs, 3, 6, 4)
            for i, shard in enumerate(expected):
                rows = [shard[j:j + 4] for j in range(0, len(shard), 4)]
                # The first record in shard zero is the unchanged read-counted edge.
                start = 1 if i == 0 else 0
                expected_bytes = b''.join(rows[:start] + [b'\x02\x00' + row[2:] for row in rows[start:]])
                self.assertEqual(Path(str(output) + '.edges.' + str(i)).read_bytes(), expected_bytes)

    def test_word_widths_orientations_multiplicity_and_order(self):
        for overlap, k in [(3, 5), (21, 29), (31, 32), (39, 59), (63, 79), (95, 99),
                            (127, 141), (159, 179), (191, 199), (223, 239), (254, 255)]:
            rng = random.Random(overlap * 1000 + k)
            dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
            edges = [canonical(dna(k + 1)) for _ in range(12)]
            shards = [[(seq, i + 7) for i, seq in enumerate(edges[:6])],
                      [(seq, i + 19) for i, seq in enumerate(edges[6:])]]
            desired = [edge[-overlap:] + dna(3) for edge in edges[:4]]
            desired += [dna(3) + edge[:overlap] for edge in edges[4:8]]
            self.fixture(overlap, k, shards, [seq[::-1] for seq in desired],
                         [dna(3 * (k + 1)), 'N' + dna(k + 1)])
            for index, threads in [('hash', 1), ('hash', 4), ('fm', 4)]:
                with self.subTest(k=k, index=index, threads=threads):
                    output, result = self.run_merge('sort-{}-{}-{}'.format(k, index, threads),
                                                    index=index, threads=threads)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertIn('no merge Bloom filters', result.stderr)
                    self.check_output(output, threads)

    def test_deferred_merge_preserves_duplicate_candidates(self):
        self.fixture(3, 5, [[('AACCGG', 7)]], [], ['AACCGGTTAACC'] * 50)
        for threads in (1, 4, 12):
            output, result = self.run_merge('deferred-' + str(threads), dedup='none', threads=threads,
                                            extra=['--contig_edge_multiplicity', '2'])
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('deferred until reduction', result.stderr)
            self.assertIn('Long-contig extraction:', result.stderr)
            rows = []
            for path in self.root.glob(output.name + '.edges.[0-9]*'):
                data = path.read_bytes()
                shard_id = int(path.suffix[1:])
                if shard_id < threads:
                    self.assertEqual(data, pack('AACCGG', 7) if shard_id == 0 else b'')
                else:
                    self.assertTrue(all(data[i:i + 2] == b'\x02\x00' for i in range(0, len(data), 4)))
                rows.extend(data[i:i + 4] for i in range(0, len(data), 4))
            expected = [pack('AACCGG', 7)]
            for seq in self.long_contigs:
                seq = seq[::-1]
                expected += [pack(canonical(seq[i:i + 6]), 2) for i in range(len(seq) - 5)]
            self.assertCountEqual(rows, expected)
            self.assertFalse(list(self.root.glob(output.name + '.sort-*')))

    def test_duplicate_base_edges_palindromes_and_empty_inputs(self):
        for shards, short, long in [([[('AACGTT', 7), ('AACCGG', 2)], [('AACCGG', 19)]],
                                     ['CAATTT', 'CAATTT'], ['ACGTACGTACGT', 'ACGTNNACGT']),
                                    ([[], []], [], []),
                                    ([[], []], ['NNNNNN'], ['AACCGGTTAACC'])]:
            self.fixture(3, 5, shards, short, long)
            output, result = self.run_merge('duplicates')
            self.assertEqual(result.returncode, 0, result.stderr)
            self.check_output(output, 4)

    def test_forced_multipass_runs_and_candidate_order_sort(self):
        rng = random.Random(38152)
        dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
        records = [(canonical(dna(30)), i % 100 + 1) for i in range(60000)]
        self.fixture(21, 29, [records], [], [dna(10000)])
        output, result = self.run_merge('multipass', memory='0.0625')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.check_output(output, 4)
        match = re.search(r'(\d+) runs, (\d+) temporary bytes written', result.stderr)
        self.assertIsNotNone(match, result.stderr)
        self.assertGreater(int(match.group(1)), 32)
        self.assertGreater(int(match.group(2)), 0)
        bloom, result = self.run_merge('bloom', dedup='bloom')
        self.assertEqual(result.returncode, 0, result.stderr)
        for i in range(4):
            self.assertEqual(Path(str(output) + '.edges.' + str(i)).read_bytes(),
                             Path(str(bloom) + '.edges.' + str(i)).read_bytes())

    def test_write_failure_cleans_temporary_runs_and_does_not_publish_metadata(self):
        rng = random.Random(71)
        dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
        self.fixture(21, 29, [[(canonical(dna(30)), 7) for _ in range(10000)]], [], [])
        def limit_output():
            signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
            resource.setrlimit(resource.RLIMIT_FSIZE, (4096, 4096))
        output, result = self.run_merge('disk-failure', threads=1, memory='0.0625', preexec_fn=limit_output)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Failed to write sorted run', result.stderr)
        self.assertFalse(Path(str(output) + '.edges.info').exists())
        self.assertFalse(list(self.root.glob(output.name + '.sort-*')))

    def test_invalid_sort_settings_and_unsupported_reports_fail(self):
        self.fixture(3, 5, [[('AACCGG', 7)]], [], [])
        for extra in [('--sort_memory_mb', '0'), ('--merge_dedup', 'invalid'), ('--human_readable', '1')]:
            output, result = self.run_merge('invalid', extra=extra)
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(Path(str(output) + '.edges.info').exists())
