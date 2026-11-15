"""Packed query-first merge against an independent sequence oracle."""
import os
import json
from pathlib import Path
import random
import struct
import subprocess
import tempfile
import unittest


CORE = os.environ.get('MEGAHIT_TEST_CORE')
BASELINE = os.environ.get('MEGAHIT_TEST_BASELINE_CORE')


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def canonical(sequence):
    sequence = sequence.upper()
    return min(sequence, reverse_complement(sequence))


def pack(sequence, multiplicity):
    words = [0] * ((len(sequence) * 2 + 16 + 31) // 32)
    for i, base in enumerate(sequence):
        words[i // 16] |= 'ACGT'.index(base) << (30 - 2 * (i % 16))
    words[-1] |= multiplicity
    return struct.pack('=' + 'I' * len(words), *words)


def oracle(shards, contigs, long_contigs, overlap, length, threads):
    base = {canonical(sequence) for shard in shards for sequence, _ in shard}
    orientations = []
    output = [[] for _ in range(threads)]
    for file_id, shard in enumerate(shards):
        for sequence, multiplicity in shard:
            output[file_id % threads].append(pack(sequence, multiplicity))
            orientations.append(sequence)
            rc = reverse_complement(sequence)
            if rc != sequence:
                orientations.append(rc)
    emitted = set()
    added = 0

    def emit(sequence, first, last):
        nonlocal added
        for i in range(first, min(last, len(sequence) - length + 1)):
            window = sequence[i:i + length].upper()
            if any(base not in 'ACGT' for base in window):
                continue
            key = canonical(window)
            if key in base or key in emitted:
                continue
            emitted.add(key)
            output[added % threads].append(pack(key, 50))
            added += 1

    for sequence in contigs:
        sequence = sequence[::-1]
        if len(sequence) < overlap:
            continue
        prefix, suffix = sequence[:overlap].upper(), sequence[-overlap:].upper()
        if any(base not in 'ACGT' for base in prefix + suffix):
            continue
        extra = len(sequence) - overlap
        for edge in orientations:
            if edge[-overlap:] == prefix:
                emit(edge + sequence[overlap:], 1, extra + 1)
        for edge in orientations:
            if edge[:overlap] == suffix:
                emit(sequence[:-overlap] + edge, 0, extra)
    for sequence in long_contigs:
        sequence = sequence[::-1]
        emit(sequence, 0, len(sequence) - length + 1)
    return [b''.join(shard) for shard in output]


@unittest.skipUnless(CORE, 'set MEGAHIT_TEST_CORE to the newly built core')
class QueryMergeTest(unittest.TestCase):
    def check_merge(self, overlap, k, shards, contigs, long_contigs, threads=3, reports=False,
                    benchmark=False):
        with tempfile.TemporaryDirectory(prefix='megahit-query-merge-') as directory:
            root = Path(directory)
            prefix = root / 'input'
            length = k + 1
            words = (length * 2 + 16 + 31) // 32
            for i, shard in enumerate(shards):
                Path(str(prefix) + '.edges.' + str(i)).write_bytes(
                    b''.join(pack(sequence, multiplicity) for sequence, multiplicity in shard))
            Path(str(prefix) + '.edges.info').write_text(
                'kmer_size {}\nwords_per_edge {}\nnum_files {}\nnum_buckets 0\n'
                'num_edges {}\nis_sorted 0\n'.format(k, words, len(shards),
                                                    sum(map(len, shards))))
            for filename, sequences in [('short.fa', contigs), ('long.fa', long_contigs)]:
                (root / filename).write_text(''.join('>c{}\n{}\n'.format(i, seq)
                                                     for i, seq in enumerate(sequences)))
            command = ['prefixbf', '--mode', 'merge', '--short_k', str(overlap),
                       '--long_k', str(k), '--long_prefix', str(prefix),
                       '--contig', str(root / 'short.fa'), '--long_contig', str(root / 'long.fa'),
                       '--num_cpu_threads', str(threads), '--max_bloom_gb', '0.01',
                       '--human_readable', '1' if reports else '0']
            expected = oracle(shards, contigs, long_contigs, overlap, length, threads)
            benchmark_counts = None
            if BASELINE:
                baseline = root / 'old'
                subprocess.run([BASELINE] + command + ['--output_prefix', str(baseline)],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            for backend in ('fm', 'hash'):
                output = root / ('new-' + backend)
                result = subprocess.run([CORE] + command + ['--output_prefix', str(output),
                                         '--merge_index', backend, '--merge_dedup',
                                         'bloom' if reports else 'sort'], stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE, check=True, text=True)
                self.assertIn('Merge index: ' + backend, result.stderr)
                if benchmark:
                    bench = list(command)
                    bench[bench.index('merge')] = 'index-bench'
                    measured = subprocess.run([CORE] + bench + ['--merge_index', backend],
                                              stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                              text=True, check=True)
                    counts = json.loads(measured.stdout)
                    for field in ('backend', 'threads', 'build_seconds', 'lookup_seconds', 'fm_bytes'):
                        counts.pop(field)
                    self.assertEqual(counts['records'], sum(map(len, shards)))
                    if benchmark_counts is not None:
                        self.assertEqual(counts, benchmark_counts)
                    benchmark_counts = counts
                for i, data in enumerate(expected):
                    self.assertEqual(Path(str(output) + '.edges.' + str(i)).read_bytes(), data)
                if BASELINE:
                    for suffix in ['.edges.info'] + ['.edges.' + str(i) for i in range(threads)]:
                        self.assertEqual(Path(str(output) + suffix).read_bytes(),
                                         Path(str(baseline) + suffix).read_bytes())
                    if reports:
                        self.assertEqual(Path(str(output) + '.contig_edges.txt').read_bytes(),
                                         Path(str(baseline) + '.contig_edges.txt').read_bytes())

    def test_word_boundaries_and_maximum_k(self):
        for overlap, k in [(3, 5), (21, 29), (21, 31), (31, 32), (32, 39), (33, 39),
                           (39, 59), (63, 64), (64, 79), (65, 99), (127, 141), (255, 255)]:
            with self.subTest(overlap=overlap, k=k):
                rng = random.Random(overlap * 1000 + k)
                dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
                edges = [canonical(dna(k + 1)) for _ in range(18)]
                shards = [[(seq, i + 1) for i, seq in enumerate(edges[:9])],
                          [(seq, i + 19) for i, seq in enumerate(edges[9:])]]
                desired = [edge[-overlap:] + dna(5) for edge in edges[:4]]
                desired += [dna(5) + reverse_complement(edge)[:overlap] for edge in edges[4:8]]
                desired += [desired[0], desired[1].lower(), 'N' + desired[2][1:],
                            edges[0][-overlap:] + 'N' + dna(overlap) + edges[1][:overlap]]
                self.check_merge(overlap, k, shards, [seq[::-1] for seq in desired],
                                 [dna(3 * (k + 1)), dna(k + 1) + 'NN' + dna(k + 1)],
                                 reports=True)

    def test_duplicate_edges_palindromes_and_repeated_queries(self):
        self.check_merge(3, 5, [[('AACGTT', 7), ('AACCGG', 9)],
                               [('AACCGG', 19), ('GGTTAA', 3)]],
                         ['CAATTT', 'CAATTT', 'CCGAAA', 'AAAAAA', 'NNNNNN'],
                         ['ACGTACGTACGT', 'AAAAAA', 'ACGTNNACGTACGT'], reports=True)

    def test_no_queries_and_empty_input_shards(self):
        self.check_merge(3, 5, [[], [('AACCGG', 9)], []], [], ['ACGTACGTACGT'])
        self.check_merge(3, 5, [[], []], ['CAATTT'], ['ACGTACGTACGT'])

    def test_invalid_end_disables_both_directions(self):
        self.check_merge(3, 5, [[('AACCGG', 9)]],
                         ['AAANCC', 'CCNAAA', 'A', 'CCG'], ['acgtnnacgtacgt'])

    def test_parallel_scan_preserves_shard_ownership_and_match_order(self):
        rng = random.Random(207)
        dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
        for overlap, k in [(21, 29), (39, 59)]:
            edges = [canonical(dna(k + 1)) for _ in range(12)]
            shards = [[(edges[(i + j) % len(edges)], j + 1) for j in range(16)]
                      for i in range(8)]
            # Empty files and fewer/more workers than files exercise stream ownership.
            shards[3] = []
            desired = [edge[-overlap:] + dna(5) for edge in edges[:6]]
            desired += [dna(5) + edge[:overlap] for edge in edges[6:]]
            for threads in [1, 2, 3, 4, 8, 12]:
                with self.subTest(overlap=overlap, threads=threads):
                    self.check_merge(overlap, k, shards, [seq[::-1] for seq in desired],
                                     [dna(4 * (k + 1)), 'A' * (k + 10)], threads=threads)

    def test_large_dictionary_and_seed_collisions(self):
        rng = random.Random(925)
        dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
        edges = [canonical(dna(30)) for _ in range(150)]
        desired = [edge[-21:] + dna(5) for edge in edges[:80]]
        desired += [dna(5) + edge[:21] for edge in edges[80:]]
        desired += [dna(21) + 'ACGTACGT' for _ in range(6000)]
        desired += ['A' * 21 + dna(5)] * 100
        self.check_merge(21, 29, [[(seq, 3) for seq in edges[:75]],
                                  [(seq, 9) for seq in edges[75:]]],
                         [seq[::-1] for seq in desired], [], threads=4, benchmark=True)

    def test_single_large_shard_parallel_search_preserves_order(self):
        rng = random.Random(3932)
        dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
        edge = canonical(dna(30))
        records = [(edge, i % 100 + 1) for i in range(100003)]
        contig = edge[-21:] + dna(3)
        self.check_merge(21, 29, [records], [contig[::-1]], [], threads=4)
