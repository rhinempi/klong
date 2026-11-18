"""Differential assembly tests against the preserved pre-optimization core."""
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
BASELINE = os.environ.get('MEGAHIT_TEST_ASSEMBLY_BASELINE')


def records(path):
    header, sequence = None, []
    with Path(path).open() as source:
        for line in source:
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(sequence)
                header, sequence = line[1:].strip(), []
            else:
                sequence.append(line.strip())
    if header is not None:
        yield header, ''.join(sequence)


def semantic_record(header, sequence):
    metadata = dict(field.split('=', 1) for field in header.split()[1:] if '=' in field)
    flag = int(metadata['flag'])
    if flag & 2:
        k = int(header.split('_', 1)[0][1:])
        core = sequence[:-k]
        sequence = min(least_rotation(core), least_rotation(reverse_complement(core)))
    else:
        sequence = min(sequence, reverse_complement(sequence))
    return '{}\t{}\t{}\t{}'.format(flag, metadata['multi'], len(sequence), sequence)


def least_rotation(sequence):
    doubled = sequence + sequence
    n, i, j, offset = len(sequence), 0, 1, 0
    while i < n and j < n and offset < n:
        a, b = doubled[i + offset], doubled[j + offset]
        if a == b:
            offset += 1
            continue
        if a > b:
            i += offset + 1
            if i == j:
                i += 1
        else:
            j += offset + 1
            if i == j:
                j += 1
        offset = 0
    start = min(i, j)
    return doubled[start:start + n]


class AssemblyNormalizationTest(unittest.TestCase):
    def test_rotation_and_reverse_complement_normalization(self):
        rng = random.Random(1751)
        for n in range(1, 40):
            sequence = ''.join(rng.choice('ACGT') for _ in range(n))
            self.assertEqual(least_rotation(sequence),
                             min(sequence[i:] + sequence[:i] for i in range(n)))
        first = semantic_record('k3_1 flag=3 multi=7.0000 len=8', 'AACGTAAC')
        rotated = semantic_record('k3_9 flag=3 multi=7.0000 len=8', 'CGTAACGT')
        self.assertEqual(first, rotated)


@unittest.skipUnless(CORE and BASELINE, 'set MEGAHIT_TEST_CORE and MEGAHIT_TEST_ASSEMBLY_BASELINE')
class AssemblyOptimizationsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.workspace = tempfile.TemporaryDirectory(prefix='assembly-optimizations-')
        cls.root = Path(cls.workspace.name)
        rng = random.Random(7391)
        dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
        genome = dna(650)
        variants = [genome[:150] + ('A' if genome[150] != 'A' else 'C') + genome[151:],
                    genome[:300] + 'ACG' + genome[300:], genome[:430] + genome[433:]]
        paths = [(genome, 30)] + [(seq, 2 + i) for i, seq in enumerate(variants)]
        paths += [(genome[70:120] + dna(13), 1), (dna(17) + genome[500:550], 1)]
        paths += [(dna(45), count) for count in (1, 2, 255, 256, 65535)]
        circular = dna(37)
        paths.append((circular + circular[:21], 7))
        edges = collections.Counter()
        for sequence, count in paths:
            for i in range(len(sequence) - 21):
                edges[canonical(sequence[i:i + 22])] += count
        prefix = cls.root / 'edges'
        Path(str(prefix) + '.edges.0').write_bytes(b''.join(pack(seq, min(count, 65535))
                                                         for seq, count in edges.items()))
        Path(str(prefix) + '.edges.info').write_text(
            'kmer_size 21\nwords_per_edge 2\nnum_files 1\nnum_buckets 0\n'
            'num_edges {}\nis_sorted 0\n'.format(len(edges)))
        cls.graph = cls.root / 'graph'
        subprocess.run([CORE, 'seq2sdbg', '--input_prefix', str(prefix), '-k', '21',
                        '--output_prefix', str(cls.graph), '--host_mem', '200000000',
                        '--mem_flag', '1', '--num_cpu_threads', '4'],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, timeout=60)

    @classmethod
    def tearDownClass(cls):
        cls.workspace.cleanup()

    def test_cleaning_modes_inference_statistics_and_outputs(self):
        for threads in (1, 4):
            for prune, depth, final in [(0, 0, False), (0, 7, False), (1, 0, False),
                                        (2, 2, False), (2, 0, True), (3, 2, True)]:
                with self.subTest(threads=threads, prune=prune, depth=depth, final=final):
                    args = ['assemble', '-s', str(self.graph), '-t', str(threads),
                            '--min_standalone', '0', '--prune_level', str(prune),
                            '--min_depth', str(depth), '--bubble_level', '2',
                            '--merge_len', '20', '--merge_similar', '0.95',
                            '--max_tip_len', '42', '--cleaning_rounds', '3',
                            '--disconnect_ratio', '0', '--careful_bubble']
                    if final:
                        args.append('--is_final_round')
                    results = []
                    for label, binary in [('before', BASELINE), ('after', CORE)]:
                        output = self.root / '{}-{}-{}-{}-{}'.format(label, threads, prune, depth, final)
                        result = subprocess.run([binary] + args + ['-o', str(output)],
                                                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                                check=True, text=True, timeout=60)
                        outputs = {}
                        for suffix in ('contigs.fa', 'final.contigs.fa', 'bubble_seq.fa', 'addi.fa'):
                            path = Path(str(output) + '.' + suffix)
                            if path.exists():
                                outputs[suffix] = collections.Counter(semantic_record(*row) for row in records(path))
                                outputs[suffix + '.info'] = Path(str(path) + '.info').read_bytes()
                        results.append((result.stderr, outputs))
                    self.assertEqual(results[0][1], results[1][1])
                    stats = lambda log: re.findall(r'Max: .*total size: \d+,', log)
                    self.assertEqual(stats(results[0][0]), stats(results[1][0]))
                    inferred = lambda log: re.findall(r'min depth set to ([\d.]+)', log)
                    if prune and depth <= 0:
                        self.assertEqual(inferred(results[0][0]), inferred(results[1][0]))
                        self.assertTrue(inferred(results[1][0]))
                    elif not prune:
                        self.assertFalse(inferred(results[1][0]))

    def test_existing_weighted_cycle_coverage_depends_on_start_edge(self):
        rng = random.Random(4011)
        cycle = ''.join(rng.choice('ACGT') for _ in range(37))
        prefix, graph = self.root / 'weighted-cycle', self.root / 'weighted-graph'
        edges = {canonical((cycle + cycle)[i:i + 22]): i + 2 for i in range(37)}
        self.assertEqual(len(edges), 37)
        Path(str(prefix) + '.edges.0').write_bytes(b''.join(pack(seq, count) for seq, count in edges.items()))
        Path(str(prefix) + '.edges.info').write_text(
            'kmer_size 21\nwords_per_edge 2\nnum_files 1\nnum_buckets 0\nnum_edges 37\nis_sorted 0\n')
        subprocess.run([CORE, 'seq2sdbg', '--input_prefix', str(prefix), '-k', '21',
                        '--output_prefix', str(graph), '--host_mem', '200000000',
                        '--mem_flag', '1', '--num_cpu_threads', '4'],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, timeout=60)
        # The unchanged loop constructor counts its start edge twice.
        allowed = {'{:.4f}'.format((sum(edges.values()) + count) / 37) for count in edges.values()}
        observed = {}
        for label, binary in [('before', BASELINE), ('after', CORE)]:
            values = set()
            for repeat in range(12):
                output = self.root / 'cycle-output'
                subprocess.run([binary, 'assemble', '-s', str(graph), '-o', str(output), '-t', '8',
                                '--prune_level', '0', '--max_tip_len', '0', '--bubble_level', '0',
                                '--merge_len', '0', '--cleaning_rounds', '1', '--disconnect_ratio', '0'],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, timeout=60)
                rows = list(records(Path(str(output) + '.contigs.fa')))
                self.assertEqual(len(rows), 1)
                header, sequence = rows[0]
                self.assertIn('flag=3 ', header)
                self.assertEqual(semantic_record(header, sequence).split('\t')[-1],
                                 min(least_rotation(cycle[::-1]), least_rotation(reverse_complement(cycle[::-1]))))
                values.add(re.search(r'multi=([\d.]+)', header)[1])
            self.assertTrue(values <= allowed, values)
            observed[label] = sorted(values)
        print('Existing weighted-cycle coverage annotations:', observed)
