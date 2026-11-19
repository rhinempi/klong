"""Reduction contig output thresholds, independently of graph pruning."""
import collections
import os
from pathlib import Path
import random
import subprocess
import tempfile
import unittest

from test_assembly_optimizations import records, semantic_record
from test_query_merge import canonical, pack

CORE = os.environ.get('MEGAHIT_TEST_CORE')


@unittest.skipUnless(CORE, 'set MEGAHIT_TEST_CORE to the built core')
class ReductionOutputTest(unittest.TestCase):
    def test_output_cutoff_excludes_next_k_and_includes_next_k_plus_one(self):
        with tempfile.TemporaryDirectory(prefix='reduction-output-') as directory:
            root = Path(directory)
            rng = random.Random(15523)
            edges = {}
            lengths = [22, 29, 30, 31, 59, 60, 80]
            for length in lengths:
                sequence = ''.join(rng.choice('ACGT') for _ in range(length))
                for i in range(length - 21):
                    edges[canonical(sequence[i:i + 22])] = 7
            prefix, graph = root / 'edges', root / 'graph'
            Path(str(prefix) + '.edges.0').write_bytes(b''.join(pack(seq, count) for seq, count in edges.items()))
            Path(str(prefix) + '.edges.info').write_text(
                'kmer_size 21\nwords_per_edge 2\nnum_files 1\nnum_buckets 0\n'
                'num_edges {}\nis_sorted 0\n'.format(len(edges)))
            subprocess.run([CORE, 'seq2sdbg', '--input_prefix', str(prefix), '-k', '21',
                            '--output_prefix', str(graph), '--host_mem', '200000000',
                            '--mem_flag', '1', '--num_cpu_threads', '2'],
                           capture_output=True, check=True, timeout=60)
            unfiltered = None
            for cutoff in (0, 30, 60, 81):
                output = root / ('out-' + str(cutoff))
                subprocess.run([CORE, 'assemble', '-s', str(graph), '-o', str(output), '-t', '2',
                                '--prune_level', '0', '--bubble_level', '0', '--merge_len', '0',
                                '--max_tip_len', '0', '--disconnect_ratio', '0',
                                '--min_contig_len', str(cutoff)],
                               capture_output=True, check=True, timeout=60)
                rows = list(records(Path(str(output) + '.contigs.fa')))
                if cutoff == 0:
                    unfiltered = rows
                    self.assertEqual(sorted(len(sequence) for _, sequence in rows), lengths)
                self.assertEqual(collections.Counter(semantic_record(*row) for row in rows),
                                 collections.Counter(semantic_record(*row) for row in unfiltered if len(row[1]) >= cutoff))
                info = Path(str(output) + '.contigs.fa.info').read_text().split()
                self.assertEqual(int(info[0]), len(rows))

    def test_negative_output_cutoff_is_rejected(self):
        result = subprocess.run([CORE, 'assemble', '-s', 'unused', '--min_contig_len', '-1'],
                                capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('cannot be negative', result.stderr)
