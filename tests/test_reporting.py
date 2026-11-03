"""Binary output must not depend on optional edge diagnostics."""
import os
from pathlib import Path
import struct
import subprocess
import tempfile
import unittest

from test_flash import load_launcher

CORE = os.environ.get('MEGAHIT_TEST_CORE')


class ReportingDefaultsTest(unittest.TestCase):
    def test_launcher_disables_reports_by_default(self):
        self.assertFalse(load_launcher().opt.human_readable_output)


@unittest.skipUnless(CORE, 'set MEGAHIT_TEST_CORE to the newly built core')
class BinaryOnlyReductionTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='megahit-reporting-')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.short = self.write_edges('short', 3, ['AACC', 'ACGT', 'TTTT'])
        self.long = self.write_edges('long', 5, ['AACCGG', 'GGTTAA'])

    def write_edges(self, name, k, sequences):
        prefix = self.root / name
        records = []
        for sequence in sequences:
            packed = 7
            for i, base in enumerate(sequence):
                packed |= 'ACGT'.index(base) << (30 - i * 2)
            records.append(struct.pack('=I', packed))
        Path(str(prefix) + '.edges.0').write_bytes(b''.join(records))
        Path(str(prefix) + '.edges.info').write_text(
            'kmer_size {}\nwords_per_edge 1\nnum_files 1\nnum_buckets 0\n'
            'num_edges {}\nis_sorted 0\n'.format(k, len(sequences)))
        return prefix

    def run_prefix(self, mode, report=None):
        output = self.root / '{}-{}'.format(mode, report)
        command = [CORE, 'prefixbf', '--mode', mode, '--short_k', '3', '--long_k', '5',
                   '--long_prefix', str(self.long), '--output_prefix', str(output),
                   '--num_cpu_threads', '1', '--max_bloom_gb', '0.01']
        if mode == 'reduce':
            command += ['--short_prefix', str(self.short)]
        else:
            contigs = self.root / 'contigs.fa'
            contigs.write_text('>contig\nCAATTT\n')
            command += ['--contig', str(contigs)]
        if report is not None:
            command += ['--human_readable', str(report)]
        subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        return output

    def test_reduction_binary_identical_with_and_without_reports(self):
        default = self.run_prefix('reduce')
        disabled = self.run_prefix('reduce', 0)
        enabled = self.run_prefix('reduce', 1)
        for suffix in ['.edges.0', '.edges.info']:
            self.assertEqual(Path(str(default) + suffix).read_bytes(),
                             Path(str(enabled) + suffix).read_bytes())
            self.assertEqual(Path(str(disabled) + suffix).read_bytes(),
                             Path(str(enabled) + suffix).read_bytes())
        self.assertEqual(Path(str(default) + '.edges.0').stat().st_size, 8)
        self.assertEqual(list(self.root.glob(default.name + '*.txt')), [])
        self.assertEqual(list(self.root.glob(disabled.name + '*.txt')), [])
        self.assertIn('AACC', Path(str(enabled) + '.removed.txt').read_text())

    def test_merge_binary_identical_with_and_without_reports(self):
        default = self.run_prefix('merge')
        enabled = self.run_prefix('merge', 1)
        for suffix in ['.edges.0', '.edges.info']:
            self.assertEqual(Path(str(default) + suffix).read_bytes(),
                             Path(str(enabled) + suffix).read_bytes())
        self.assertEqual(list(self.root.glob(default.name + '*.txt')), [])
        self.assertTrue(Path(str(enabled) + '.contig_edges.txt').is_file())

    def test_count_writer_is_binary_only_unless_explicitly_enabled(self):
        reads = self.root / 'reads.fa'
        reads.write_text('>read\n' + 'AACGTTGCAACCTGATCGTA' * 10 + '\n')
        library = self.root / 'reads.lib'
        library.write_text(str(reads) + '\nse ' + str(reads) + '\n')
        subprocess.run([CORE, 'buildlib', str(library), str(library), '1'],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        outputs = []
        for report in [None, '1']:
            output = self.root / 'count-{}'.format(report)
            env = os.environ.copy()
            env.pop('MEGAHIT_HUMAN_READABLE_OUTPUT', None)
            if report is not None:
                env['MEGAHIT_HUMAN_READABLE_OUTPUT'] = report
            subprocess.run([CORE, 'count', '-k', '21', '-m', '1', '--host_mem', '100000000',
                            '--mem_flag', '0', '--output_prefix', str(output),
                            '--num_cpu_threads', '1', '--read_lib_file', str(library)],
                           env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            outputs.append(output)
        self.assertEqual(list(self.root.glob(outputs[0].name + '*.txt')), [])
        self.assertTrue(list(self.root.glob(outputs[1].name + '*.txt')))
        for suffix in ['.edges.0', '.edges.info']:
            self.assertEqual(Path(str(outputs[0]) + suffix).read_bytes(),
                             Path(str(outputs[1]) + suffix).read_bytes())
