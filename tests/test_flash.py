"""Run with MEGAHIT_TEST_CORE=/path/to/megahit_core python3 -m unittest discover -s tests."""

import bz2
import gzip
import importlib.machinery
import importlib.util
import io
import os
from pathlib import Path
import random
import shutil
import struct
import subprocess
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
CORE = os.environ.get('MEGAHIT_TEST_CORE')


def load_launcher():
    loader = importlib.machinery.SourceFileLoader('megahit_launcher', str(ROOT / 'src/megahit'))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def fastq(records):
    return ''.join('@{}\n{}\n+\n{}\n'.format(name, seq, 'I' * len(seq))
                   for name, seq in records)


def sequences(path):
    with gzip.open(path, 'rt') as source:
        return source.read().splitlines()[1::2]


def binary_sequences(path):
    result = []
    with open(path, 'rb') as source:
        while True:
            header = source.read(4)
            if not header:
                return result
            length, = struct.unpack('=I', header)
            count = (length + 15) // 16
            words = struct.unpack('=' + 'I' * count, source.read(4 * count))
            result.append(''.join('ACGT'[(words[i // 16] >> (30 - 2 * (i % 16))) & 3]
                                  for i in range(length)))


@unittest.skipUnless(CORE, 'set MEGAHIT_TEST_CORE to the newly built core executable')
class FlashIntegrationTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='megahit-flash-test-')
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        self.launcher = load_launcher()
        self.opt = self.launcher.opt
        self.opt.out_dir = str(self.directory)
        self.opt.temp_dir = str(self.directory / 'tmp')
        Path(self.opt.temp_dir).mkdir()
        self.opt.merge_pairs = True
        self.opt.num_cpu_threads = 1
        self.opt.flash_max_mismatch_density = 0
        self.launcher.software_info.megahit_core = str(Path(CORE).resolve())
        self.launcher.check_point.set_file(str(self.directory / 'checkpoints.txt'))
        self.launcher.run_sub_command = lambda cmd, *args, **kwargs: subprocess.run(
            cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        rng = random.Random(42)
        self.fragment = ''.join(rng.choice('ACGT') for _ in range(100))
        complement = str.maketrans('ACGT', 'TGCA')
        self.left = self.fragment[:70]
        self.right = self.fragment[-70:].translate(complement)[::-1]
        self.single = self.directory / 'single.fastq'
        self.single.write_text(fastq([('single', 'GATTACA')]))
        self.opt.se = [str(self.single)]

    def paired_inputs(self, left=None, right=None):
        f1, f2 = self.directory / 'r1.fastq', self.directory / 'r2.fastq'
        f1.write_text(fastq(left if left is not None else [('merge/1', self.left), ('keep/1', 'A' * 70)]))
        f2.write_text(fastq(right if right is not None else [('merge/2', self.right), ('keep/2', 'C' * 70)]))
        self.opt.pe1, self.opt.pe2 = [str(f1)], [str(f2)]
        return f1, f2

    def merged_path(self, i=0):
        return self.directory / 'tmp/flash' / ('library.{}.merged.fa.gz'.format(i))

    def test_binary_library_is_byte_identical_with_merging(self):
        self.paired_inputs()
        self.launcher.create_library_file()
        self.launcher.build_library()
        binary = Path(self.opt.read_lib_path + '.bin').read_bytes()
        metadata = Path(self.opt.read_lib_path + '.lib_info').read_bytes()
        before = (self.opt.pe1[:], self.opt.pe2[:], self.opt.pe12[:], self.opt.se[:])
        self.launcher.merge_paired_reads()
        self.assertEqual(sequences(self.merged_path()), [self.fragment])
        self.assertEqual((self.opt.pe1, self.opt.pe2, self.opt.pe12, self.opt.se), before)
        self.launcher.create_library_file()
        self.launcher.build_library()
        self.assertEqual(Path(self.opt.read_lib_path + '.bin').read_bytes(), binary)
        self.assertEqual(Path(self.opt.read_lib_path + '.lib_info').read_bytes(), metadata)
        self.assertEqual(binary_sequences(self.opt.read_lib_path + '.bin'),
                         [self.left, self.right, 'A' * 70, 'C' * 70, 'GATTACA'])
        info = Path(self.opt.read_lib_path + '.lib_info').read_text().splitlines()
        self.assertEqual(info[0].split(), ['287', '5'])
        self.assertEqual([line.split()[-1] for line in info[2::2]], ['1', '0'])

    def test_compressed_interleaved_input(self):
        content = fastq([('merge/1', self.left), ('merge/2', self.right),
                         ('keep/1', 'A' * 70), ('keep/2', 'C' * 70)])
        for suffix, opener in [('.gz', gzip.open), ('.bz2', bz2.open)]:
            path = self.directory / ('interleaved.fastq' + suffix)
            with opener(path, 'wt') as dest:
                dest.write(content)
            self.opt.pe1, self.opt.pe2 = [], []
            self.opt.pe12, self.opt.se = [str(path)], [str(self.single)]
            self.launcher.merge_paired_reads()
            self.assertEqual(sequences(self.merged_path()), [self.fragment])
            self.assertEqual(self.opt.pe12, [str(path)])
            self.assertEqual(self.opt.se, [str(self.single)])
            self.assertEqual(list(Path(self.opt.temp_dir).glob('flash-*')), [])

    def test_all_merged_and_none_merged(self):
        for left, right, expected in [(self.left, self.right, [self.fragment]),
                                      ('A' * 70, 'C' * 70, [])]:
            self.paired_inputs([('read/1', left)], [('read/2', right)])
            self.opt.se = []
            self.launcher.merge_paired_reads()
            self.assertEqual(len(self.opt.pe1), 1)
            self.assertEqual(self.opt.se, [])
            self.assertEqual(sequences(self.merged_path()), expected)
            self.launcher.create_library_file()
            self.launcher.build_library()
            self.assertEqual(binary_sequences(self.opt.read_lib_path + '.bin'),
                             [left, right])

    def test_resume_reuses_outputs_without_duplicate_reads(self):
        f1, f2 = self.paired_inputs()
        self.launcher.merge_paired_reads()
        self.launcher.create_library_file()
        expected = (self.opt.pe1[:], self.opt.pe2[:], self.opt.se[:])
        self.opt.pe1, self.opt.pe2, self.opt.se = [str(f1)], [str(f2)], [str(self.single)]
        self.launcher.check_point._current_checkpoint = 0
        self.launcher.check_point.load_for_continue()
        with mock.patch.object(self.launcher, 'run_sub_command', side_effect=AssertionError('FLASH reran')):
            self.launcher.merge_paired_reads()
        self.assertEqual((self.opt.pe1, self.opt.pe2, self.opt.se), expected)
        self.launcher.create_library_file()
        self.launcher.build_library()
        self.assertEqual(binary_sequences(self.opt.read_lib_path + '.bin'),
                         [self.left, self.right, 'A' * 70, 'C' * 70, 'GATTACA'])

    def test_resume_with_missing_merge_product_rebuilds(self):
        self.paired_inputs()
        self.launcher.merge_paired_reads()
        self.merged_path().unlink()
        self.launcher.check_point._logged_checkpoint = 0
        self.launcher.merge_paired_reads()
        self.assertEqual(sequences(self.merged_path()), [self.fragment])

    def test_resume_after_build_still_exports_merged_reads(self):
        self.paired_inputs()
        self.launcher.check_point._logged_checkpoint = 1
        self.launcher.merge_paired_reads()
        self.assertEqual(sequences(self.merged_path()), [self.fragment])

    def test_merging_disabled_leaves_inputs_unchanged(self):
        self.paired_inputs()
        self.opt.merge_pairs = False
        before = (self.opt.pe1[:], self.opt.pe2[:], self.opt.se[:])
        self.launcher.merge_paired_reads()
        self.assertEqual((self.opt.pe1, self.opt.pe2, self.opt.se), before)
        self.assertFalse((self.directory / 'tmp/flash').exists())

    def test_bad_fastq_has_no_completion_manifest(self):
        f1, _ = self.paired_inputs()
        f1.write_text('>not-fastq\nACGT\n')
        with self.assertRaises(subprocess.CalledProcessError):
            self.launcher.merge_paired_reads()
        self.assertFalse((self.directory / 'tmp/flash/merged.json').exists())

    def test_failed_fresh_merge_invalidates_previous_manifest(self):
        f1, f2 = self.paired_inputs()
        self.launcher.merge_paired_reads()
        self.opt.pe1, self.opt.pe2, self.opt.se = [str(f1)], [str(f2)], []
        f1.write_text('>not-fastq\nACGT\n')
        with self.assertRaises(subprocess.CalledProcessError):
            self.launcher.merge_paired_reads()
        self.assertFalse((self.directory / 'tmp/flash/merged.json').exists())

    def test_unequal_mate_counts_fail(self):
        self.paired_inputs([('read/1', self.left)], [])
        with self.assertRaises(subprocess.CalledProcessError):
            self.launcher.merge_paired_reads()
        self.assertFalse((self.directory / 'tmp/flash/merged.json').exists())

    def test_multiple_libraries(self):
        f1, f2 = self.paired_inputs()
        self.opt.pe1, self.opt.pe2 = [str(f1)] * 2, [str(f2)] * 2
        self.launcher.merge_paired_reads()
        self.assertEqual(len(self.opt.pe1), 2)
        self.assertEqual(self.opt.se, [str(self.single)])
        self.assertEqual([sequences(self.merged_path(i)) for i in range(2)], [[self.fragment]] * 2)

    def test_workflow_forwards_compression_thread_count(self):
        self.paired_inputs()
        self.opt.num_cpu_threads = 2
        with mock.patch.object(self.launcher, 'fastq_to_fasta_gz',
                               wraps=self.launcher.fastq_to_fasta_gz) as convert:
            self.launcher.merge_paired_reads()
        self.assertEqual(convert.call_args.kwargs['num_threads'], 2)
        self.assertEqual(sequences(self.merged_path()), [self.fragment])

    def test_compression_failure_does_not_publish_export_or_manifest(self):
        self.paired_inputs()
        def fail(source, output, **kwargs):
            Path(output).write_bytes(b'partial gzip')
            raise self.launcher.Usage('compression failed')
        with mock.patch.object(self.launcher, 'fastq_to_fasta_gz', side_effect=fail):
            with self.assertRaises(self.launcher.Usage):
                self.launcher.merge_paired_reads()
        self.assertFalse(self.merged_path().exists())
        self.assertFalse(Path(str(self.merged_path()) + '.tmp').exists())
        self.assertFalse((self.directory / 'tmp/flash/merged.json').exists())

    def test_cleanup_preserves_merged_results(self):
        self.paired_inputs()
        self.launcher.merge_paired_reads()
        junk = Path(self.opt.temp_dir) / 'discard.bin'
        junk.write_bytes(b'temporary')
        self.launcher.remove_temporary_files()
        self.assertFalse(junk.exists())
        self.assertEqual(sequences(self.merged_path()), [self.fragment])

    def test_custom_tmp_dir_keeps_output_under_output_directory(self):
        self.paired_inputs()
        custom = self.directory / 'custom-tmp'
        custom.mkdir()
        self.opt.temp_dir = str(custom)
        self.launcher.merge_paired_reads()
        self.launcher.remove_temporary_files()
        self.assertFalse(custom.exists())
        self.assertEqual(sequences(self.merged_path()), [self.fragment])


class FastaConversionTest(unittest.TestCase):
    def setUp(self):
        self.launcher = load_launcher()
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.source = Path(self.temp.name) / 'merged.fastq'
        self.output = Path(self.temp.name) / 'merged.fa.gz'
        self.source.write_bytes(b'@read description\nACGT\n+\nIIII\n')

    def test_missing_pigz_falls_back_to_python_gzip(self):
        with mock.patch.object(self.launcher, 'find_bundled_pigz', return_value=None), \
                mock.patch.object(self.launcher.subprocess, 'Popen') as popen:
            self.launcher.fastq_to_fasta_gz(str(self.source), str(self.output), num_threads=4)
        popen.assert_not_called()
        self.assertEqual(gzip.decompress(self.output.read_bytes()), b'>read description\nACGT\n')

    @unittest.skipUnless(shutil.which('pigz'), 'pigz is not installed')
    def test_pigz_preserves_content_and_uses_requested_threads(self):
        self.source.write_bytes(b'@read description\r\nACGT\r\n+\r\nIIII\r\n'
                                b'@second\nGATTACA\n+\nIIIIIII\n')
        with mock.patch.object(self.launcher, 'find_bundled_pigz', return_value=shutil.which('pigz')), \
                mock.patch.object(self.launcher.subprocess, 'Popen',
                               wraps=subprocess.Popen) as popen:
            self.launcher.fastq_to_fasta_gz(str(self.source), str(self.output), num_threads=4)
        command = popen.call_args.args[0]
        self.assertEqual(command[command.index('-p') + 1], '4')
        self.assertEqual(gzip.decompress(self.output.read_bytes()),
                         b'>read description\nACGT\n>second\nGATTACA\n')
        self.assertTrue(popen.call_args.kwargs['stdin'] is subprocess.PIPE)

    @unittest.skipUnless(shutil.which('pigz'), 'pigz is not installed')
    def test_pigz_handles_empty_and_rejects_truncated_input(self):
        self.source.write_bytes(b'')
        with mock.patch.object(self.launcher, 'find_bundled_pigz', return_value=shutil.which('pigz')):
            self.launcher.fastq_to_fasta_gz(str(self.source), str(self.output), num_threads=2)
            self.assertEqual(gzip.decompress(self.output.read_bytes()), b'')
            self.source.write_bytes(b'@read\nACGT\n+\n')
            with self.assertRaises(self.launcher.Usage):
                self.launcher.fastq_to_fasta_gz(str(self.source), str(self.output), num_threads=2)
        self.assertFalse(self.output.exists())

    def fake_compressor(self, stream=None, return_code=2):
        process = mock.Mock()
        process.stdin = io.BytesIO() if stream is None else stream
        process.wait.return_value = return_code
        process.poll.return_value = return_code

        def launch(*args, **kwargs):
            kwargs['stderr'].write(b'compression error')
            return process

        self.launcher.find_bundled_pigz = mock.Mock(return_value='/fake/pigz')
        self.launcher.subprocess = mock.Mock(Popen=mock.Mock(side_effect=launch), PIPE=subprocess.PIPE)
        return process

    def test_compressor_failure_is_reported_and_output_removed(self):
        process = self.fake_compressor()
        with self.assertRaises(self.launcher.Usage) as failure:
            self.launcher.fastq_to_fasta_gz(str(self.source), str(self.output), num_threads=2)
        self.assertIn('compression error', failure.exception.msg)
        self.assertFalse(self.output.exists())
        self.assertTrue(process.stdin.closed)
        process.wait.assert_called()

    def test_compressor_startup_failure_removes_partial_output(self):
        self.launcher.find_bundled_pigz = mock.Mock(return_value='/missing/pigz')
        self.launcher.subprocess = mock.Mock(Popen=mock.Mock(side_effect=OSError('cannot execute')),
                                            PIPE=subprocess.PIPE)
        with self.assertRaises(self.launcher.Usage) as failure:
            self.launcher.fastq_to_fasta_gz(str(self.source), str(self.output), num_threads=2)
        self.assertIn('Cannot start pigz', failure.exception.msg)
        self.assertFalse(self.output.exists())

    def test_broken_pipe_is_reported_and_child_reaped(self):
        class BrokenStream(io.BytesIO):
            def write(self, data):
                raise BrokenPipeError('compressor stopped')
        process = self.fake_compressor(stream=BrokenStream())
        with self.assertRaises(self.launcher.Usage):
            self.launcher.fastq_to_fasta_gz(str(self.source), str(self.output), num_threads=2)
        self.assertFalse(self.output.exists())
        process.wait.assert_called()

    def test_invalid_input_kills_and_reaps_compressor(self):
        process = self.fake_compressor(return_code=None)
        self.source.write_bytes(b'@read\nACGT\n+\nIII\n')
        with self.assertRaises(self.launcher.Usage):
            self.launcher.fastq_to_fasta_gz(str(self.source), str(self.output), num_threads=2)
        process.kill.assert_called_once()
        process.wait.assert_called()
        self.assertFalse(self.output.exists())

    def test_only_bundled_executable_is_selected(self):
        self.launcher.__file__ = str(Path(self.temp.name) / 'megahit')
        pigz = Path(self.temp.name) / 'pigz'
        self.assertIsNone(self.launcher.find_bundled_pigz())
        pigz.write_bytes(b'placeholder')
        pigz.chmod(0o644)
        self.assertIsNone(self.launcher.find_bundled_pigz())
        pigz.chmod(0o755)
        self.assertEqual(self.launcher.find_bundled_pigz(), str(pigz))

    def test_missing_bundle_ignores_system_pigz(self):
        self.launcher.__file__ = str(Path(self.temp.name) / 'megahit')
        with mock.patch.object(self.launcher.shutil, 'which', return_value='/usr/bin/pigz') as which, \
                mock.patch.object(self.launcher.subprocess, 'Popen') as popen:
            self.launcher.fastq_to_fasta_gz(str(self.source), str(self.output), num_threads=4)
        which.assert_not_called()
        popen.assert_not_called()
        self.assertEqual(gzip.decompress(self.output.read_bytes()), b'>read description\nACGT\n')

    def test_fast_compression_preserves_headers_and_sequences(self):
        launcher = load_launcher()
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'merged.fastq'
            output = Path(directory) / 'merged.fa.gz'
            source.write_bytes(b'@read1 description\r\nACGT\r\n+\r\nIIII\r\n'
                               b'@read2\nGATTACA\n+\nIIIIIII\n')
            with mock.patch.object(launcher.gzip, 'open', wraps=gzip.open) as opener:
                launcher.fastq_to_fasta_gz(str(source), str(output))
                self.assertEqual(opener.call_args.kwargs['compresslevel'], 1)
            self.assertEqual(gzip.decompress(output.read_bytes()),
                             b'>read1 description\nACGT\n>read2\nGATTACA\n')

    def test_empty_input_is_a_valid_empty_gzip(self):
        launcher = load_launcher()
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'empty.fastq'
            output = Path(directory) / 'empty.fa.gz'
            source.write_bytes(b'')
            launcher.fastq_to_fasta_gz(str(source), str(output))
            self.assertEqual(gzip.decompress(output.read_bytes()), b'')

    def test_invalid_fastq_is_rejected(self):
        launcher = load_launcher()
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'invalid.fastq'
            output = Path(directory) / 'invalid.fa.gz'
            for data in [b'read\nACGT\n+\nIIII\n', b'@read\nACGT\n+\nIII\n']:
                source.write_bytes(data)
                with self.assertRaises(launcher.Usage):
                    launcher.fastq_to_fasta_gz(str(source), str(output))


class FlashOptionsTest(unittest.TestCase):
    def test_options_parse(self):
        launcher = load_launcher()
        launcher.parse_option(['--merge-pairs', '--flash-min-overlap', '20',
                               '--flash-max-overlap', '100',
                               '--flash-max-mismatch-density', '0.1',
                               '--flash-phred-offset', '64'])
        self.assertTrue(launcher.opt.merge_pairs)
        self.assertEqual(launcher.opt.flash_min_overlap, 20)
        self.assertEqual(launcher.opt.flash_max_overlap, 100)
        self.assertEqual(launcher.opt.flash_max_mismatch_density, 0.1)
        self.assertEqual(launcher.opt.flash_phred_offset, 64)

    def test_invalid_parameters(self):
        for key, value in [('flash_min_overlap', 0), ('flash_max_overlap', 2),
                           ('flash_max_mismatch_density', float('nan')),
                           ('flash_max_mismatch_density', 1.1), ('flash_phred_offset', 42)]:
            launcher = load_launcher()
            setattr(launcher.opt, key, value)
            with self.subTest(key=key, value=value), self.assertRaises(launcher.Usage):
                launcher.check_and_correct_option()


if __name__ == '__main__':
    unittest.main()
