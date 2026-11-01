import random
from pathlib import Path
import tempfile
import unittest

from test_flash import load_launcher


class ContigEndExtractionTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='megahit-contig-ends-')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.launcher = load_launcher()
        self.opt = self.launcher.opt
        self.opt.out_dir = str(self.root)
        self.opt.temp_dir = str(self.root / 'tmp')
        self.source = self.root / 'intermediate_contigs/k59.contigs.fa'
        self.source.parent.mkdir()
        self.output = self.root / 'tmp/k59/k59.contig_ends.fa'

    def sequence(self, length):
        rng = random.Random(length)
        return ''.join(rng.choice('ACGT') for _ in range(length))

    def records(self):
        lines = self.output.read_text().splitlines()
        return list(zip(lines[::2], lines[1::2]))

    def test_threshold_and_exact_end_sequences(self):
        for length in [1, 289, 290, 579, 580, 581, 1000]:
            with self.subTest(length=length):
                seq = self.sequence(length)
                header = '>k59_42 flag=0 multi=12 len={}'.format(length)
                self.source.write_text(header + '\n' + seq + '\n')
                output = self.launcher.extract_contig_ends(59)
                self.assertEqual(output, str(self.output))
                if length < 580:
                    self.assertEqual(self.records(), [(header, seq)])
                else:
                    self.assertEqual(self.records(), [
                        ('>k59_42|start source_length={} region=1-290'.format(length), seq[:290]),
                        ('>k59_42|end source_length={} region={}-{}'.format(length, length - 289, length), seq[-290:])])
                self.assertFalse(Path(output + '.tmp').exists())

    def test_multiple_multiline_records_and_no_final_newline(self):
        long_seq = self.sequence(1000)
        self.source.write_bytes(('\r\n>short description\r\nAC\r\nGT\r\n\r\n'
                                 '>long description\r\n' + long_seq[:500] + '\r\n' + long_seq[500:]).encode())
        self.launcher.extract_contig_ends(59)
        self.assertEqual(self.records(), [
            ('>short description', 'ACGT'),
            ('>long|start source_length=1000 region=1-290', long_seq[:290]),
            ('>long|end source_length=1000 region=711-1000', long_seq[-290:])])

    def test_empty_input_creates_empty_fasta(self):
        self.source.write_text('')
        self.launcher.extract_contig_ends(59)
        self.assertEqual(self.output.read_text(), '')

    def test_identical_export_preserves_timestamp_for_mapping_resume(self):
        self.source.write_text('>short\nACGT\n')
        self.launcher.extract_contig_ends(59)
        before = self.output.stat().st_mtime_ns
        self.launcher.extract_contig_ends(59)
        self.assertEqual(self.output.stat().st_mtime_ns, before)

    def test_missing_input_reports_error(self):
        with self.assertRaises(self.launcher.Usage):
            self.launcher.extract_contig_ends(59)
        self.assertFalse(self.output.exists())

    def test_invalid_fasta_is_not_published(self):
        for content in ['ACGT\n', '>\nACGT\n', '>   \nACGT\n']:
            with self.subTest(content=content):
                self.source.write_text(content)
                with self.assertRaises(self.launcher.Usage):
                    self.launcher.extract_contig_ends(59)
                self.assertFalse(self.output.exists())
                self.assertFalse(Path(str(self.output) + '.tmp').exists())

    def test_failure_preserves_previous_export(self):
        self.source.write_text('>short\nACGT\n')
        self.launcher.extract_contig_ends(59)
        before = self.output.read_bytes()
        self.source.write_text('>\nACGT\n')
        with self.assertRaises(self.launcher.Usage):
            self.launcher.extract_contig_ends(59)
        self.assertEqual(self.output.read_bytes(), before)

    def test_nonpositive_flank_length_is_rejected(self):
        for length in [0, -1]:
            with self.subTest(length=length), self.assertRaises(self.launcher.Usage):
                self.launcher.extract_contig_ends(59, flank_length=length)

    def test_export_uses_output_tmp_with_custom_scratch(self):
        self.opt.temp_dir = str(self.root / 'scratch')
        Path(self.opt.temp_dir).mkdir()
        self.source.write_text('>short\nACGT\n')
        self.launcher.extract_contig_ends(59)
        self.launcher.remove_temporary_files()
        self.assertEqual(self.records(), [('>short', 'ACGT')])
        self.assertFalse(Path(self.opt.temp_dir).exists())

    def test_cleanup_retains_export_without_merge_pairs(self):
        self.assertFalse(self.opt.merge_pairs)
        self.source.write_text('>short\nACGT\n')
        self.launcher.extract_contig_ends(59)
        junk = self.output.parent / '79.edges.0'
        junk.write_bytes(b'temporary')
        old_graph = self.root / 'tmp/k39'
        old_graph.mkdir()
        (old_graph / '59.sdbg.0').write_bytes(b'temporary')
        self.launcher.remove_temporary_files()
        self.assertEqual(self.records(), [('>short', 'ACGT')])
        self.assertFalse(junk.exists())
        self.assertFalse(old_graph.exists())

    def test_cleanup_without_exports_keeps_original_behavior(self):
        Path(self.opt.temp_dir).mkdir()
        (Path(self.opt.temp_dir) / 'discard.bin').write_bytes(b'temporary')
        self.launcher.remove_temporary_files()
        self.assertFalse(Path(self.opt.temp_dir).exists())

    def test_custom_reduction_and_next_k_paths(self):
        source = self.source.with_name('k39.contigs.fa')
        source.write_text('>k39_1\nACGT\n')
        output = self.launcher.extract_contig_ends(39)
        self.assertEqual(output, str(self.root / 'tmp/k39/k39.contig_ends.fa'))
        self.assertEqual(Path(output).read_text(), '>k39_1\nACGT\n')


if __name__ == '__main__':
    unittest.main()
