import csv
import json
import os
from pathlib import Path
import random
import subprocess
import tempfile
import unittest
from unittest import mock

from test_flash import load_launcher

CORE = os.environ.get('MEGAHIT_TEST_CORE')


def sam_record(name, sequence, cigar, reference='ref', position=1, flag=0):
    return '{}\t{}\t{}\t{}\t60\t{}\t*\t0\t0\t{}\t*\n'.format(
        name, flag, reference, position, cigar, sequence)


def fasta_sequences(path):
    lines = Path(path).read_text().splitlines()
    return {header[1:].split()[0]: seq for header, seq in zip(lines[::2], lines[1::2])}


@unittest.skipUnless(CORE, 'set MEGAHIT_TEST_CORE to the consensus-enabled core')
class NativeConsensusTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='megahit-consensus-')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.reference = self.root / 'ends.fa'
        self.output = self.root / 'extended.fa'
        self.stats = self.root / 'stats.tsv'

    def run_engine(self, records, references=None, extra_sams=None, check=True):
        references = references or {'ref': 'ACGT'}
        self.reference.write_text(''.join('>{}\n{}\n'.format(name, seq) for name, seq in references.items()))
        header = ''.join('@SQ\tSN:{}\tLN:{}\n'.format(name, len(seq)) for name, seq in references.items())
        sam = self.root / 'reads.sam'
        sam.write_text(header + ''.join(records))
        sams = [str(sam)]
        for i, rows in enumerate(extra_sams or []):
            path = self.root / 'extra.{}.sam'.format(i)
            path.write_text(header + ''.join(rows))
            sams.append(str(path))
        return subprocess.run([CORE, 'endconsensus', str(self.reference), str(self.output),
                               str(self.stats)] + sams, capture_output=True, text=True, check=check)

    def test_five_A_beats_four_G(self):
        rows = [sam_record(str(i), 'ACGTA', '4M1S') for i in range(5)]
        rows += [sam_record('g' + str(i), 'ACGTG', '4M1S') for i in range(4)]
        self.run_engine(rows)
        self.assertEqual(fasta_sequences(self.output), {'ref': 'ACGTA'})
        with self.stats.open() as source:
            stats = list(csv.DictReader(source, delimiter='\t'))[0]
        self.assertEqual(stats['primary_reads'], '9')
        self.assertEqual(stats['right_added'], '1')

    def test_left_and_right_consensus_orientation(self):
        self.run_engine([sam_record('r', 'GAACGTTG', '2S4M2S')])
        self.assertEqual(fasta_sequences(self.output)['ref'], 'GAACGTTG')

    def test_only_outward_side_of_extracted_fragments_is_extended(self):
        references = {'x|start': 'ACGT', 'x|end': 'ACGT', 'whole': 'ACGT'}
        rows = [sam_record(name, 'GAACGTTG', '2S4M2S', reference=name) for name in references]
        self.run_engine(rows, references)
        self.assertEqual(fasta_sequences(self.output),
                         {'x|start': 'GAACGT', 'x|end': 'ACGTTG', 'whole': 'GAACGTTG'})

    def test_soft_clip_projection_from_inset_alignment(self):
        self.run_engine([sam_record('r', 'GAACGTTG', '3S2M3S', position=2)])
        self.assertEqual(fasta_sequences(self.output)['ref'], 'GAACGTTG')

    def test_internal_insertion_does_not_shift_right_projection(self):
        self.run_engine([sam_record('r', 'GAACNGTTG', '2S2M1I2M2S')])
        self.assertEqual(fasta_sequences(self.output)['ref'], 'GAACGTTG')

    def test_deletion_and_reference_skip_preserve_original_sequence(self):
        for operation in ['D', 'N']:
            with self.subTest(operation=operation):
                self.run_engine([sam_record('r', 'GAACTATG', '2S2M1{}2M2S'.format(operation))],
                                {'ref': 'ACGTA'})
                self.assertEqual(fasta_sequences(self.output)['ref'], 'GAACGTATG')

    def test_reverse_flag_seq_is_not_reverse_complemented_again(self):
        self.run_engine([sam_record('r', 'ACGTA', '4M1S', flag=16)])
        self.assertEqual(fasta_sequences(self.output)['ref'], 'ACGTA')

    def test_secondary_supplementary_and_unmapped_do_not_vote(self):
        rows = [sam_record('primary', 'ACGTA', '4M1S')]
        rows += [sam_record(str(flag), 'ACGTG', '4M1S', flag=flag) for flag in [4, 256, 2048]]
        self.run_engine(rows)
        self.assertEqual(fasta_sequences(self.output)['ref'], 'ACGTA')

    def test_tie_stops_extension_at_that_position(self):
        self.run_engine([sam_record('a', 'ACGTAG', '4M2S'), sam_record('b', 'ACGTAT', '4M2S')])
        self.assertEqual(fasta_sequences(self.output)['ref'], 'ACGTA')

    def test_tie_at_boundary_adds_nothing(self):
        self.run_engine([sam_record('a', 'ACGTA', '4M1S'), sam_record('b', 'ACGTG', '4M1S')])
        self.assertEqual(fasta_sequences(self.output), {})
        self.assertEqual(Path(str(self.output) + '.info').read_text(), '0 0\n')

    def test_N_does_not_vote_and_gap_stops_extension(self):
        self.run_engine([sam_record('r', 'GNAACGT', '3S4M')])
        self.assertEqual(fasta_sequences(self.output)['ref'], 'AACGT')

    def test_reference_is_not_polished_from_mismatching_reads(self):
        self.run_engine([sam_record('r', 'TTTTA', '4X1S')])
        self.assertEqual(fasta_sequences(self.output)['ref'], 'ACGTA')

    def test_hard_clips_do_not_invent_bases(self):
        self.run_engine([sam_record('r', 'GAACGTTG', '3H2S4M2S4H')])
        self.assertEqual(fasta_sequences(self.output)['ref'], 'GAACGTTG')

    def test_multiple_sam_libraries_pool_votes(self):
        self.run_engine([sam_record('a', 'ACGTA', '4M1S')],
                        extra_sams=[[sam_record('g1', 'ACGTG', '4M1S'), sam_record('g2', 'ACGTG', '4M1S')]])
        self.assertEqual(fasta_sequences(self.output)['ref'], 'ACGTG')

    def test_header_only_sam_emits_no_unchanged_references(self):
        self.run_engine([], {'one': 'ACGT', 'two': 'NNACGT'})
        self.assertEqual(fasta_sequences(self.output), {})
        self.assertEqual(Path(str(self.output) + '.info').read_text(), '0 0\n')

    def test_only_extended_records_and_exact_info_counts(self):
        self.run_engine([sam_record('r', 'ACGTA', '4M1S')], {'ref': 'ACGT', 'unchanged': 'GGGG'})
        self.assertEqual(fasta_sequences(self.output), {'ref': 'ACGTA'})
        self.assertIn(' flag=0 multi=1 ', self.output.read_text())
        self.assertEqual(Path(str(self.output) + '.info').read_text(), '1 5\n')

    def test_malformed_alignment_is_rejected_before_output(self):
        cases = [('4M0', 'ACGT', 1), ('4M', 'ACGTA', 1), ('0S4M', 'ACGT', 1),
                 ('2M1S2M', 'ACGTA', 1), ('2M', 'AC', 4), ('4M', 'ACGT', 0),
                 ('2S4D2S', 'ACGT', 1)]
        for cigar, seq, pos in cases:
            with self.subTest(cigar=cigar, pos=pos):
                result = self.run_engine([sam_record('bad', seq, cigar, position=pos)], check=False)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(self.output.exists())

    def test_unknown_reference_is_rejected(self):
        result = self.run_engine([sam_record('r', 'ACGT', '4M', reference='wrong')], check=False)
        self.assertNotEqual(result.returncode, 0)


@unittest.skipUnless(CORE, 'set MEGAHIT_TEST_CORE to the consensus-enabled core')
class ConsensusWorkflowTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='megahit-consensus-workflow-')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.launcher = load_launcher()
        self.opt = self.launcher.opt
        self.opt.out_dir = str(self.root)
        self.opt.temp_dir = str(self.root / 'tmp')
        self.opt.merge_pairs = True
        self.launcher.software_info.megahit_core = str(Path(CORE).resolve())
        self.launcher.run_sub_command = lambda cmd, *args, **kwargs: subprocess.run(
            cmd, capture_output=True, text=True, check=True)
        self.reference = Path(self.launcher.contig_ends_path(59))
        self.reference.parent.mkdir(parents=True)
        self.reference.write_text('>ref\nACGT\n')
        self.mapping = self.reference.parent / 'flank_mapping'
        self.mapping.mkdir()
        self.sam = self.mapping / 'library.0.sam'
        self.sam.write_text('@SQ\tSN:ref\tLN:4\n' + sam_record('a', 'ACGTA', '4M1S'))
        (self.mapping / 'mapping.json').write_text(json.dumps({
            'configuration': {'reference': self.launcher.file_signature(str(self.reference))},
            'outputs': [str(self.sam)]}))
        self.output = self.reference.parent / 'k59.contig_ends.extended.fa'

    def test_workflow_and_cleanup(self):
        self.assertEqual(self.launcher.extend_contig_ends(59), str(self.output))
        self.launcher.remove_temporary_files()
        self.assertEqual(fasta_sequences(self.output)['ref'], 'ACGTA')
        self.assertTrue((self.mapping / 'consensus.tsv').is_file())
        self.assertTrue((self.mapping / 'consensus.json').is_file())

    def test_resume_reuses_consensus(self):
        self.launcher.extend_contig_ends(59)
        self.launcher.check_point._logged_checkpoint = 100
        with mock.patch.object(self.launcher, 'run_sub_command', side_effect=AssertionError('Consensus reran')):
            self.launcher.extend_contig_ends(59)

    def test_changed_sam_rebuilds_consensus(self):
        self.launcher.extend_contig_ends(59)
        with self.sam.open('a') as dest:
            dest.write(sam_record('g1', 'ACGTG', '4M1S') + sam_record('g2', 'ACGTG', '4M1S'))
        self.launcher.check_point._logged_checkpoint = 100
        self.launcher.extend_contig_ends(59)
        self.assertEqual(fasta_sequences(self.output)['ref'], 'ACGTG')

    def test_truncated_output_is_rebuilt(self):
        self.launcher.extend_contig_ends(59)
        self.output.write_text('')
        self.launcher.check_point._logged_checkpoint = 100
        self.launcher.extend_contig_ends(59)
        self.assertEqual(fasta_sequences(self.output)['ref'], 'ACGTA')

    def test_empty_consensus_is_a_valid_cached_result(self):
        self.sam.write_text('@SQ\tSN:ref\tLN:4\n')
        self.launcher.extend_contig_ends(59)
        self.assertEqual(self.output.read_text(), '')
        self.assertEqual(Path(str(self.output) + '.info').read_text(), '0 0\n')
        self.launcher.check_point._logged_checkpoint = 100
        with mock.patch.object(self.launcher, 'run_sub_command', side_effect=AssertionError('Consensus reran')):
            self.launcher.extend_contig_ends(59)

    def test_missing_info_is_rebuilt(self):
        self.launcher.extend_contig_ends(59)
        Path(str(self.output) + '.info').unlink()
        self.launcher.check_point._logged_checkpoint = 100
        self.launcher.extend_contig_ends(59)
        self.assertEqual(Path(str(self.output) + '.info').read_text(), '1 5\n')

    def test_local_stage_publishes_consensus_and_info(self):
        self.launcher.check_point.set_file(str(self.root / 'checkpoints.txt'))
        with mock.patch.object(self.launcher, 'extract_contig_ends'), mock.patch.object(self.launcher, 'map_merged_reads'):
            self.launcher.local_assemble(59, 79)
        local = self.root / 'intermediate_contigs/k59.local.fa'
        self.assertEqual(local.read_bytes(), self.output.read_bytes())
        self.assertEqual(Path(str(local) + '.info').read_text(), '1 5\n')

    def test_local_stage_without_merged_reads_clears_legacy_local(self):
        self.opt.merge_pairs = False
        self.launcher.check_point.set_file(str(self.root / 'checkpoints.txt'))
        local = self.root / 'intermediate_contigs/k59.local.fa'
        local.parent.mkdir()
        local.write_text('>legacy\nACGT\n')
        Path(str(local) + '.info').write_text('1 4\n')
        with mock.patch.object(self.launcher, 'run_sub_command', side_effect=AssertionError('Legacy assembler called')):
            self.launcher.local_assemble(59, 79)
        self.assertEqual(local.read_text(), '')
        self.assertEqual(Path(str(local) + '.info').read_text(), '0 0\n')

    def test_legacy_v1_completion_record_does_not_reuse_old_outputs(self):
        self.launcher.extend_contig_ends(59)
        marker = self.mapping / 'consensus.json'
        state = json.loads(marker.read_text())
        state['algorithm'] = 'outward-plurality-v1'
        marker.write_text(json.dumps(state))
        self.launcher.check_point._logged_checkpoint = 100
        with mock.patch.object(self.launcher, 'run_sub_command', wraps=self.launcher.run_sub_command) as run:
            self.launcher.extend_contig_ends(59)
        self.assertEqual(run.call_count, 1)

    def test_graph_builder_accepts_consensus_local_fasta_and_info(self):
        rng = random.Random(71)
        sequence = ''.join(rng.choice('ACGT') for _ in range(300))
        self.reference.write_text('>ref\n' + sequence + '\n')
        self.sam.write_text('@SQ\tSN:ref\tLN:300\n' +
                            sam_record('r', sequence[-50:] + 'A', '50M1S', position=251))
        (self.mapping / 'mapping.json').write_text(json.dumps({
            'configuration': {'reference': self.launcher.file_signature(str(self.reference))},
            'outputs': [str(self.sam)]}))
        self.opt.memory = 100000000
        self.opt.num_cpu_threads = 2
        self.opt.keep_tmp_files = True
        self.launcher.check_point.set_file(str(self.root / 'checkpoints.txt'))
        with mock.patch.object(self.launcher, 'extract_contig_ends'), mock.patch.object(self.launcher, 'map_merged_reads'):
            self.launcher.local_assemble(59, 79)
        original = self.root / 'intermediate_contigs/k59.contigs.fa'
        original.write_text('>original flag=0 multi=1\n' + sequence + '\n')
        Path(str(original) + '.info').write_text('1 300\n')
        bubble = original.with_name('k59.bubble_seq.fa')
        bubble.write_text('')
        Path(str(bubble) + '.info').write_text('0 0\n')
        self.launcher.build_graph(79, 59)
        self.assertEqual((self.root / 'intermediate_contigs/k59.local.fa.info').read_text(), '1 301\n')
        self.assertTrue((self.root / 'tmp/k79/79.sdbg_info').is_file())
        assembled = self.root / 'assembled'
        subprocess.run([CORE, 'assemble', '-s', str(self.root / 'tmp/k79/79'), '-o', str(assembled),
                        '-t', '2', '--prune_level', '0', '--bubble_level', '0', '--merge_len', '0',
                        '--max_tip_len', '0', '--min_standalone', '0'],
                       capture_output=True, text=True, check=True)
        sequences = list(fasta_sequences(str(assembled) + '.contigs.fa').values())
        expected = sequence + 'A'
        reverse = expected.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
        self.assertTrue(expected in sequences or reverse in sequences)

    def test_no_local_excludes_existing_local_file_from_graph_command(self):
        self.opt.no_local = True
        self.opt.memory = 100000000
        self.opt.num_cpu_threads = 2
        self.launcher.check_point.set_file(str(self.root / 'checkpoints.txt'))
        local = self.root / 'intermediate_contigs/k59.local.fa'
        local.parent.mkdir()
        local.write_text('>old flag=0 multi=1\nACGT\n')
        Path(str(local) + '.info').write_text('1 4\n')
        graph_dir = self.root / 'tmp/k79'
        graph_dir.mkdir()
        (graph_dir / '79.edges.0').write_bytes(b'edge')
        with mock.patch.object(self.launcher, 'run_sub_command') as run:
            self.launcher.build_graph(79, 59)
        self.assertNotIn('--local_contig', run.call_args[0][0])

    def test_stale_reference_is_rejected(self):
        self.reference.write_text('>ref\nACGTA\n')
        with self.assertRaises(self.launcher.Usage):
            self.launcher.extend_contig_ends(59)

    def test_bad_sam_does_not_publish_consensus(self):
        self.sam.write_text(sam_record('bad', 'ACGT', '3M'))
        with self.assertRaises(subprocess.CalledProcessError):
            self.launcher.extend_contig_ends(59)
        self.assertFalse(self.output.exists())
        self.assertFalse((self.mapping / 'consensus.json').exists())
        self.assertFalse(Path(str(self.output) + '.tmp').exists())


if __name__ == '__main__':
    unittest.main()
