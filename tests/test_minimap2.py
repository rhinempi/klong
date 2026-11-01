"""Run with MEGAHIT_TEST_CORE set to the FLASH/Minimap2-enabled core."""
import gzip
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


@unittest.skipUnless(CORE, 'set MEGAHIT_TEST_CORE to the newly built core')
class Minimap2IntegrationTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='megahit-minimap2-test-')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.launcher = load_launcher()
        self.opt = self.launcher.opt
        self.opt.out_dir = str(self.root)
        self.opt.temp_dir = str(self.root / 'tmp')
        self.opt.merge_pairs = True
        self.opt.num_cpu_threads = 2
        self.launcher.software_info.megahit_core = str(Path(CORE).resolve())
        self.launcher.run_sub_command = lambda cmd, *args, **kwargs: subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        self.reference = self.root / 'tmp/k59/k59.contig_ends.fa'
        self.reference.parent.mkdir(parents=True)
        rng = random.Random(103)
        self.sequence = ''.join(rng.choice('ACGT') for _ in range(2000))
        self.reference.write_text('>reference\n' + self.sequence + '\n')
        self.flash = self.root / 'tmp/flash'
        self.flash.mkdir(parents=True)
        self.query = self.flash / 'library.0.merged.fa.gz'
        with gzip.open(self.query, 'wt') as dest:
            dest.write('>merged\n' + self.sequence[300:550] + '\n>unmapped\n' + 'N' * 250 + '\n')
        self.merged_manifest = self.flash / 'merged.json'
        self.merged_manifest.write_text(json.dumps([str(self.query)]))
        self.output_dir = self.root / 'tmp/k59/flank_mapping'
        self.sam = self.output_dir / 'library.0.sam'
        self.manifest = self.output_dir / 'mapping.json'

    def test_version_and_real_sam_with_requested_options(self):
        self.launcher.map_merged_reads(59)
        metadata = json.loads(self.manifest.read_text())
        config = metadata['configuration']
        self.assertEqual(config['version'], '2.31-r1302')
        self.assertEqual(config['parameters'], ['-x', 'sr', '-a', '-Y', '-t', '2', '-I', '55G'])
        self.assertEqual(config['reference']['path'], str(self.reference))
        lines = self.sam.read_text().splitlines()
        self.assertTrue(any(line.startswith('@SQ\tSN:reference\tLN:2000') for line in lines))
        records = [line.split('\t') for line in lines if not line.startswith('@')]
        mapped = [r for r in records if r[0] == 'merged' and not int(r[1]) & 4]
        self.assertEqual(len(mapped), 1)
        self.assertEqual(mapped[0][2:4], ['reference', '301'])
        self.assertEqual(mapped[0][5], '250M')
        self.assertTrue(any(r[0] == 'unmapped' and int(r[1]) & 4 for r in records))
        self.assertFalse((self.output_dir / 'library.0.sam.tmp').exists())

    def test_multiple_libraries_are_mapped_separately(self):
        self.merged_manifest.write_text(json.dumps([str(self.query)] * 2))
        self.launcher.map_merged_reads(59)
        self.assertTrue((self.output_dir / 'library.1.sam').is_file())
        for path in self.output_dir.glob('*.sam'):
            records = [line for line in path.read_text().splitlines() if not line.startswith('@')]
            self.assertTrue(all(not int(line.split('\t')[1]) & 1 for line in records))

    def test_resume_reuses_completed_mapping(self):
        self.launcher.map_merged_reads(59)
        self.launcher.check_point._logged_checkpoint = 100
        with mock.patch.object(self.launcher, 'run_sub_command', side_effect=AssertionError('Mapper reran')):
            self.launcher.map_merged_reads(59)

    def test_resume_rebuilds_missing_sam(self):
        self.launcher.map_merged_reads(59)
        self.sam.unlink()
        self.launcher.check_point._logged_checkpoint = 100
        self.launcher.map_merged_reads(59)
        self.assertTrue(self.sam.is_file())

    def test_resume_invalidates_changed_reference(self):
        self.launcher.map_merged_reads(59)
        old = json.loads(self.manifest.read_text())
        self.reference.write_text('>reference\n' + self.sequence + 'A\n')
        self.launcher.check_point._logged_checkpoint = 100
        self.launcher.map_merged_reads(59)
        self.assertNotEqual(old, json.loads(self.manifest.read_text()))
        self.assertIn('LN:2001', self.sam.read_text())

    def test_failure_does_not_publish_partial_sam(self):
        def fail(cmd, *args, **kwargs):
            Path(cmd[cmd.index('-o') + 1]).write_text('partial SAM')
            raise SystemExit(2)
        self.launcher.run_sub_command = fail
        with self.assertRaises(SystemExit):
            self.launcher.map_merged_reads(59)
        self.assertFalse(self.sam.exists())
        self.assertFalse(self.manifest.exists())
        self.assertFalse((self.output_dir / 'library.0.sam.tmp').exists())

    def test_missing_query_is_reported(self):
        self.query.unlink()
        with self.assertRaises(self.launcher.Usage):
            self.launcher.map_merged_reads(59)

    def test_missing_reference_is_reported(self):
        self.reference.unlink()
        with self.assertRaises(self.launcher.Usage):
            self.launcher.map_merged_reads(59)

    def test_empty_reference_is_skipped(self):
        self.reference.write_text('')
        self.launcher.map_merged_reads(59)
        self.assertFalse(self.sam.exists())

    def test_empty_query_produces_header_only_sam(self):
        with gzip.open(self.query, 'wt'):
            pass
        self.launcher.map_merged_reads(59)
        self.assertTrue(all(line.startswith('@') for line in self.sam.read_text().splitlines()))

    def test_no_merged_libraries_is_skipped(self):
        self.merged_manifest.unlink()
        self.launcher.map_merged_reads(59)
        self.assertFalse(self.output_dir.exists())

    def test_merging_disabled_is_skipped(self):
        self.opt.merge_pairs = False
        self.launcher.map_merged_reads(59)
        self.assertFalse(self.output_dir.exists())

    def test_cleanup_preserves_mapping_only(self):
        self.launcher.map_merged_reads(59)
        junk = self.output_dir.parent / 'discard.bin'
        junk.write_bytes(b'temporary')
        self.launcher.remove_temporary_files()
        self.assertTrue(self.sam.is_file())
        self.assertTrue(self.query.is_file())
        self.assertFalse(junk.exists())

    def test_custom_tmp_keeps_mapping_under_output_dir(self):
        custom = self.root / 'scratch'
        custom.mkdir()
        self.opt.temp_dir = str(custom)
        self.launcher.map_merged_reads(59)
        self.launcher.remove_temporary_files()
        self.assertTrue(self.sam.is_file())
        self.assertFalse(custom.exists())


class MappingWorkflowOrderTest(unittest.TestCase):
    def run_workflow(self, k_list, merge_pairs=True, no_local=False):
        with tempfile.TemporaryDirectory() as directory:
            launcher = load_launcher()
            launcher.opt.out_dir = directory
            launcher.opt.temp_dir = directory
            launcher.opt.k_list = k_list
            launcher.opt.k_min, launcher.opt.k_max = k_list[0], k_list[-1]
            launcher.opt.memory = 100000000
            launcher.opt.keep_tmp_files = True
            launcher.opt.merge_pairs = merge_pairs
            launcher.opt.no_local = no_local
            launcher.check_point.set_file(str(Path(directory) / 'checkpoints.txt'))
            events = []
            patches = {name: mock.Mock() for name in [
                'check_bin', 'parse_option', 'setup_output_dir', 'setup_logger',
                'check_and_correct_option', 'check_reads', 'cpu_dispatch',
                'create_library_file', 'build_library', 'merge_paired_reads',
                'count_and_build_read_graph', 'reduce_edges_for_k',
                'build_graph_from_edge_prefix', 'merge_contig_edges_into_next_k',
                'build_graph', 'merge_final', 'remove_temporary_files', 'current_edge_prefix']}
            patches['set_max_k_by_lib'] = mock.Mock(return_value=False)
            patches['build_first_graph'] = mock.Mock(return_value=True)
            for name in ['assemble', 'map_merged_reads', 'extract_contig_ends', 'iterate']:
                patches[name] = mock.Mock(side_effect=lambda *args, _name=name: events.append((_name, args)))
            def extend(cur_k):
                events.append(('extend_contig_ends', (cur_k,)))
                output = Path(directory) / ('k{}.extended.fa'.format(cur_k))
                output.write_text('>extended flag=0 multi=1\nACGTA\n')
                Path(str(output) + '.info').write_text('1 5\n')
                return str(output)
            patches['extend_contig_ends'] = mock.Mock(side_effect=extend)
            original_local = launcher.local_assemble
            def local(*args):
                events.append(('local_assemble', args))
                original_local(*args)
            patches['local_assemble'] = mock.Mock(side_effect=local)
            with mock.patch.multiple(launcher, **patches), mock.patch.object(launcher.opt, 'dump'):
                launcher.main(['megahit'])
            return events

    def test_k59_mapping_precedes_flanking(self):
        events = self.run_workflow([21, 39, 59, 79])
        mapping_index = events.index(('map_merged_reads', (59,)))
        self.assertEqual(events[mapping_index - 2], ('local_assemble', (59, 79)))
        self.assertEqual(events[mapping_index - 1], ('extract_contig_ends', (59,)))
        self.assertEqual(events[mapping_index + 1], ('extend_contig_ends', (59,)))
        self.assertEqual(events[mapping_index + 2], ('iterate', (59, 20)))
        self.assertLess(mapping_index, events.index(('iterate', (59, 20))))

    def test_custom_last_reduction_k_is_used(self):
        events = self.run_workflow([21, 39, 61])
        self.assertIn(('map_merged_reads', (39,)), events)
        self.assertIn(('extract_contig_ends', (39,)), events)

    def test_original_only_workflow_also_uses_consensus_local(self):
        events = self.run_workflow([61, 79])
        self.assertIn(('map_merged_reads', (61,)), events)
        self.assertIn(('extend_contig_ends', (61,)), events)

    def test_no_next_k_does_not_extract_contig_ends(self):
        self.assertFalse(any(name == 'extract_contig_ends' for name, _ in self.run_workflow([21, 39, 59], merge_pairs=False)))

    def test_consensus_runs_at_every_flanking_transition(self):
        events = self.run_workflow([59, 79, 99, 119])
        for cur_k in [59, 79, 99]:
            self.assertIn(('extract_contig_ends', (cur_k,)), events)
            self.assertIn(('map_merged_reads', (cur_k,)), events)
            self.assertIn(('extend_contig_ends', (cur_k,)), events)
            self.assertLess(events.index(('extend_contig_ends', (cur_k,))),
                            events.index(('iterate', (cur_k, 20))))

    def test_no_local_disables_local_publication_at_every_transition(self):
        events = self.run_workflow([59, 79, 99], no_local=True)
        self.assertFalse(any(name == 'local_assemble' for name, _ in events))


if __name__ == '__main__':
    unittest.main()
