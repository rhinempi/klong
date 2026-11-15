"""Native runtime integration. Python is the test harness, never the workload."""
import bz2
import collections
import gzip
import hashlib
import json
import os
from pathlib import Path
import random
import re
import shutil
import subprocess
import tempfile
import unittest

from test_flash import binary_sequences, fastq

NATIVE = os.environ.get('MEGAHIT_TEST_NATIVE')
ZLIB_ONLY = os.environ.get('MEGAHIT_TEST_NATIVE_ZLIB')


def fasta(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as source:
        lines = source.read().splitlines()
    return dict(zip(lines[::2], lines[1::2]))


@unittest.skipUnless(NATIVE, 'set MEGAHIT_TEST_NATIVE to the native megahit executable')
class NativeWorkflowTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.workspace = tempfile.TemporaryDirectory(prefix='native workflow ')
        cls.root = Path(cls.workspace.name)
        cls.binary = cls.root / 'megahit'
        shutil.copy2(NATIVE, cls.binary)
        cls.environment = os.environ.copy()
        cls.environment['PATH'] = ''
        cls.environment['PYTHONHOME'] = '/no-python-installed'
        cls.environment['PYTHONPATH'] = '/no-python-installed'
        rng = random.Random(9182)
        cls.genome = ''.join(rng.choice('ACGT') for _ in range(1400))
        complement = str.maketrans('ACGT', 'TGCA')
        cls.left, cls.right, cls.merged = [], [], []
        for i in range(0, len(cls.genome) - 159, 8):
            fragment = cls.genome[i:i + 160]
            cls.left.append(('r{}/1'.format(i), fragment[:100]))
            cls.right.append(('r{}/2'.format(i), fragment[-100:].translate(complement)[::-1]))
            cls.merged.append(fragment)
        cls.first, cls.second = cls.root / 'mate one.fq.gz', cls.root / 'mate two.fq.bz2'
        with gzip.open(cls.first, 'wt') as output:
            output.write(fastq(cls.left))
        with bz2.open(cls.second, 'wt') as output:
            output.write(fastq(cls.right))

    @classmethod
    def tearDownClass(cls):
        cls.workspace.cleanup()

    def invoke(self, arguments, check=True, binary=None):
        return subprocess.run([str(binary or self.binary)] + list(map(str, arguments)),
                              env=self.environment, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              text=True, check=check, timeout=120)

    def assemble(self, label, extra=(), binary=None, reads=None, keep=True, k_list='21,39,59,79'):
        output = self.root / label
        inputs = reads or ['-1', self.first, '-2', self.second]
        self.invoke(inputs + (['--keep-tmp-files'] if keep else []) +
                    (['--k-list', k_list] if k_list else []) + ['--merge-pairs',
                              '--min-contig-len', '80', '-t', '2', '-m', '200000000',
                              '-o', output] + list(extra), binary=binary)
        return output

    def test_native_elf_and_direct_core_dispatch(self):
        with self.binary.open('rb') as source:
            self.assertEqual(source.read(4), b'\x7fELF')
        info = json.loads(self.invoke(['--native-info']).stdout)
        self.assertFalse(info['python_required'])
        self.assertTrue(info['core_linked'])
        self.assertIn('native', self.invoke(['--version']).stdout)
        self.assertEqual(self.invoke(['--core', 'kmax']).stdout.strip(), '255')
        self.assertIn('FLASH', self.invoke(['flash', '--version']).stdout)
        self.assertTrue(self.invoke(['minimap2', '--version']).stdout.strip())

    def test_native_gzip_bzip2_and_concatenated_inputs(self):
        for extension, encoder in [('gz', gzip.compress), ('bz2', bz2.compress)]:
            source = self.root / ('concat.' + extension)
            target = self.root / ('concat.' + extension + '.plain')
            source.write_bytes(encoder(b'first\n') + encoder(b'second\n'))
            self.invoke(['--native-decompress', source, target])
            self.assertEqual(target.read_bytes(), b'first\nsecond\n')
        corrupt = self.root / 'corrupt.gz'
        corrupt.write_bytes(b'\x1f\x8bnot-gzip')
        self.assertNotEqual(self.invoke(['--native-decompress', corrupt, self.root / 'bad'], check=False).returncode, 0)

    def test_paired_workflow_without_python_and_resume(self):
        output = self.assemble('paired output')
        self.assertTrue((output / 'done').exists())
        options = json.loads((output / 'options.json').read_text())
        self.assertEqual(options['workflow_driver'], 'native-cpp-v1')
        self.assertEqual(options['merge_dedup'], 'sort')
        self.assertEqual(options['merge_sort_memory'], 128)
        self.assertFalse(options['flanking_consensus_only'])
        expected = [sequence for pair in zip(self.left, self.right) for _, sequence in pair]
        self.assertEqual(binary_sequences(output / 'tmp/reads.lib.bin'), expected)
        merged = fasta(output / 'tmp/flash/library.0.merged.fa.gz')
        self.assertEqual(collections.Counter(merged.values()), collections.Counter(self.merged))
        log = (output / 'log').read_text()
        self.assertIn('no merge Bloom filters', log)
        assembled = re.findall(r'Assemble contigs from SdBG for k = (\d+)', log)
        order = list(map(int, assembled))
        expected_order = [21, 39, 59, 21, 39, 59, 79]
        # A perfectly resolved synthetic genome can legitimately terminate
        # early when neither iterative edges nor consensus overhangs remain.
        self.assertEqual(order, expected_order[:len(order)])
        self.assertGreaterEqual(len(order), 4)
        self.assertNotIn('--core","local"', log)
        for phase, k in [('', 59), ('flanking/', 21), ('flanking/', 39), ('flanking/', 59)]:
            local = output / ('intermediate_contigs/' + phase + 'k{}.local.fa'.format(k))
            if k != 21 and phase and not local.exists():
                continue
            extended = output / ('tmp/' + phase + 'k{}/k{}.contig_ends.extended.fa'.format(k, k))
            self.assertEqual(local.read_bytes(), extended.read_bytes())
            self.assertEqual(Path(str(local) + '.info').read_bytes(), Path(str(extended) + '.info').read_bytes())
            self.assertTrue((extended.parent / 'flank_mapping/library.0.sam').exists())
        final = output / 'final.contigs.fa'
        self.assertTrue(fasta(final))
        before = final.read_bytes()
        (output / 'done').unlink()
        final.unlink()
        self.invoke(['--continue', '-o', output])
        self.assertEqual(final.read_bytes(), before)
        self.assertIn('Reuse native checkpoint', (output / 'log').read_text())
        self.invoke(['--continue', '-o', output])

    def test_interleaved_gzip_workflow_and_cleanup(self):
        interleaved = self.root / 'interleaved.fq.gz'
        with gzip.open(interleaved, 'wt') as output:
            output.write(fastq([record for pair in zip(self.left, self.right) for record in pair]))
        output = self.assemble('interleaved', reads=['--12', interleaved])
        self.assertEqual(collections.Counter(fasta(output / 'tmp/flash/library.0.merged.fa.gz').values()),
                         collections.Counter(self.merged))
        # Completed-run continuation does not attempt to rebuild removed scratch.
        self.invoke(['--continue', '-o', output])

    def test_hash_merge_backend_remains_selectable(self):
        output = self.assemble('hash merge backend', extra=['--merge-index', 'hash'])
        self.assertEqual(json.loads((output / 'options.json').read_text())['merge_index'], 'hash')
        self.assertIn('Merge index: hash', (output / 'log').read_text())
        self.assertTrue((output / 'done').exists())

    def test_fm_merge_backend_integrates_with_native_workflow(self):
        output = self.assemble('fm merge backend', extra=['--merge-index', 'fm'])
        self.assertEqual(json.loads((output / 'options.json').read_text())['merge_index'], 'fm')
        self.assertIn('Merge index: fm', (output / 'log').read_text())
        self.assertTrue((output / 'done').exists())

    def test_exact_sort_dedup_integrates_with_native_workflow_and_resume(self):
        output = self.assemble('sort merge dedup', extra=['--merge-dedup', 'sort', '--merge-sort-memory', '0.125'])
        options = json.loads((output / 'options.json').read_text())
        self.assertEqual(options['merge_dedup'], 'sort')
        self.assertEqual(options['merge_sort_memory'], 0.125)
        self.assertIn('no merge Bloom filters', (output / 'log').read_text())
        self.assertFalse(list((output / 'tmp').rglob('*.sort-*')))
        final = output / 'final.contigs.fa'
        before = final.read_bytes()
        final.unlink()
        (output / 'done').unlink()
        self.invoke(['--continue', '-o', output])
        self.assertEqual(final.read_bytes(), before)

    def test_bloom_dedup_and_legacy_resume_preserve_backend(self):
        output = self.assemble('bloom merge dedup', extra=['--merge-dedup', 'bloom'])
        options_path = output / 'options.json'
        options = json.loads(options_path.read_text())
        self.assertEqual(options['merge_dedup'], 'bloom')
        self.invoke(['--continue', '-o', output])
        self.assertEqual(json.loads(options_path.read_text())['merge_dedup'], 'bloom')
        del options['merge_dedup']
        for key in ('reduction_k_list', 'flanking_k_list', 'legacy_k_list',
                    'reduction_k_list_set', 'flanking_k_list_set'):
            options.pop(key, None)
        options_path.write_text(json.dumps(options))
        final = output / 'final.contigs.fa'
        before = final.read_bytes()
        final.unlink()
        (output / 'done').unlink()
        self.invoke(['--continue', '-o', output])
        self.assertEqual(json.loads(options_path.read_text())['merge_dedup'], 'bloom')
        self.assertEqual(final.read_bytes(), before)

    def test_default_phase_k_lists_and_resume(self):
        output = self.assemble('default phase lists', k_list=None)
        options = json.loads((output / 'options.json').read_text())
        self.assertEqual(options['reduction_k_list'], [21, 29, 39, 59])
        self.assertEqual(options['flanking_k_list'], [21, 31, 41, 61, 81, 101, 141])
        log = (output / 'log').read_text()
        self.assertIn('Reduction k list: [21,29,39,59]', log)
        self.assertIn('Flanking k list: [21,31,41,61,81,101,141]', log)
        final = output / 'final.contigs.fa'
        before = final.read_bytes()
        final.unlink()
        (output / 'done').unlink()
        self.invoke(['--continue', '-o', output])
        self.assertEqual(final.read_bytes(), before)

    def test_consensus_only_flanking_retains_backbone_and_ignores_stale_read_edges(self):
        output = self.assemble('consensus only flanking', k_list=None,
                               extra=['--flanking-consensus-only'])
        options = json.loads((output / 'options.json').read_text())
        self.assertTrue(options['flanking_consensus_only'])
        def verify_commands():
            commands = [json.loads(line[line.index('["--core"'):])
                        for line in (output / 'log').read_text().splitlines()
                        if 'command ' in line and '["--core"' in line]
            self.assertFalse(any(args[1] == 'iterate' for args in commands))
            graphs = [args for args in commands if args[1] == 'seq2sdbg'
                      and '/flanking/' in args[args.index('--output_prefix') + 1]]
            self.assertGreaterEqual(len(graphs), 2)
            for args in graphs:
                self.assertNotIn('--input_prefix', args)
                self.assertIn('--contig', args)
                self.assertIn('--local_contig', args)
            self.assertTrue(any(args[args.index('-k') + 1] == '141' for args in graphs))
        verify_commands()
        self.assertFalse(list((output / 'tmp/flanking').rglob('*.edges.*')))
        # A restarted graph must not pick up old iterative edges even if present.
        prefix = output / 'tmp/flanking/k31/31'
        Path(str(prefix) + '.edges.0').write_bytes(b'not valid read-derived edges')
        Path(str(prefix) + '.edges.info').write_text('not valid edge metadata')
        Path(str(prefix) + '.sdbg_info').unlink()
        (output / 'done').unlink()
        self.invoke(['--continue', '-o', output])
        self.assertTrue(json.loads((output / 'options.json').read_text())['flanking_consensus_only'])
        verify_commands()
        self.assertTrue(fasta(output / 'final.contigs.fa'))

    def test_independent_phase_lists_and_different_starting_k(self):
        output = self.assemble('independent phase lists', extra=[
            '--reduction-k-list', '23,33,53', '--flanking-k-list', '31,41,61,81'])
        options = json.loads((output / 'options.json').read_text())
        self.assertEqual(options['reduction_k_list'], [23, 33, 53])
        self.assertEqual(options['flanking_k_list'], [31, 41, 61, 81])
        log = (output / 'log').read_text()
        self.assertIn('Reduce edges for k = 23 using k = 33', log)
        self.assertIn('Restarting consensus flanking assembly at k = 31 from k = 53', log)
        self.assertTrue((output / 'tmp/k23/23.edges.info').exists())
        self.assertTrue((output / 'tmp/flanking/k31/31.sdbg_info').exists())
        self.assertFalse((output / 'tmp/k21').exists())

    def test_phase_override_is_independent_of_shared_option_order(self):
        for index, args in enumerate([
                ['--reduction-k-list', '23,33,53', '--k-list', '21,39,59,79'],
                ['--k-list', '21,39,59,79', '--flanking-k-list', '31,41,61']]):
            output = self.assemble('phase override order ' + str(index), extra=args, k_list=None)
            options = json.loads((output / 'options.json').read_text())
            self.assertEqual(options['reduction_k_list'], [23, 33, 53] if index == 0 else [21, 39, 59])
            self.assertEqual(options['flanking_k_list'], [21, 39, 59, 79] if index == 0 else [31, 41, 61])

    def test_reduction_keeps_counted_edges_without_preliminary_graphs(self):
        output = self.assemble('edges-only reduction', extra=['-t', '1'])
        commands = [json.loads(line[line.index('["--core"'):])
                    for line in (output / 'log').read_text().splitlines()
                    if 'command ' in line and '["--core"' in line]
        graphs = [args for args in commands if args[1] == 'seq2sdbg']
        reductions = [args for args in commands if args[1:4] == ['prefixbf', '--mode', 'reduce']]
        self.assertTrue(reductions)
        for args in reductions:
            self.assertEqual(args[args.index('--host_mem') + 1], '200000000')
        for k in (39, 59):
            prefix = output / 'tmp/k{}/read_count/{}'.format(k, k)
            self.assertTrue(Path(str(prefix) + '.edges.info').exists())
            self.assertTrue(list(prefix.parent.glob(str(k) + '.edges.[0-9]*')))
            self.assertFalse(Path(str(prefix) + '.sdbg_info').exists())
            self.assertFalse(any(args[args.index('--output_prefix') + 1] == str(prefix)
                                 for args in graphs))
        for k in (21, 39, 59):
            self.assertTrue((output / 'tmp/k{}/assembly_graph/{}.sdbg_info'.format(k, k)).exists())

        # Native checkpoints from before this optimization may contain the
        # retired stages; resume must not require or rebuild their graphs.
        checkpoint = output / 'checkpoints.native.json'
        state = json.loads(checkpoint.read_text())
        state.update({'read-count-graph-k39': True, 'read-count-graph-k59': True})
        checkpoint.write_text(json.dumps(state))
        final = output / 'final.contigs.fa'
        before = final.read_bytes()
        final.unlink()
        (output / 'done').unlink()
        previous_log = (output / 'log').read_text()
        self.invoke(['--continue', '-o', output])
        self.assertEqual(final.read_bytes(), before)
        state = json.loads(checkpoint.read_text())
        self.assertFalse(any(key.startswith('read-count-graph-') for key in state))
        self.assertNotIn('Build graph for k =', (output / 'log').read_text()[len(previous_log):])

    @unittest.skipUnless(ZLIB_ONLY, 'set MEGAHIT_TEST_NATIVE_ZLIB to test the no-pigz executable')
    def test_native_zlib_fallback(self):
        binary = self.root / 'megahit-zlib'
        shutil.copy2(ZLIB_ONLY, binary)
        self.assertFalse(json.loads(self.invoke(['--native-info'], binary=binary).stdout)['bundled_pigz'])
        output = self.assemble('zlib-only', binary=binary)
        self.assertIn('native C zlib', (output / 'log').read_text())
        self.assertEqual(collections.Counter(fasta(output / 'tmp/flash/library.0.merged.fa.gz').values()),
                         collections.Counter(self.merged))

    def test_no_local_original_only_and_single_end(self):
        source = self.root / 'single end.fa'
        source.write_text(''.join('>s{}\n{}\n'.format(i, seq) for i, (_, seq) in enumerate(self.left)))
        output = self.root / 'single-only'
        self.invoke(['-r', source, '--no-local', '--k-list', '61,79', '--min-contig-len', '80',
                     '-t', '2', '-m', '200000000', '-o', output])
        self.assertTrue((output / 'done').exists())
        self.assertFalse((output / 'tmp/reads.lib.bin').exists())
        self.assertFalse((output / 'tmp/flash').exists())
        self.assertNotIn('Restarting consensus', (output / 'log').read_text())

    def test_export_preservation_and_custom_scratch_cleanup(self):
        scratch = self.root / 'custom scratch'
        output = self.assemble('cleaned output', extra=['--tmp-dir', scratch], keep=False)
        options = json.loads((output / 'options.json').read_text())
        self.assertFalse(Path(options['temp_dir']).exists())
        self.assertTrue((output / 'tmp/flash/library.0.merged.fa.gz').exists())
        self.assertTrue((output / 'tmp/k59/k59.contig_ends.extended.fa').exists())
        self.assertTrue((output / 'tmp/flanking/k21/flank_mapping/library.0.sam').exists())
        self.assertTrue((output / 'done').exists())

    def test_corrupt_compressed_input_does_not_publish_completion(self):
        corrupt = self.root / 'bad-input.fq.gz'
        corrupt.write_bytes(b'\x1f\x8binvalid data')
        output = self.root / 'failed input'
        result = self.invoke(['-r', corrupt, '--k-list', '21', '-m', '200000000', '-t', '1', '-o', output], check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse((output / 'done').exists())

    def test_failed_rerun_invalidates_old_checkpoint_tail(self):
        output = self.assemble('resume invalidation')
        (output / 'done').unlink()
        (output / 'final.contigs.fa').unlink()
        library = output / 'tmp/reads.lib.bin'
        library.unlink()
        blocker = output / 'tmp/reads.lib.tmp'
        blocker.mkdir()
        failed = self.invoke(['--continue', '-o', output], check=False)
        self.assertNotEqual(failed.returncode, 0)
        state = json.loads((output / 'checkpoints.native.json').read_text())
        self.assertFalse(state.get('read-library', False))
        blocker.rmdir()
        library.write_bytes(b'partial data from an interrupted stage')
        self.invoke(['--continue', '-o', output])
        expected = [sequence for pair in zip(self.left, self.right) for _, sequence in pair]
        self.assertEqual(binary_sequences(library), expected)
        self.assertTrue((output / 'done').exists())

    def test_rejects_legacy_checkpoints_and_bad_options(self):
        for options in [['--k-step', '0'], ['--k-list', '21,21'], ['--flash-min-overlap', '0'],
                        ['--merge-index', 'invalid'], ['--merge-dedup', 'invalid'],
                        ['--merge-sort-memory', '0'], ['--reduction-k-list', '21,61'],
                        ['--flanking-k-list', '21,63'], ['--flanking-k-list', '21,31,31'],
                        ['--flanking-k-list', '31,21'], ['--reduction-k-list', '22'],
                        ['--flanking-k-list', ',']]:
            result = self.invoke(['-r', self.first, '-o', self.root / ('bad-' + options[0][2:])] + options, check=False)
            self.assertNotEqual(result.returncode, 0)
        legacy = self.root / 'legacy'
        legacy.mkdir()
        (legacy / 'options.json').write_text('{}')
        result = self.invoke(['--continue', '-o', legacy], check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Legacy Python checkpoints', result.stderr)
