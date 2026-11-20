"""Native runtime integration. Python is the test harness, never the workload."""
import bz2
import collections
import csv
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

    def assemble(self, label, extra=(), binary=None, reads=None, keep=True, k_list='21,39,59,79', explicit_merge=True):
        output = self.root / label
        inputs = reads or ['-1', self.first, '-2', self.second]
        self.invoke(inputs + (['--keep-tmp-files'] if keep else []) +
                    (['--k-list', k_list] if k_list else []) + (['--merge-pairs'] if explicit_merge else []) + [
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
        self.assertEqual(collections.Counter(binary_sequences(output / 'tmp/reads.lib.bin')),
                         collections.Counter(self.merged))
        merged = fasta(output / 'tmp/flash/library.0.merged.fa.gz')
        self.assertEqual(collections.Counter(merged.values()), collections.Counter(self.merged))
        log = (output / 'log').read_text()
        commands = [json.loads(line[line.index('["--core"'):]) for line in log.splitlines()
                    if 'command ' in line and '["--core"' in line]
        stages = [args[1] for args in commands]
        self.assertLess(stages.index('flash'), stages.index('buildlib'))
        self.assertLess(stages.index('buildlib'), stages.index('count'))
        consensus = [args for args in commands if args[1] == 'endconsensus']
        self.assertTrue(consensus)
        self.assertTrue(all(args[args.index('-t') + 1] == str(options['num_cpu_threads']) for args in consensus))
        self.assertIn('Reduced-edge exact dedup:', log)
        assembled = re.findall(r'Assemble contigs from SdBG for k = (\d+)', log)
        order = list(map(int, assembled))
        expected_order = [59, 39, 21, 21, 39, 59, 79]
        # A perfectly resolved synthetic genome can legitimately terminate
        # early when neither iterative edges nor consensus overhangs remain.
        self.assertEqual(order, expected_order[:len(order)])
        self.assertGreaterEqual(len(order), 4)
        self.assertNotIn('--core","local"', log)
        for phase, k in [('', 21), ('flanking/', 21), ('flanking/', 39), ('flanking/', 59)]:
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
        self.assertEqual(collections.Counter(binary_sequences(output / 'tmp/reads.lib.bin')),
                         collections.Counter(self.merged))
        # Completed-run continuation does not attempt to rebuild removed scratch.
        self.invoke(['--continue', '-o', output])

    def test_flash_preprocessing_is_enabled_by_default(self):
        output = self.assemble('default flash preprocessing', explicit_merge=False)
        self.assertTrue(json.loads((output / 'options.json').read_text())['merge_pairs'])
        self.assertEqual(collections.Counter(binary_sequences(output / 'tmp/reads.lib.bin')),
                         collections.Counter(self.merged))
        self.assertTrue((output / 'tmp/flash/assembly_inputs.json').exists())

    def test_merged_unmerged_and_explicit_single_end_reads_are_retained(self):
        rng = random.Random(3192)
        unmerged = [''.join(rng.choice('ACGT') for _ in range(100)) for _ in range(2)]
        first, second = self.root / 'mixed.1.fq.gz', self.root / 'mixed.2.fq.gz'
        for path, records in [(first, self.left + [('not-merged/1', unmerged[0])]),
                              (second, self.right + [('not-merged/2', unmerged[1])])]:
            with gzip.open(path, 'wt') as target:
                target.write(fastq(records))
        single = self.root / 'explicit-single.fa'
        single.write_text('>single\n' + self.genome[:140] + '\n')
        output = self.assemble('mixed merged and unmerged', reads=['-1', first, '-2', second, '-r', single],
                               extra=['--flash-min-overlap', '30'])
        expected = self.merged + unmerged + [self.genome[:140]]
        self.assertEqual(collections.Counter(binary_sequences(output / 'tmp/reads.lib.bin')),
                         collections.Counter(expected))
        exports = output / 'tmp/flash'
        inputs = json.loads((exports / 'assembly_inputs.json').read_text())
        self.assertEqual(len(inputs), 3)
        self.assertEqual(list(fasta(exports / 'library.0.unmerged.1.fa.gz').values()), [unmerged[0]])
        self.assertEqual(list(fasta(exports / 'library.0.unmerged.2.fa.gz').values()), [unmerged[1]])
        info = (output / 'tmp/reads.lib.lib_info').read_text().splitlines()
        self.assertTrue(all(line.split()[-1] == '0' for line in info[2::2]))
        mapping = json.loads((output / 'tmp/k21/flank_mapping/mapping.json').read_text())
        self.assertEqual([entry['path'] for entry in mapping['configuration']['queries']], inputs)
        self.assertEqual(len(mapping['outputs']), 3)

    def test_all_unmerged_pairs_still_build_library_and_map_for_consensus(self):
        # The synthetic pairs overlap by only 40 bp, so requiring 90 prevents merging.
        output = self.assemble('all unmerged', extra=['--flash-min-overlap', '90', '--flash-max-overlap', '100'])
        expected = [seq for pair in zip(self.left, self.right) for _, seq in pair]
        self.assertEqual(collections.Counter(binary_sequences(output / 'tmp/reads.lib.bin')),
                         collections.Counter(expected))
        self.assertEqual(fasta(output / 'tmp/flash/library.0.merged.fa.gz'), {})
        inputs = json.loads((output / 'tmp/flash/assembly_inputs.json').read_text())
        self.assertEqual(len(inputs), 2)
        self.assertTrue(all('.unmerged.' in path for path in inputs))
        mapping = json.loads((output / 'tmp/k21/flank_mapping/mapping.json').read_text())
        self.assertEqual([entry['path'] for entry in mapping['configuration']['queries']], inputs)

    def test_no_merge_pairs_preserves_original_mates(self):
        output = self.assemble('original mates opt out', extra=['--no-merge-pairs'])
        expected = [seq for pair in zip(self.left, self.right) for _, seq in pair]
        self.assertEqual(binary_sequences(output / 'tmp/reads.lib.bin'), expected)
        self.assertFalse((output / 'tmp/flash/assembly_inputs.json').exists())
        self.assertFalse(json.loads((output / 'options.json').read_text())['merge_pairs'])

    def test_old_raw_read_library_checkpoint_is_rebuilt(self):
        output = self.assemble('raw read library upgrade')
        state_path = output / 'checkpoints.native.json'
        state = json.loads(state_path.read_text())
        state['flash-merged-reads'] = state.pop('flash-assembly-inputs-v2')
        state_path.write_text(json.dumps(state))
        (output / 'tmp/reads.lib.bin').write_bytes(b'old raw read library')
        (output / 'done').unlink()
        (output / 'final.contigs.fa').unlink()
        self.invoke(['--continue', '-o', output])
        self.assertEqual(collections.Counter(binary_sequences(output / 'tmp/reads.lib.bin')),
                         collections.Counter(self.merged))
        state = json.loads(state_path.read_text())
        self.assertIn('flash-assembly-inputs-v2', state)
        self.assertNotIn('flash-merged-reads', state)

    def test_hash_merge_setting_does_not_enable_overlap_search(self):
        output = self.assemble('hash merge backend', extra=['--merge-index', 'hash'])
        self.assertEqual(json.loads((output / 'options.json').read_text())['merge_index'], 'hash')
        self.assertNotIn('Merge k =', (output / 'log').read_text())
        self.assertTrue((output / 'done').exists())

    def test_fm_merge_setting_does_not_enable_overlap_search(self):
        output = self.assemble('fm merge backend', extra=['--merge-index', 'fm'])
        self.assertEqual(json.loads((output / 'options.json').read_text())['merge_index'], 'fm')
        self.assertNotIn('Merge k =', (output / 'log').read_text())
        self.assertTrue((output / 'done').exists())

    def test_exact_sort_dedup_integrates_with_native_workflow_and_resume(self):
        output = self.assemble('sort merge dedup', extra=['--merge-dedup', 'sort', '--merge-sort-memory', '0.125'])
        options = json.loads((output / 'options.json').read_text())
        self.assertEqual(options['merge_dedup'], 'sort')
        self.assertEqual(options['merge_sort_memory'], 0.125)
        self.assertIn('Reduced-edge exact dedup:', (output / 'log').read_text())
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
        self.assertFalse(options['flanking_end_anchors'])
        self.assertFalse(options['reduction_end_anchors'])
        log = (output / 'log').read_text()
        self.assertIn('Reduction k list: [21,29,39,59]', log)
        self.assertIn('Flanking k list: [21,31,41,61,81,101,141]', log)
        self.assertFalse(list((output / 'intermediate_contigs').rglob('*.anchor-threaded.fa')))
        self.assertFalse(list((output / 'intermediate_contigs').rglob('*.anchor-graph-paths.tsv')))
        commands = [json.loads(line[line.index('["--core"'):]) for line in log.splitlines()
                    if 'command ' in line and '["--core"' in line]
        graphs = [args for args in commands if args[1] == 'seq2sdbg'
                  and '/flanking/' in args[args.index('--output_prefix') + 1]]
        self.assertTrue(graphs)
        full_inputs = [output / 'intermediate_contigs/k21.contigs.fa']
        full_inputs += [output / f'intermediate_contigs/flanking/k{k}.contigs.fa'
                        for k in options['flanking_k_list'][:-1]]
        self.assertEqual([args[args.index('--contig') + 1] for args in graphs],
                         list(map(str, full_inputs[:len(graphs)])))
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
                i = options['flanking_k_list'].index(int(args[args.index('-k') + 1]))
                source = (output / 'intermediate_contigs/k21.contigs.fa' if i == 0 else
                          output / f"intermediate_contigs/flanking/k{options['flanking_k_list'][i - 1]}.contigs.fa")
                self.assertEqual(args[args.index('--contig') + 1], str(source))
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

    def test_end_anchors_restore_both_orientations_and_reject_ambiguous_joins(self):
        rng = random.Random(381)
        dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
        rc = lambda seq: seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
        records = {name: dna(length) for name, length in [
            ('forward', 110), ('reverse', 108), ('short', 38), ('duplicated', 105),
            ('missing', 120), ('loop', 105)]}
        end, k = 24, 17
        def stitched(name):
            seq = records[name]
            return seq if len(seq) <= 2 * end else seq[:end] + seq[-end:]
        left, right = dna(7), dna(8)
        raw_records = [
            left + stitched('forward') + right,
            rc(left + stitched('reverse') + right),
            left + stitched('short') + right,
            left + stitched('duplicated'), stitched('duplicated') + right,
            left + stitched('loop') + right,
        ]
        original, raw, result = [self.root / name for name in ('anchor-original.fa', 'anchor-raw.fa', 'anchor-restored.fa')]
        original.write_text(''.join(f'>{name} flag={3 if name == "loop" else 0} multi=3 len={len(seq)}\n{seq}\n'
                                    for name, seq in records.items()))
        raw.write_text(''.join(f'>unitig{i}\n{seq}\n' for i, seq in enumerate(raw_records)))
        self.invoke(['--native-restore-anchors', original, raw, result, end, k, 1])
        expected_bytes = result.read_bytes()
        self.invoke(['--native-restore-anchors', original, raw, result, end, k, 4])
        self.assertEqual(result.read_bytes(), expected_bytes)
        restored = {header.split()[0][1:]: seq for header, seq in fasta(result).items()}
        for name in ('forward', 'reverse', 'short'):
            records[name] = left + records[name] + right
        self.assertEqual(restored, records)
        self.assertEqual(Path(str(result) + '.info').read_text().split(),
                         [str(len(records)), str(sum(map(len, records.values())))])

    def test_multi_anchor_restoration_inserts_full_contigs_and_bridges_once(self):
        rng = random.Random(676)
        dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
        rc = lambda seq: seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
        end, k = 33, 31
        records = {name: dna(100) for name in
                   ('b', 'c', 'a', 'before', 'ambiguous', 'after', 'gap_left', 'gap_right',
                    'overlap_left', 'overlap_right')}
        records['overlap_right'] = records['overlap_left'][-17:] + records['overlap_right'][17:]
        def stitched(name):
            seq = records[name]
            return seq[:end] + seq[-end:]
        left, right, bridge = dna(8), dna(12), dna(9)
        unitigs = [
            left + stitched('a') + bridge + rc(stitched('b')) + stitched('c') + right,
            stitched('before') + stitched('ambiguous') + stitched('after'),
            rc(stitched('ambiguous')),
            stitched('gap_left') + 'N' + stitched('gap_right'),
            stitched('overlap_left') + stitched('overlap_right')[17:],
        ]
        original, raw, result = [self.root / name for name in ('chain-original.fa', 'chain-raw.fa', 'chain-restored.fa')]
        original.write_text(''.join(f'>{name} flag=0 multi=3 len={len(seq)}\n{seq}\n'
                                    for name, seq in records.items()))
        raw.write_text(''.join(f'>unitig{i}\n{seq}\n' for i, seq in enumerate(unitigs)))
        joined = left + records.pop('a') + bridge + rc(records['b']) + records.pop('c') + right
        records['b'] = joined
        for threads in (1, 8):
            self.invoke(['--native-restore-anchors', original, raw, result, end, k, threads])
            restored = fasta(result)
            self.assertEqual({header.split()[0][1:]: seq for header, seq in restored.items()}, records)
            self.assertEqual(next(iter(restored)), f'>b flag=0 multi=3 len={len(joined)}')
            self.assertEqual(Path(str(result) + '.info').read_text().split(),
                             [str(len(records)), str(sum(map(len, records.values())))])

    def test_packed_anchor_word_boundaries_and_global_uniqueness(self):
        rng = random.Random(728)
        dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
        original, raw, result = [self.root / name for name in
                                 ('packed-original.fa', 'packed-unitigs.fa', 'packed-restored.fa')]
        for k in (32, 141, 255):
            with self.subTest(k=k):
                end = max(k, 35)
                records = {name: dna(2 * end + 25) for name in ('unique', 'duplicate', 'same_seed')}
                def anchor(name):
                    return records[name][:end] + records[name][-end:]
                left, right = dna(3), dna(5)
                changed = anchor('same_seed')
                changed = changed[:-1] + ('A' if changed[-1] != 'A' else 'T')
                unitigs = [dna(25) + 'N' + left + anchor('unique') + right + 'N' + dna(7),
                           left + anchor('duplicate'),
                           left + changed + right]
                # Put the second placement in a later scan batch and on the opposite strand.
                unitigs += [dna(60) for _ in range(260)]
                unitigs += [(anchor('duplicate') + right).translate(str.maketrans('ACGT', 'TGCA'))[::-1]]
                original.write_text(''.join(f'>{name} flag=0 multi=3\n{seq}\n' for name, seq in records.items()))
                raw.write_text(''.join(f'>unitig{i}\n{seq}\n' for i, seq in enumerate(unitigs)))
                for threads in (1, 8):
                    self.invoke(['--native-restore-anchors', original, raw, result, end, k, threads])
                    restored = {header.split()[0][1:]: seq for header, seq in fasta(result).items()}
                    self.assertEqual(restored, dict(records, unique=left + records['unique'] + right))

    def test_short_seeds_count_substitutions_and_report_partial_anchors(self):
        rng = random.Random(1570)
        dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
        rc = lambda seq: seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
        records = {name: dna(120) for name in ('forward', 'reverse', 'fragment', 'terminal', 'tied', 'missing')}

        def mutated(name, *positions):
            seq = list(records[name][:40] + records[name][-40:])
            for pos in positions:
                seq[pos] = 'C' if seq[pos] == 'A' else 'A'
            return ''.join(seq)

        left, right = dna(9), dna(13)
        unitigs = [left + mutated('forward', 27, 44) + right,
                   rc(left + mutated('reverse', 27, 44) + right),
                   mutated('fragment', 27)[12:54],
                   left + mutated('terminal', 4) + right,
                   left + mutated('tied', 27), mutated('tied', 27) + right]
        original, raw, result, report = [self.root / ('short-seeds' + suffix) for suffix in
                                          ('-source.fa', '-unitigs.fa', '-restored.fa', '.tsv')]
        original.write_text(''.join(f'>{name} flag=0 multi=3\n{seq}\n' for name, seq in records.items()))
        raw.write_text(''.join(f'>unitig{i}\n{seq}\n' for i, seq in enumerate(unitigs)))
        args = ['--native-restore-anchors', original, raw, result, 40, 21]
        for threads, tolerance in ((1, 0), (4, 1), (1, 2), (8, 2)):
            self.invoke(args + [threads, tolerance, report])
            expected = dict(records)
            if tolerance == 2:
                for name in ('forward', 'reverse'):
                    expected[name] = left + records[name] + right
            self.assertEqual({h.split()[0][1:]: s for h, s in fasta(result).items()}, expected)
            with report.open() as source:
                hits = {row['source_contig']: row for row in csv.DictReader(source, delimiter='\t')}
            for name, kind, aligned, mm, unaligned in [
                ('forward', 'full', '80', '2', '0'), ('reverse', 'full', '80', '2', '0'),
                ('fragment', 'fragment', '42', '1', '38'), ('terminal', 'full', '80', '1', '0'),
                ('missing', 'no_hit', '0', '.', '80'),
            ]:
                self.assertEqual(tuple(hits[name][key] for key in ('match', 'aligned_nt', 'mismatches', 'unaligned_nt')),
                                 (kind, aligned, mm, unaligned))
            self.assertEqual(hits['reverse']['strand'], '-')

    def test_cleaned_graph_anchor_paths_cross_branches_and_reject_pruned_edges(self):
        rng = random.Random(1275)
        dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
        rc = lambda seq: seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
        left, middle, right = dna(30), dna(20), dna(30)
        other = dna(20)
        other = ('C' if middle[0] == 'A' else 'A') + other[1:-1] + ('C' if middle[-1] == 'A' else 'A')
        tip = left + next(b for b in 'ACGT' if b not in (middle[0], other[0])) + dna(5)
        stitched = left + middle + right
        original = stitched[:40] + dna(30) + stitched[40:]
        left_flank, right_flank = dna(9), dna(11)
        cycle = dna(37)
        inputs = self.root / 'graph path inputs.fa'
        records = [left_flank + stitched + right_flank + base + dna(22) for base in ('A', 'C')]
        records += [left + other + right, tip, cycle + cycle[:21]]
        inputs.write_text(''.join(f'>input{i} flag=0 multi=10\n{s}\n' for i, s in enumerate(records)))
        Path(str(inputs) + '.info').write_text(f'{len(records)} {sum(map(len, records))}\n')
        graph = self.root / 'graph path sdbg'
        self.invoke(['--core', 'seq2sdbg', '--contig', inputs, '--output_prefix', graph,
                     '--host_mem', 200000000, '--num_cpu_threads', 2, '-k', 21, '--kmer_from', 21])

        queries = self.root / 'graph path queries.fa'
        queries.write_text(f'>forward\n{original}\n>reverse\n{rc(original)}\n>tip\n{tip}\n'
                           f'>absent\n{dna(40)}\n>short\n{left[:21]}\n'
                           f'>ambiguous\n{left}N\n>loop flag=3\n{stitched}\n'
                           f'>cycle_path\n{cycle[:31]}\n')
        report = self.root / 'graph path hits.tsv'
        threaded = self.root / 'threaded graph paths.fa'
        args = ['--core', 'assemble', '-s', graph, '-o', self.root / 'graph path unitigs',
                '--prune_level', 0, '--bubble_level', 0, '--disconnect_ratio', 0,
                '--anchor_source', queries, '--anchor_end_length', 40, '--anchor_path_report', report]

        self.invoke(args + ['-t', 1, '--max_tip_len', 0])
        with report.open() as source:
            hits = {row['source_contig']: row for row in csv.DictReader(source, delimiter='\t')}
        self.assertEqual(hits['tip']['match'], 'full')
        self.assertFalse(threaded.exists())
        args += ['--anchor_threaded_contigs', threaded]
        self.invoke(args + ['-t', 1, '--max_tip_len', 8])
        cleaned = report.read_bytes()
        threaded_bytes = threaded.read_bytes()
        self.invoke(args + ['-t', 8, '--max_tip_len', 8])
        self.assertEqual(report.read_bytes(), cleaned)
        self.assertEqual(threaded.read_bytes(), threaded_bytes)
        with report.open() as source:
            hits = {row['source_contig']: row for row in csv.DictReader(source, delimiter='\t')}
        for name in ('forward', 'reverse'):
            self.assertEqual(tuple(hits[name][k] for k in ('match', 'matched_prefix_nt', 'branch_crossings')),
                             ('full', '80', '2'))
        self.assertEqual(hits['tip']['match'], 'missing')
        self.assertEqual(hits['tip']['matched_prefix_nt'], '30')
        self.assertEqual(hits['absent']['matched_prefix_nt'], '0')
        self.assertEqual(hits['short']['match'], 'too_short')
        self.assertEqual(hits['ambiguous']['match'], 'non_acgt')
        self.assertEqual(hits['loop']['match'], 'circular')
        paths = {h.split()[0]: (h, s) for h, s in fasta(threaded).items()}
        self.assertEqual(set(paths), {'>forward', '>reverse', '>cycle_path'})
        expected = left_flank + stitched + right_flank
        self.assertEqual(paths['>forward'][1], expected)
        self.assertEqual(paths['>reverse'][1], rc(expected))
        self.assertIn('left_extension=9 right_extension=11', paths['>forward'][0])
        self.assertIn('left_extension=11 right_extension=9', paths['>reverse'][0])
        self.assertEqual(len(paths['>cycle_path'][1]), 37 + 21)
        self.assertTrue(paths['>cycle_path'][1].endswith(cycle[:31]))

        queries.write_text(f'>forward flag=0 multi=10 len={len(original)}\n{original}\n')
        self.invoke(args + ['-t', 8, '--max_tip_len', 8])
        restored = self.root / 'restored graph path.fa'
        self.invoke(['--native-restore-anchors', queries, threaded, restored, 40, 21, 8])
        self.assertEqual(list(fasta(restored).values()), [left_flank + original + right_flank])

    def test_reduction_excludes_unanchored_unitigs_even_if_they_survive_assembly(self):
        isolated = ''.join(random.Random(5540).choices('ACGT', k=48))
        reads = self.root / 'standalone reduction input.fa'
        reads.write_text(''.join(f'>r{i}\n{seq}\n' for i, seq in
                                 enumerate(self.merged + [isolated] * 4)))
        output = self.assemble('standalone reduction output', reads=['-r', reads],
                               explicit_merge=False, k_list=None,
                               extra=['--reduction-k-list', '31,41,51', '--flanking-k-list', '31',
                                      '--reduction-end-anchors', '--no-anchor-graph-thread',
                                      '--no-local', '--max-tip-len', '0'])
        self.assertFalse((output / 'intermediate_contigs/k41.anchor-threaded.fa').exists())
        self.assertFalse((output / 'intermediate_contigs/k41.anchor-graph-paths.tsv').exists())
        seqs = (isolated, isolated.translate(str.maketrans('ACGT', 'TGCA'))[::-1])
        raw = fasta(output / 'tmp/k41/assembly_graph/41.anchor-unitigs.contigs.fa')
        self.assertTrue(any(s in seqs for s in raw.values()))
        for k in (41, 31):
            records = fasta(output / f'intermediate_contigs/k{k}.contigs.fa')
            self.assertFalse(any(s in seqs for s in records.values()))

        checkpoints = output / 'checkpoints.native.json'
        state = json.loads(checkpoints.read_text())
        state['restore-reduction-unitigs-k41'] = state.pop('restore-reduction-anchors-only-k41')
        checkpoints.write_text(json.dumps(state))
        k41 = output / 'intermediate_contigs/k41.contigs.fa'
        with k41.open('a') as contigs:
            contigs.write(f'>reduced_k41_0 flag=1 multi=4.0000 len=48\n{isolated}\n')
        (output / 'done').unlink()
        self.invoke(['--continue', '-o', output])
        self.assertFalse(any(s in seqs for s in fasta(k41).values()))
        self.assertIn('restore-reduction-anchors-only-k41', json.loads(checkpoints.read_text()))

    def test_mixed_flanking_restore_routes_whole_and_stitched_anchors(self):
        rng = random.Random(8127)
        dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
        short, boundary, long = dna(579), dna(580), dna(581)
        stitched = long[:290] + long[-290:]
        left, right = dna(7), dna(10)
        source, threaded, raw, restored, report = [self.root / ('mix-' + suffix) for suffix in
                                                   ('source.fa', 'threaded.fa', 'raw.fa', 'restored.fa', 'hits.tsv')]
        source.write_text(f'>short\n{short}\n>boundary\n{boundary}\n>long\n{long}\n')
        raw.write_text(f'>raw_short\n{left + short + right + stitched}\n'
                       f'>raw_boundary\n{left + boundary + right}\n')
        threaded.write_text(f'>threaded_long\n{right + stitched + left}\n'
                            f'>threaded_short\n{short}\n>threaded_boundary\n{boundary}\n')
        for mismatches in (0, 1):
            self.invoke(['--native-restore-anchors', source, threaded, restored, 290, 21, 4, mismatches, report, raw])
            self.assertEqual({h.split()[0]: s for h, s in fasta(restored).items()},
                             {'>short': left + short + right, '>boundary': left + boundary + right,
                              '>long': right + long + left})
            with report.open() as hits:
                matches = {r['source_contig']: r['unitig'] for r in csv.DictReader(hits, delimiter='\t')
                           if r['match'] == 'full'}
            self.assertEqual(matches, {'short': 'raw_short', 'boundary': 'raw_boundary', 'long': 'threaded_long'})

    def test_threaded_restoration_prefers_owned_path_over_nested_copies(self):
        rng = random.Random(52124)
        dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
        rc = lambda seq: seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
        names = ('nested', 'divergent', 'beyond', 'repeat', 'unowned')
        originals = {name: dna(120) for name in names}
        anchors = {name: seq[:40] + seq[-40:] for name, seq in originals.items()}
        left, right = dna(6), dna(9)
        owned = {name: left + seq + right for name, seq in anchors.items()}
        records = [('nested_copy', owned['nested'] + 'AGT'),
                   ('nested_left', 'TGA' + owned['nested']),
                   ('nested_reverse', rc('TGA' + owned['nested'] + 'AGT')),
                   ('divergent_copy', left + anchors['divergent'] + ('A' if right[0] != 'A' else 'C') + right[1:]),
                   ('beyond_copy_a', owned['beyond'] + 'A'),
                   ('beyond_copy_c', owned['beyond'] + 'C'),
                   ('unowned_duplicate', owned['unowned'])]
        for name in names:
            header = name + ' anchor_nt=80 left_extension=6' if name != 'unowned' else name
            seq = owned[name] + anchors[name] if name == 'repeat' else owned[name]
            records.append((header, seq))
        source, threaded, raw, restored = [self.root / ('nested-' + name + '.fa') for name in
                                           ('source', 'threaded', 'raw', 'restored')]
        source.write_text(''.join(f'>{name}\n{seq}\n' for name, seq in originals.items()))
        threaded.write_text(''.join(f'>{header}\n{seq}\n' for header, seq in records))
        raw.write_text('')
        for threads, mismatches in ((1, -1), (8, -1), (4, 1)):
            self.invoke(['--native-restore-anchors', source, threaded, restored, 40, 21, threads, mismatches, '', raw])
            self.assertEqual({h.split()[0][1:]: seq for h, seq in fasta(restored).items()},
                             dict(originals, nested=left + originals['nested'] + right))

    def test_best_hits_restore_multiple_anchors_and_keep_overlap_loser(self):
        rng = random.Random(87654)
        dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
        rc = lambda seq: seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
        anchors = {name: dna(80) for name in ('a', 'b', 'c', 'd', 'e', 'missing')}
        anchors['e'] = anchors['d'][-30:] + anchors['e'][30:]
        originals = {name: seq[:40] + dna(35) + seq[40:] for name, seq in anchors.items()}
        flank, bridge = dna(6), dna(11)
        records = [('alternative_a', anchors['a'] + dna(21)),
                   ('alternative_b', dna(12) + rc(anchors['b'])),
                   ('chain', flank + anchors['a'] + bridge + rc(anchors['b']) + anchors['c'] + flank),
                   ('overlap', anchors['d'] + anchors['e'][30:])]
        source, unitigs, result, report = [self.root / ('best-' + name) for name in
                                            ('source.fa', 'raw.fa', 'result.fa', 'matches.tsv')]
        source.write_text(''.join(f'>{name} flag=0\n{seq}\n' for name, seq in originals.items()))
        expected = {'>a': flank + originals['a'] + bridge + rc(originals['b']) + originals['c'] + flank,
                    '>d': originals['d'] + anchors['e'][30:40], '>e': originals['e'],
                    '>missing': originals['missing']}
        for threads, rows in ((1, records), (8, reversed(records))):
            unitigs.write_text(''.join(f'>{name}\n{seq}\n' for name, seq in rows))
            self.invoke(['--native-restore-anchors', source, unitigs, result, 40, 21, threads, 0, report, '', 'best'])
            self.assertEqual({h.split()[0]: seq for h, seq in fasta(result).items()}, expected)
            with report.open() as file:
                full = {r['source_contig']: r['unitig'] for r in csv.DictReader(file, delimiter='\t')
                        if r['match'] == 'full'}
            self.assertEqual([full[name] for name in ('a', 'b', 'c')], ['chain'] * 3)

    def test_gapped_anchor_best_identity_and_cutoff(self):
        rng = random.Random(15147)
        dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
        rc = lambda seq: seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
        def substitute(seq, positions):
            seq = list(seq)
            for p in positions:
                seq[p] = 'A' if seq[p] != 'A' else 'C'
            return ''.join(seq)

        anchor = dna(80)
        original = anchor[:40] + dna(40) + anchor[40:]
        left, right = dna(8), dna(12)
        strong = substitute(anchor, range(14, 66, 6))
        strong = strong[:40] + 'GGG' + strong[40:]
        weak = substitute(anchor, range(9, 72, 2))
        source, targets, restored, report = [self.root / ('gapped-' + name) for name in
                                             ('source.fa', 'raw.fa', 'restored.fa', 'hits.tsv')]
        source.write_text(f'>original\n{original}\n')
        targets.write_text(f'>weak\n{left + weak + right}\n>strong\n{rc(left + strong + right)}\n')
        args = ['--native-restore-anchors', source, targets, restored, 40, 21]
        for threads in (1, 8):
            self.invoke(args + [threads, -1, report, '', 'best', 30])
            with report.open() as matches:
                hit = next(csv.DictReader(matches, delimiter='\t'))
            self.assertEqual(hit['unitig'], 'strong')
            self.assertEqual(hit['strand'], '-')
            self.assertEqual(hit['target_span_nt'], '83')
            self.assertGreaterEqual(int(hit['gap_nt']), 3)
            self.assertGreater(float(hit['identity_pct']), 80)
            self.assertEqual(list(fasta(restored).values()), [left + original + right])

        low_identity = 'A' * 80
        source.write_text(f'>low_identity\n{low_identity}\n')
        targets.write_text('>low\n' + 'A' * 28 + 'C' * 52 + '\n')
        for threshold in (30, 36):
            self.invoke(args + [1, -1, report, '', 'best', threshold])
            with report.open() as matches:
                hit = next(csv.DictReader(matches, delimiter='\t'))
            self.assertEqual(hit['match'], 'full' if threshold == 30 else 'no_hit')
            self.assertEqual(hit['identity_pct'], '35' if threshold == 30 else '.')

    def test_gapped_anchor_option_writes_flanking_alignment_report(self):
        output = self.assemble('gapped flanking', extra=['--flanking-end-anchors', '--anchor-min-identity', '30'],
                               k_list='21,39')
        self.assertEqual(json.loads((output / 'options.json').read_text())['anchor_min_identity'], 30)
        report = output / 'intermediate_contigs/flanking/k21.anchor-alignments.tsv'
        with report.open() as source:
            self.assertIn('identity_pct', next(csv.DictReader(source, delimiter='\t')))
        self.assertTrue(fasta(output / 'final.contigs.fa'))

    def test_end_anchor_flanking_preserves_original_contigs_and_uses_configured_edges(self):
        for consensus_only in (False, True):
            extra = ['--reduction-k-list', '31,41,51', '--flanking-k-list', '21,31,51,91,131,141',
                     '--flanking-end-anchors']
            if consensus_only:
                extra += ['--flanking-consensus-only', '--anchor-max-mismatches', '0', '--anchor-graph-thread']
            else:
                extra += ['--anchor-graph-report']
            output = self.assemble('end anchors ' + str(consensus_only), extra=extra, k_list=None)
            options = json.loads((output / 'options.json').read_text())
            self.assertTrue(options['flanking_end_anchors'])
            self.assertEqual(options['anchor_max_mismatches'], 0 if consensus_only else -1)
            self.assertEqual(options['anchor_graph_thread'], consensus_only)
            original = fasta(output / 'intermediate_contigs/k31.contigs.fa')
            anchors = fasta(output / 'tmp/flanking/k21/21.anchors.fa')
            self.assertCountEqual(anchors.values(), [s if len(s) <= 580 else s[:290] + s[-290:]
                                                      for header, s in original.items() if 'flag=3' not in header])
            previous = {h.split()[0]: s for h, s in original.items()}
            for k in (21, 31, 51, 91, 131, 141):
                report = output / f'intermediate_contigs/flanking/k{k}.anchor-matches.tsv'
                self.assertEqual(report.exists(), consensus_only)
                graph_report = output / f'intermediate_contigs/flanking/k{k}.anchor-graph-paths.tsv'
                self.assertTrue(graph_report.exists())
                threaded = output / f'intermediate_contigs/flanking/k{k}.anchor-threaded.fa'
                self.assertEqual(threaded.exists(), consensus_only)
                restored = {h.split()[0]: s for h, s in
                            fasta(output / f'intermediate_contigs/flanking/k{k}.contigs.fa').items()}
                self.assertLessEqual(set(restored), set(previous))
                self.assertTrue(all(any(s in new or
                                        s.translate(str.maketrans('ACGT', 'TGCA'))[::-1] in new
                                        for new in restored.values()) for s in previous.values()))
                previous = restored
            self.assertEqual((output / 'tmp/flanking/k31/31.edges.info').exists(), not consensus_only)
            self.assertTrue(fasta(output / 'final.contigs.fa'))
            state = json.loads((output / 'checkpoints.native.json').read_text())
            stage = 'restore-mixed-flanking-anchors-v2-k21' if consensus_only else 'restore-best-hit-flanking-anchors-k21'
            self.assertTrue(any(key.endswith(stage) for key in state))
            if consensus_only:
                expected = output / 'expected threaded k21.fa'
                self.invoke(['--native-restore-anchors', output / 'intermediate_contigs/k31.contigs.fa',
                             output / 'intermediate_contigs/flanking/k21.anchor-threaded.fa', expected,
                             290, 21, 2, -1, '', output / 'tmp/flanking/k21/21.anchor-unitigs.contigs.fa'])
                self.assertEqual((output / 'intermediate_contigs/flanking/k21.contigs.fa').read_bytes(),
                                 expected.read_bytes())
            else:
                expected = output / 'expected best hit k21.fa'
                self.invoke(['--native-restore-anchors', output / 'intermediate_contigs/k31.contigs.fa',
                             output / 'tmp/flanking/k21/21.anchor-unitigs.contigs.fa', expected,
                             290, 21, 2, -1, '', '', 'best'])
                self.assertEqual((output / 'intermediate_contigs/flanking/k21.contigs.fa').read_bytes(),
                                 expected.read_bytes())
            (output / 'done').unlink()
            (output / 'final.contigs.fa').unlink()
            if consensus_only:
                threaded.unlink()
            self.invoke(['--continue', '-o', output])
            self.assertTrue((output / 'done').exists())
            self.assertTrue(graph_report.exists())
            self.assertEqual(threaded.exists(), consensus_only)

    def test_reduction_end_anchors_restore_contigs_and_seed_either_flanking_mode(self):
        complement = str.maketrans('ACGT', 'TGCA')
        for anchored_flanking in (False, True):
            extra = ['--reduction-k-list', '21,31,41,51', '--flanking-k-list', '21,31',
                     '--reduction-end-anchors', '--anchor-graph-thread']
            if anchored_flanking:
                extra += ['--flanking-end-anchors', '--flanking-consensus-only']
            output = self.assemble('reduction anchors ' + str(anchored_flanking), extra=extra, k_list=None)
            self.assertTrue(json.loads((output / 'options.json').read_text())['reduction_end_anchors'])
            log = (output / 'log').read_text()
            commands = [json.loads(line[line.index('["--core"'):]) for line in log.splitlines()
                        if 'command ' in line and '["--core"' in line]
            graphs = [args for args in commands if args[1] == 'seq2sdbg']
            self.assertNotIn('--contig', graphs[0])
            self.assertEqual(graphs[0][graphs[0].index('--input_prefix') + 1],
                             str(output / 'tmp/k51/read_count/51'))
            previous = list(fasta(output / 'intermediate_contigs/k51.contigs.fa').values())
            for args, k, source_k in zip(graphs[1:4], (41, 31, 21), (51, 41, 31)):
                self.assertEqual(args[args.index('--input_prefix') + 1], str(output / f'tmp/k{k}/reduced_edges/{k}'))
                self.assertEqual(args[args.index('--kmer_from') + 1], str(source_k))
                anchor_file = output / f'tmp/k{k}/assembly_graph/{k}.anchors.fa'
                self.assertEqual(args[args.index('--contig') + 1], str(anchor_file))
                self.assertTrue(all(len(s) <= 200 for s in fasta(anchor_file).values()))
                self.assertNotIn('--bubble', args)
                self.assertNotIn('--local_contig', args)
                restored = list(fasta(output / f'intermediate_contigs/k{k}.contigs.fa').values())
                self.assertTrue((output / f'intermediate_contigs/k{k}.anchor-graph-paths.tsv').exists())
                self.assertTrue((output / f'intermediate_contigs/k{k}.anchor-threaded.fa').exists())
                for seq in previous:
                    reverse = seq.translate(complement)[::-1]
                    self.assertTrue(any(seq in existing or reverse in existing for existing in restored))
                previous = restored
            self.assertEqual((output / 'intermediate_contigs/k21.bubble_seq.fa').read_text(), '')
            seed = graphs[4]
            self.assertEqual(seed[seed.index('--contig') + 1], str(output / (
                             'tmp/flanking/k21/21.anchors.fa' if anchored_flanking else
                             'intermediate_contigs/k21.contigs.fa')))
            self.assertTrue(fasta(output / 'final.contigs.fa'))
            (output / 'done').unlink()
            (output / 'final.contigs.fa').unlink()
            self.invoke(['--continue', '-o', output])
            self.assertIn('Reuse native checkpoint: restore-reduction-anchors-only-k21',
                          (output / 'log').read_text())

    def test_old_consensus_checkpoint_invalidates_mapping_and_downstream_stages(self):
        output = self.assemble('consensus upgrade')
        checkpoint = output / 'checkpoints.native.json'
        state = json.loads(checkpoint.read_text())
        prefix = 'consensus-local-read-support-primary-580-all-flash-v8-k'
        old_prefix = 'consensus-local-read-support-primary-290-v6-k'
        for key in list(state):
            if prefix in key:
                state[key.replace(prefix, old_prefix)] = state.pop(key)
        checkpoint.write_text(json.dumps(state))
        mapping = output / 'tmp/k21/flank_mapping'
        marker = mapping / 'mapping.json'
        cached = json.loads(marker.read_text())
        cached['configuration']['parameters'] = [param for param in cached['configuration']['parameters']
                                                  if param != '--secondary=no'] + ['--secondary=yes', '--secondary-seq']
        marker.write_text(json.dumps(cached))
        (mapping / 'library.0.sam').write_text('stale secondary-enabled alignment')
        (output / 'done').unlink()
        (output / 'final.contigs.fa').unlink()
        self.invoke(['--continue', '-o', output])
        self.assertTrue((output / 'done').exists())
        self.assertIn('--secondary=no', json.loads(marker.read_text())['configuration']['parameters'])
        self.assertTrue((mapping / 'library.0.sam').read_text().startswith('@'))
        self.assertEqual(json.loads((mapping / 'consensus.json').read_text())['algorithm'],
                         'outward-plurality-read-support-primary-580-all-flash-v8')
        state = json.loads(checkpoint.read_text())
        self.assertTrue(any(prefix in key for key in state))
        self.assertFalse(any(old_prefix in key for key in state))

    def test_native_contig_ends_split_at_580_bases(self):
        output = self.assemble('native ends threshold', extra=['--flanking-consensus-only'])
        rng = random.Random(9183)
        sequences = {n: ''.join(rng.choice('ACGT') for _ in range(n)) for n in (289, 290, 291, 579, 580, 581)}
        source = output / 'intermediate_contigs/k21.contigs.fa'
        source.write_text(''.join('>boundary_{} flag=0 multi=5 len={}\n{}\n'.format(n, n, seq)
                                  for n, seq in sequences.items()))
        Path(str(source) + '.info').write_text('{} {}\n'.format(len(sequences), sum(sequences)))
        checkpoint = output / 'checkpoints.native.json'
        state = json.loads(checkpoint.read_text())
        state = {key: value for key, value in state.items() if 'consensus-local-' not in key}
        checkpoint.write_text(json.dumps(state))
        (output / 'done').unlink()
        (output / 'final.contigs.fa').unlink()
        self.invoke(['--continue', '-o', output])
        actual = {header[1:].split()[0]: seq
                  for header, seq in fasta(output / 'tmp/k21/k21.contig_ends.fa').items()}
        expected = {}
        for n, seq in sequences.items():
            if n < 580:
                expected['boundary_' + str(n)] = seq
            else:
                expected['boundary_' + str(n) + '|start'] = seq[:290]
                expected['boundary_' + str(n) + '|end'] = seq[-290:]
        self.assertEqual(actual, expected)

    def test_independent_phase_lists_and_different_starting_k(self):
        output = self.assemble('independent phase lists', extra=[
            '--reduction-k-list', '23,33,53', '--flanking-k-list', '31,41,61,81'])
        options = json.loads((output / 'options.json').read_text())
        self.assertEqual(options['reduction_k_list'], [23, 33, 53])
        self.assertEqual(options['flanking_k_list'], [31, 41, 61, 81])
        log = (output / 'log').read_text()
        self.assertIn('Reduce edges for k = 23 using k = 33', log)
        self.assertIn('Restarting consensus flanking assembly at k = 31 from k = 23', log)
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
        output = self.assemble('edges-only reduction', extra=[
            '-t', '1', '--reduction-k-list', '21,31,41,51', '--flanking-consensus-only'])
        commands = [json.loads(line[line.index('["--core"'):])
                    for line in (output / 'log').read_text().splitlines()
                    if 'command ' in line and '["--core"' in line]
        graphs = [args for args in commands if args[1] == 'seq2sdbg']
        reductions = [args for args in commands if args[1:4] == ['prefixbf', '--mode', 'reduce']]
        merges = [args for args in commands if args[1:4] == ['prefixbf', '--mode', 'merge']]
        self.assertFalse(merges)
        self.assertFalse(list((output / 'tmp').glob('k*/merged_edges')))
        self.assertEqual(len(reductions), 3)
        for args, k, longer in zip(reductions, (21, 31, 41), (31, 41, 51)):
            source = 'tmp/k21/21' if k == 21 else f'tmp/k{k}/read_count/{k}'
            self.assertEqual(args[args.index('--short_prefix') + 1], str(output / source))
            self.assertEqual(args[args.index('--long_prefix') + 1],
                             str(output / f'tmp/k{longer}/read_count/{longer}'))
            self.assertEqual(args[args.index('--host_mem') + 1], '200000000')
            self.assertIn('--reduce_dedup', args)
            self.assertNotIn('--read_counted_shards', args)
        self.assertLess(commands.index(reductions[-1]), commands.index(graphs[0]))
        for args, k, source in zip(graphs[:4], (51, 41, 31, 21), (0, 51, 41, 31)):
            edges = 'read_count' if source == 0 else 'reduced_edges'
            self.assertEqual(args[args.index('--input_prefix') + 1],
                             str(output / f'tmp/k{k}/{edges}/{k}'))
            self.assertEqual(args[args.index('--kmer_from') + 1], str(source))
            self.assertEqual(args[args.index('-k') + 1], str(k))
            if source:
                self.assertEqual(args[args.index('--contig') + 1],
                                 str(output / f'intermediate_contigs/k{source}.contigs.fa'))
            else:
                self.assertNotIn('--contig', args)
            for flag in ('--bubble', '--addi_contig', '--local_contig'):
                self.assertNotIn(flag, args)
        assembles = [args for args in commands if args[1] == 'assemble']
        self.assertEqual([Path(args[args.index('-o') + 1]).name for args in assembles[:4]],
                         ['k51', 'k41', 'k31', 'k21'])
        self.assertTrue(all('--min_contig_len' not in args for args in assembles[:4]))
        end_of_descent = commands.index(assembles[3])
        self.assertTrue(all(args[1] not in ('minimap2', 'iterate', 'endconsensus')
                            for args in commands[:end_of_descent]))
        self.assertIn(str(output / 'intermediate_contigs/k21.contigs.fa'), graphs[4])
        for k in (31, 41, 51):
            prefix = output / 'tmp/k{}/read_count/{}'.format(k, k)
            self.assertTrue(Path(str(prefix) + '.edges.info').exists())
            self.assertTrue(list(prefix.parent.glob(str(k) + '.edges.[0-9]*')))
            self.assertFalse(Path(str(prefix) + '.sdbg_info').exists())
            self.assertFalse(any(args[args.index('--output_prefix') + 1] == str(prefix)
                                 for args in graphs))
        for k in (21, 31, 41, 51):
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

    def test_old_forward_merge_checkpoints_rebuild_reverse_reduction(self):
        output = self.assemble('reverse reduction upgrade')
        path = output / 'checkpoints.native.json'
        state = json.loads(path.read_text())
        for key in list(state):
            if key.startswith('reduce-counted-v4-k'):
                state[key.replace('reduce-counted-v4-k', 'reduce-dedup-v3-k')] = state.pop(key)
        stale = output / 'tmp/k39/merged_edges'
        stale.mkdir()
        (stale / '39.edges.info').write_text('obsolete; this must not be read')
        path.write_text(json.dumps(state))
        (output / 'done').unlink()
        (output / 'final.contigs.fa').unlink()
        self.invoke(['--continue', '-o', output])
        state = json.loads(path.read_text())
        self.assertIn('assemble-reverse-reduction-k21', state)
        self.assertIn('reverse-reduction-graph-k39-to-k21', state)
        self.assertIn('reduce-counted-v4-k21-k39', state)
        self.assertNotIn('reduce-dedup-v3-k21-k39', state)
        self.assertTrue(fasta(output / 'intermediate_contigs/k21.contigs.fa'))

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
        self.assertTrue((output / 'tmp/k21/k21.contig_ends.extended.fa').exists())
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
        self.assertEqual(collections.Counter(binary_sequences(library)), collections.Counter(self.merged))
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
