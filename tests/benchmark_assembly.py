"""Before/after assembly timing with exact, order-independent output comparison."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import statistics
import subprocess
import time

from test_assembly_optimizations import records, semantic_record


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def fingerprint(path, work, threads):
    unsorted, ordered = work / 'records.unsorted', work / 'records.sorted'
    count, bases = 0, 0
    with unsorted.open('w') as output:
        for header, sequence in records(path):
            output.write(semantic_record(header, sequence) + '\n')
            count += 1
            bases += len(sequence)
    subprocess.run(['sort', '-S', '128M', '--parallel=' + str(threads), '-T', str(work),
                    '-o', str(ordered), str(unsorted)], env=dict(os.environ, LC_ALL='C'), check=True)
    result = {'records': count, 'bases': bases, 'semantic_sha256': digest(ordered),
              'info': Path(str(path) + '.info').read_text().strip()}
    stable = hashlib.sha256()
    loop_rows, coverage = [], []
    with ordered.open('rb') as source:
        for line in source:
            flag, multi, remainder = line.split(b'\t', 2)
            if int(flag) & 2:
                loop_rows.append(flag + b'\t*\t' + remainder)
                coverage.append([remainder.decode().strip(), multi.decode()])
            else:
                stable.update(line)
    for line in sorted(loop_rows):
        stable.update(line)
    result['sequence_and_nonloop_metadata_sha256'] = stable.hexdigest()
    result['loop_coverage'] = sorted(coverage)
    unsorted.unlink()
    ordered.unlink()
    return result


def stable_outputs(outputs):
    fields = ('records', 'bases', 'info', 'sequence_and_nonloop_metadata_sha256')
    return {name: {key: data[key] for key in fields} for name, data in outputs.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--before', required=True, type=Path)
    parser.add_argument('--after', required=True, type=Path)
    parser.add_argument('--graph', required=True, type=Path)
    parser.add_argument('--work-dir', required=True, type=Path)
    parser.add_argument('--profile', choices=['reduction', 'cleaning'], default='reduction')
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--repeats', type=int, default=1)
    parser.add_argument('--allow-loop-coverage-variation', action='store_true',
                        help='permit only start-edge-dependent circular coverage differences; still report them')
    parser.add_argument('--recheck', action='store_true', help='verify saved outputs without rerunning assembly')
    args = parser.parse_args()
    if args.threads < 1 or args.repeats < 1:
        parser.error('threads and repeats must be positive')
    if not args.recheck:
        args.work_dir.mkdir(parents=True, exist_ok=False)
    shared = ['assemble', '-s', str(args.graph.resolve()), '-t', str(args.threads),
              '--disconnect_ratio', '0.1', '--low_local_ratio', '0.2', '--merge_similar', '0.95']
    if args.profile == 'reduction':
        shared += ['--min_standalone', '0', '--prune_level', '0', '--merge_len', '0',
                   '--cleaning_rounds', '1', '--bubble_level', '0', '--max_tip_len', '0', '--min_depth', '0']
    else:
        shared += ['--min_standalone', '300', '--prune_level', '2', '--merge_len', '20',
                   '--cleaning_rounds', '5', '--bubble_level', '2', '--max_tip_len', '-1',
                   '--min_depth', '0', '--careful_bubble']
    report = {'graph': str(args.graph.resolve()), 'profile': args.profile,
              'threads': args.threads, 'repeats': args.repeats,
              'binary_sha256': {label: digest(binary) for label, binary in
                                [('before', args.before), ('after', args.after)]}, 'runs': []}
    if args.recheck:
        saved = json.loads((args.work_dir / 'results.json').read_text())
        for key in ('graph', 'profile', 'threads', 'repeats', 'binary_sha256'):
            if saved[key] != report[key]:
                raise RuntimeError('Recheck configuration differs: ' + key)
        report = saved
    report['allow_loop_coverage_variation'] = args.allow_loop_coverage_variation
    reference = None
    for repeat in range(args.repeats):
        runs = [('before', args.before), ('after', args.after)]
        if repeat % 2:
            runs.reverse()
        for label, binary in runs:
            name = '{}-{}'.format(label, repeat + 1)
            output = args.work_dir / name
            command = [str(binary.resolve())] + shared + ['-o', str(output.resolve())]
            if args.recheck:
                row = next(row for row in report['runs'] if row['label'] == label and row['repeat'] == repeat + 1)
                elapsed = row['seconds']
            else:
                print('Starting ' + name, flush=True)
                started = time.monotonic()
                with (args.work_dir / (name + '.log')).open('w') as log:
                    subprocess.run(['/usr/bin/time', '-v', '-o', str(args.work_dir / (name + '.time'))]
                                   + command, stdout=log, stderr=subprocess.STDOUT, check=True)
                elapsed = time.monotonic() - started
            print('{}: {:.3f}s; checking outputs'.format(name, elapsed), flush=True)
            outputs = {}
            for suffix in ('contigs.fa', 'final.contigs.fa', 'bubble_seq.fa', 'addi.fa'):
                path = Path(str(output) + '.' + suffix)
                if path.exists():
                    outputs[suffix] = fingerprint(path, args.work_dir, args.threads)
            if reference is None:
                reference = outputs
            timing = (args.work_dir / (name + '.time')).read_text()
            log = (args.work_dir / (name + '.log')).read_text()
            row = {'label': label, 'repeat': repeat + 1, 'seconds': elapsed, 'command': command,
                   'max_rss_kib': int(re.search(r'Maximum resident set size \(kbytes\): (\d+)', timing)[1]),
                   'outputs': outputs, 'matches_reference': outputs == reference,
                   'sequence_and_nonloop_metadata_match': stable_outputs(outputs) == stable_outputs(reference),
                   'statistics': re.findall(r'Max: .*total size: \d+,', log)}
            if args.recheck:
                index = next(i for i, old in enumerate(report['runs'])
                             if old['label'] == label and old['repeat'] == repeat + 1)
                report['runs'][index] = row
            else:
                report['runs'].append(row)
            (args.work_dir / 'results.json').write_text(json.dumps(report, indent=2) + '\n')
            if outputs != reference and not (args.allow_loop_coverage_variation and row['sequence_and_nonloop_metadata_match']):
                raise RuntimeError('Assembly output differs: ' + name)
            print(name + (': exact normalized output match' if outputs == reference else
                          ': sequences/nonloop metadata match; circular coverage differs (see results.json)'), flush=True)
    report['median_seconds'] = {label: statistics.median(row['seconds'] for row in report['runs']
                                                       if row['label'] == label)
                                for label in ('before', 'after')}
    report['speedup'] = report['median_seconds']['before'] / report['median_seconds']['after']
    (args.work_dir / 'results.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'median_seconds': report['median_seconds'], 'speedup': report['speedup']}), flush=True)


if __name__ == '__main__':
    main()
