"""Benchmark serial and threaded consensus with byte-exact output checks."""
import argparse
import hashlib
import json
from pathlib import Path
import random
import statistics
import subprocess
import time


def sha256(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def synthetic(directory, reads, references):
    rng = random.Random(12701)
    sequences = [''.join(rng.choice('ACGT') for _ in range(290)) for _ in range(references)]
    names = ['ref' + str(i) + ('|start' if i % 3 == 0 else '|end' if i % 3 == 1 else '')
             for i in range(references)]
    fasta, sam = directory / 'ends.fa', directory / 'reads.sam'
    with fasta.open('w') as output:
        for name, sequence in zip(names, sequences):
            output.write('>' + name + '\n' + sequence + '\n')
    with sam.open('w', buffering=1024 * 1024) as output:
        for name in names:
            output.write('@SQ\tSN:' + name + '\tLN:290\n')
        for i in range(reads):
            ref = i % references
            cigar, sequence = '290M', sequences[ref]
            if i % 4:
                cigar, sequence = '10S290M40S', 'ACGTACGTAC' + sequence + ('G' if i % 9 == 0 else 'A') * 40
            row = 'q{}\t0\t{}\t1\t60\t{}\t*\t0\t0\t{}\t*\n'.format(i, names[ref], cigar, sequence)
            output.write(row)
            if i % 100 == 0:
                output.write(row)
    return fasta, [sam]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--before', required=True, type=Path)
    parser.add_argument('--after', required=True, type=Path)
    parser.add_argument('--work-dir', required=True, type=Path)
    parser.add_argument('--reference', type=Path)
    parser.add_argument('--sam', nargs='+', type=Path)
    parser.add_argument('--synthetic-reads', type=int)
    parser.add_argument('--synthetic-references', type=int, default=4096)
    parser.add_argument('--threads', default='1,2,4,8')
    parser.add_argument('--repeats', type=int, default=3)
    args = parser.parse_args()
    threads = [int(n) for n in args.threads.split(',')]
    if not threads or min(threads) < 1 or args.repeats < 1:
        parser.error('threads and repeats must be positive')
    if args.synthetic_reads is not None:
        if args.synthetic_reads < 1 or args.synthetic_references < 1 or args.reference or args.sam:
            parser.error('synthetic counts must be positive, without reference/SAM arguments')
    elif not args.reference or not args.sam:
        parser.error('supply either synthetic reads or a reference and SAM files')
    args.work_dir.mkdir(parents=True, exist_ok=False)
    reference, sams = args.reference, args.sam
    if args.synthetic_reads is not None:
        print('Generating synthetic SAM', flush=True)
        reference, sams = synthetic(args.work_dir, args.synthetic_reads, args.synthetic_references)
    report = {'reference': str(reference.resolve()), 'sam': [str(path.resolve()) for path in sams],
              'sam_bytes': sum(path.stat().st_size for path in sams),
              'synthetic_reads': args.synthetic_reads, 'repeats': args.repeats,
              'binary_sha256': {'before': sha256(args.before), 'after': sha256(args.after)}, 'runs': []}
    expected = None
    configurations = [('before', args.before, None)] + [('after-t' + str(n), args.after, n) for n in threads]
    for repeat in range(args.repeats):
        order = configurations if repeat % 2 == 0 else list(reversed(configurations))
        for label, binary, workers in order:
            name = '{}-{}'.format(label, repeat + 1)
            output, stats = args.work_dir / (name + '.fa'), args.work_dir / (name + '.tsv')
            command = [str(binary.resolve()), 'endconsensus']
            if workers is not None:
                command += ['-t', str(workers)]
            command += [str(reference.resolve()), str(output.resolve()), str(stats.resolve())]
            command += [str(path.resolve()) for path in sams]
            measurement = args.work_dir / (name + '.time')
            started = time.monotonic()
            with (args.work_dir / (name + '.log')).open('w') as log:
                subprocess.run(['/usr/bin/time', '-f', '%e %U %S %M', '-o', str(measurement)] + command,
                               stdout=log, stderr=subprocess.STDOUT, check=True)
            seconds = time.monotonic() - started
            outputs = {'fasta': sha256(output), 'info': sha256(str(output) + '.info'), 'tsv': sha256(stats)}
            if expected is None:
                expected = outputs
            timing = measurement.read_text().split()
            row = {'label': label, 'repeat': repeat + 1, 'seconds': seconds,
                   'user_seconds': float(timing[1]), 'system_seconds': float(timing[2]),
                   'rss_kib': int(timing[3]), 'command': command, 'outputs': outputs,
                   'matches_baseline': outputs == expected}
            report['runs'].append(row)
            (args.work_dir / 'results.json').write_text(json.dumps(report, indent=2) + '\n')
            if outputs != expected:
                raise RuntimeError('Output mismatch: ' + name)
            print('{}: {:.4f}s, exact output match'.format(name, seconds), flush=True)
    report['median_seconds'] = {label: statistics.median(row['seconds'] for row in report['runs']
                                                       if row['label'] == label) for label, _, _ in configurations}
    report['speedup_vs_before'] = {label: report['median_seconds']['before'] / seconds
                                  for label, seconds in report['median_seconds'].items()}
    (args.work_dir / 'results.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'median_seconds': report['median_seconds'],
                     'speedup_vs_before': report['speedup_vs_before']}), flush=True)


if __name__ == '__main__':
    main()
