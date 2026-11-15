"""Optional Python harness; the measured index implementations are native."""
import argparse
import json
from pathlib import Path
import random
import statistics
import subprocess

from test_query_merge import canonical, pack


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--core', required=True)
    parser.add_argument('--work-dir', required=True)
    parser.add_argument('--long-prefix')
    parser.add_argument('--contig')
    parser.add_argument('--short-k', type=int, default=21)
    parser.add_argument('--long-k', type=int, default=29)
    parser.add_argument('--queries', type=int, default=100000)
    parser.add_argument('--edges', type=int, default=500000)
    parser.add_argument('--threads', default='1,8')
    parser.add_argument('--repeats', type=int, default=3)
    args = parser.parse_args()
    root = Path(args.work_dir)
    root.mkdir(parents=True, exist_ok=True)
    if bool(args.long_prefix) != bool(args.contig):
        parser.error('--long-prefix and --contig must be supplied together')
    if args.queries < 1 or args.edges < 1 or args.repeats < 1 or args.short_k < 1 or args.long_k <= args.short_k:
        parser.error('positive sizes and short-k < long-k are required')
    if args.long_prefix:
        prefix, contig = args.long_prefix, args.contig
    else:
        rng = random.Random(71453)
        dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
        length = args.long_k + 1
        sequences = [dna(args.short_k + min(5, args.long_k - args.short_k))
                     for _ in range(args.queries)]
        contig, prefix = root / 'contigs.fa', root / 'reads'
        with contig.open('w') as output:
            for i, sequence in enumerate(sequences):
                output.write('>q{}\n{}\n'.format(i, sequence[::-1]))
        with Path(str(prefix) + '.edges.0').open('wb') as output:
            for i in range(args.edges):
                sequence = (dna(length - args.short_k) + sequences[i % len(sequences)][:args.short_k]
                            if i % 4 == 0 else dna(length))
                output.write(pack(canonical(sequence), 7))
        words = (2 * length + 16 + 31) // 32
        Path(str(prefix) + '.edges.info').write_text(
            'kmer_size {}\nwords_per_edge {}\nnum_files 1\nnum_buckets 0\n'
            'num_edges {}\nis_sorted 0\n'.format(args.long_k, words, args.edges))
    runs, expected = [], None
    for threads in map(int, args.threads.split(',')):
        for repeat in range(args.repeats):
            for backend in ('hash', 'fm'):
                memory = root / 'peak-rss.txt'
                command = ['/usr/bin/time', '-f', '%M', '-o', str(memory), args.core,
                           'prefixbf', '--mode', 'index-bench', '--long_prefix', str(prefix),
                           '--contig', str(contig), '--short_k', str(args.short_k),
                           '--long_k', str(args.long_k), '-t', str(threads), '--merge_index', backend]
                result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        text=True, check=True)
                row = json.loads(result.stdout)
                signature = {key: row[key] for key in ('records', 'orientations', 'prefix_hits',
                                                       'suffix_hits', 'prefix_keys', 'suffix_keys')}
                if expected is not None and signature != expected:
                    raise RuntimeError('Index backends disagree: {}'.format(row))
                expected = signature
                row.update(repeat=repeat + 1, peak_mib=int(memory.read_text().strip()) / 1024)
                runs.append(row)
                print(json.dumps(row), flush=True)
    summary = []
    for threads in sorted({row['threads'] for row in runs}):
        for backend in ('hash', 'fm'):
            selected = [row for row in runs if row['threads'] == threads and row['backend'] == backend]
            summary.append(dict(backend=backend, threads=threads,
                                build_seconds=statistics.median(row['build_seconds'] for row in selected),
                                lookup_seconds=statistics.median(row['lookup_seconds'] for row in selected),
                                peak_mib=statistics.median(row['peak_mib'] for row in selected)))
    report = dict(inputs=dict(long_prefix=str(prefix), contig=str(contig)), counts=expected,
                  runs=runs, summary=summary)
    (root / 'results.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(dict(summary=summary), indent=2))


if __name__ == '__main__':
    main()
