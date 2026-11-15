"""Measure complete native merge runs, including exact-sort temporary I/O."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import statistics
import subprocess


def fingerprint(prefix):
    hashes = {}
    for path in sorted(prefix.parent.glob(prefix.name + '.edges.*')):
        digest = hashlib.sha256()
        with path.open('rb') as source:
            for block in iter(lambda: source.read(1 << 20), b''):
                digest.update(block)
        hashes[path.name[len(prefix.name):]] = digest.hexdigest()
    return hashes


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--core', required=True)
    parser.add_argument('--work-dir', required=True)
    parser.add_argument('--long-prefix', required=True)
    parser.add_argument('--contig', required=True)
    parser.add_argument('--long-contig')
    parser.add_argument('--short-k', type=int, default=21)
    parser.add_argument('--long-k', type=int, default=29)
    parser.add_argument('--threads', default='1,8')
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--sort-memory', type=float, default=128)
    parser.add_argument('--memory-bytes', type=int, default=4000000000)
    args = parser.parse_args()
    core_sha256 = hashlib.sha256(Path(args.core).read_bytes()).hexdigest()
    root = Path(args.work_dir)
    root.mkdir(parents=True, exist_ok=True)
    rows, references = [], {}
    for threads in map(int, args.threads.split(',')):
        for repeat in range(1, args.repeats + 1):
            for backend in ('bloom', 'sort'):
                label = '{}-t{}-r{}'.format(backend, threads, repeat)
                prefix, timing = root / label, root / (label + '.time')
                command = ['/usr/bin/time', '-f', '%e %U %S %M %I %O', '-o', str(timing), args.core,
                           'prefixbf', '--mode', 'merge', '--long_prefix', args.long_prefix,
                           '--contig', args.contig, '--short_k', str(args.short_k), '--long_k', str(args.long_k),
                           '--output_prefix', str(prefix), '-t', str(threads), '--human_readable', '0',
                           '--merge_index', 'hash', '--merge_dedup', backend, '--sort_memory_mb', str(args.sort_memory),
                           '--host_mem', str(args.memory_bytes)]
                if args.long_contig:
                    command += ['--long_contig', args.long_contig]
                print('Start ' + label, flush=True)
                log_path = root / (label + '.log')
                with log_path.open('w') as log:
                    subprocess.run(command, stdout=subprocess.PIPE, stderr=log, text=True, check=True)
                stderr = log_path.read_text()
                elapsed, user, system, rss, reads, writes = timing.read_text().split()
                hashes = fingerprint(prefix)
                if backend == 'bloom' and threads not in references:
                    references[threads] = hashes
                row = dict(backend=backend, threads=threads, repeat=repeat, seconds=float(elapsed),
                           user_seconds=float(user), system_seconds=float(system), peak_mib=int(rss) / 1024,
                           fs_input_blocks=int(reads), fs_output_blocks=int(writes),
                           matches_bloom=hashes == references[threads], output_sha256=hashes)
                for line in Path(str(prefix) + '.edges.info').read_text().splitlines():
                    fields = line.split()
                    if fields and fields[0] == 'num_edges':
                        row['output_edges'] = int(fields[1])
                        break
                match = re.search(r'(\d+) candidates, (\d+) new edges, (\d+) duplicate candidates, '
                                  r'(\d+) runs, (\d+) temporary bytes written, (\d+) temporary bytes read', stderr)
                if match:
                    row.update(zip(('candidates', 'new_edges', 'duplicates', 'runs',
                                    'temporary_bytes_written', 'temporary_bytes_read'), map(int, match.groups())))
                if list(root.glob(prefix.name + '.sort-*')):
                    raise RuntimeError('Sort workspace was not cleaned up')
                rows.append(row)
                temporary = root / 'results.partial.json.tmp'
                temporary.write_text(json.dumps(dict(inputs=vars(args), core_sha256=core_sha256, runs=rows), indent=2) + '\n')
                temporary.replace(root / 'results.partial.json')
                print(json.dumps({key: value for key, value in row.items() if key != 'output_sha256'}), flush=True)
                for path in root.glob(prefix.name + '.edges.*'):
                    path.unlink()
    summary = []
    for threads in sorted({row['threads'] for row in rows}):
        for backend in ('bloom', 'sort'):
            selected = [row for row in rows if row['threads'] == threads and row['backend'] == backend]
            summary.append(dict(backend=backend, threads=threads,
                                seconds=statistics.median(row['seconds'] for row in selected),
                                peak_mib=statistics.median(row['peak_mib'] for row in selected),
                                all_outputs_match=all(row['matches_bloom'] for row in selected),
                                temporary_bytes_written=statistics.median(row.get('temporary_bytes_written', 0)
                                                                          for row in selected)))
    (root / 'results.json').write_text(json.dumps(dict(inputs=vars(args), core_sha256=core_sha256,
                                                     runs=rows, summary=summary), indent=2) + '\n')
    print(json.dumps(dict(summary=summary), indent=2))


if __name__ == '__main__':
    main()
