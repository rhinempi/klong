"""Prepare a real k21->k29 merge fixture without FLASH or flanking assembly."""
import argparse
import json
from pathlib import Path
import subprocess
import time


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--core', required=True)
    parser.add_argument('--reads', required=True)
    parser.add_argument('--work-dir', required=True)
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--memory-bytes', type=int, default=24000000000)
    args = parser.parse_args()
    root, reads = Path(args.work_dir).resolve(), Path(args.reads).resolve()
    root.mkdir(parents=True, exist_ok=True)
    signature = dict(path=str(reads), size=reads.stat().st_size, mtime_ns=reads.stat().st_mtime_ns)
    manifest = root / 'preparation.json'
    config = dict(core=str(Path(args.core).resolve()), input=signature,
                  threads=args.threads, memory_bytes=args.memory_bytes, min_count=2)
    if manifest.exists():
        state = json.loads(manifest.read_text())
        if state['config'] != config:
            raise RuntimeError('Preparation parameters or input changed; use a fresh directory')
    else:
        state = dict(config=config, stages={})
    visited, resume_valid = {}, True
    def save():
        temporary = manifest.with_suffix('.json.tmp')
        temporary.write_text(json.dumps(state, indent=2) + '\n')
        temporary.replace(manifest)
    def stage(name, arguments, outputs):
        nonlocal resume_valid
        if resume_valid and name in state['stages'] and all(Path(str(path)).exists() for path in outputs):
            visited[name] = state['stages'][name]
            print('Reuse ' + name, flush=True)
            return
        resume_valid = False
        state['stages'] = dict(visited)
        save()
        command = [args.core] + list(map(str, arguments))
        print('Start {}: {}'.format(name, json.dumps(command)), flush=True)
        started = time.monotonic()
        with (root / (name + '.log')).open('w') as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        if not all(Path(str(path)).exists() for path in outputs):
            raise RuntimeError('Missing output from ' + name)
        state['stages'][name] = dict(seconds=time.monotonic() - started, command=command)
        visited[name] = state['stages'][name]
        save()
        print('Finished {} in {:.1f} seconds'.format(name, state['stages'][name]['seconds']), flush=True)
    library = root / 'reads.lib'
    if not library.exists():
        library.write_text(str(reads) + '\ninterleaved ' + str(reads) + '\n')
    stage('library', ['buildlib', library, library, args.threads], [str(library) + '.bin', str(library) + '.lib_info'])
    prefixes = {}
    for k in (21, 29):
        directory = root / ('k' + str(k)) / 'read_count'
        directory.mkdir(parents=True, exist_ok=True)
        prefix = directory / str(k)
        prefixes[k] = prefix
        stage('count-k' + str(k), ['count', '-k', k, '-m', 2, '--host_mem', args.memory_bytes,
                                  '--mem_flag', 1, '--output_prefix', prefix, '--num_cpu_threads', args.threads,
                                  '--read_lib_file', library], [str(prefix) + '.edges.info'])
    directory = root / 'k21/reduced_edges'
    directory.mkdir(parents=True, exist_ok=True)
    reduced = directory / '21'
    stage('reduce-k21', ['prefixbf', '--mode', 'reduce', '--short_prefix', prefixes[21],
                         '--long_prefix', prefixes[29], '--output_prefix', reduced,
                         '--short_k', 21, '--long_k', 29, '-t', args.threads, '--human_readable', 0,
                         '--max_bloom_gb', 60, '--reserve_gb', 1, '--host_mem', args.memory_bytes],
          [str(reduced) + '.edges.info'])
    directory = root / 'k21/assembly_graph'
    directory.mkdir(parents=True, exist_ok=True)
    graph = directory / '21'
    stage('graph-k21', ['seq2sdbg', '--host_mem', args.memory_bytes, '--mem_flag', 1,
                        '--output_prefix', graph, '--num_cpu_threads', args.threads, '-k', 21,
                        '--kmer_from', 0, '--input_prefix', reduced], [str(graph) + '.sdbg_info'])
    directory = root / 'intermediate_contigs'
    directory.mkdir(parents=True, exist_ok=True)
    contigs = directory / 'k21'
    stage('assemble-k21', ['assemble', '-s', graph, '-o', contigs, '-t', args.threads,
                           '--min_standalone', 0, '--prune_level', 0, '--merge_len', 0,
                           '--merge_similar', 0.95, '--cleaning_rounds', 1, '--disconnect_ratio', 0.1,
                           '--low_local_ratio', 0.2, '--min_depth', 0, '--bubble_level', 0,
                           '--max_tip_len', 0], [str(contigs) + '.contigs.fa'])
    short_path = Path(str(contigs) + '.lt_k29_plus1.contigs.fa')
    long_path = Path(str(contigs) + '.ge_k29_plus1.contigs.fa')
    counts = {'short': 0, 'long': 0}
    with Path(str(contigs) + '.contigs.fa').open() as source, short_path.open('w') as short, long_path.open('w') as long:
        header, parts = None, []
        def emit():
            if header is None:
                return
            sequence = ''.join(parts)
            small = len(sequence) < 30
            (short if small else long).write(header + '\n' + sequence + '\n')
            counts['short' if small else 'long'] += 1
        for line in source:
            line = line.rstrip('\r\n')
            if line.startswith('>'):
                emit()
                header, parts = line, []
            else:
                parts.append(line)
        emit()
    state['merge_inputs'] = dict(long_prefix=str(prefixes[29]), contig=str(short_path),
                                 long_contig=str(long_path), contig_counts=counts)
    save()
    print(json.dumps(state['merge_inputs'], indent=2), flush=True)


if __name__ == '__main__':
    main()
