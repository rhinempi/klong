"""Build pigz without making an optional compressor failure fatal to MEGAHIT."""
import os
from pathlib import Path
import subprocess
import sys


def build_pigz(source, build_dir, cmake='cmake', compiler=None,
               zlib_include=None, zlib_library=None, static=False):
    build_dir = Path(build_dir)
    binary = build_dir / ('pigz.exe' if os.name == 'nt' else 'pigz')
    try:
        build_dir.mkdir(parents=True, exist_ok=True)
        # A failed rebuild must not package an older executable left behind.
        if binary.exists():
            binary.unlink()
        configure = [cmake, '-S', str(source), '-B', str(build_dir),
                     '-DCMAKE_BUILD_TYPE=Release', '-DSTATIC_BUILD=' + ('ON' if static else 'OFF')]
        for name, value in [('CMAKE_C_COMPILER', compiler), ('ZLIB_INCLUDE_DIR', zlib_include),
                            ('ZLIB_LIBRARY', zlib_library)]:
            if value:
                configure.append('-D{}={}'.format(name, value))
        commands = [configure, [cmake, '--build', str(build_dir), '--config', 'Release',
                                '--target', 'pigz']]
        for command in commands:
            subprocess.run(command, check=True, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, timeout=180)
        subprocess.run([str(binary), '--version'], check=True, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, timeout=10)
        print('Bundled pigz built: {}'.format(binary))
        return binary
    except (OSError, subprocess.SubprocessError) as error:
        if binary.exists():
            try:
                binary.unlink()
            except OSError:
                pass
        details = getattr(error, 'stdout', None) or b''
        if isinstance(details, bytes):
            details = details.decode('utf-8', errors='replace')
        print('WARNING: optional pigz build failed; packaging Python gzip fallback.\n{}\n{}'.format(
            error, details[-2000:]), file=sys.stderr)
        return None
