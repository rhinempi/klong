"""The production CMake pigz generator must work without Python subprocesses."""
from pathlib import Path
import os
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(shutil.which('cmake'), 'cmake unavailable')
class NativePigzBuildTest(unittest.TestCase):
    @unittest.skipUnless(os.environ.get('MEGAHIT_TEST_LEGACY_CMAKE'), 'legacy CMake not configured')
    def test_native_pigz_generator_with_cmake_310(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / 'pigz.cpp'
            subprocess.run([os.environ['MEGAHIT_TEST_LEGACY_CMAKE'], '-DENABLE_PIGZ=ON',
                            '-DPIGZ_SOURCE=' + str(ROOT / 'src/pigz'),
                            '-DPIGZ_BUILD_DIR=' + str(root / 'build'),
                            '-DOUTPUT=' + str(output), '-P', str(ROOT / 'cmake/embed_native_pigz.cmake')],
                           capture_output=True, check=True)
            self.assertNotIn('bundled_pigz_size = 0', output.read_text())
            self.assertTrue((root / 'build/pigz').exists())

    def test_optional_compiler_failure_generates_native_zlib_fallback(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            binary = root / 'pigz'
            binary.write_bytes(b'stale executable')
            output = root / 'pigz.cpp'
            result = subprocess.run(['cmake', '-DENABLE_PIGZ=ON',
                                     '-DPIGZ_SOURCE=' + str(ROOT / 'src/pigz'),
                                     '-DPIGZ_BUILD_DIR=' + str(root), '-DC_COMPILER=/missing/compiler',
                                     '-DOUTPUT=' + str(output), '-P', str(ROOT / 'cmake/embed_native_pigz.cmake')],
                                    capture_output=True, text=True, check=True)
            self.assertIn('native zlib gzip fallback', result.stderr)
            self.assertIn('bundled_pigz_size = 0', output.read_text())
            self.assertFalse(binary.exists())

    def test_disabled_compressor_does_not_need_source_or_python(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'pigz.cpp'
            subprocess.run(['cmake', '-DENABLE_PIGZ=OFF', '-DOUTPUT=' + str(output),
                            '-P', str(ROOT / 'cmake/embed_native_pigz.cmake')],
                           capture_output=True, check=True)
            self.assertIn('bundled_pigz_size = 0', output.read_text())
