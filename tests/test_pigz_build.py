"""Bundled pigz packaging must be optional and independent of system PATH."""
import contextlib
import gzip
import importlib.machinery
import importlib.util
import io
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

from test_flash import load_launcher

ROOT = Path(__file__).resolve().parents[1]
LEGACY_CMAKE = os.environ.get('MEGAHIT_TEST_LEGACY_CMAKE')


def load_file(name, path):
    loader = importlib.machinery.SourceFileLoader(name, str(path))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


class OptionalPigzFailureTest(unittest.TestCase):
    def test_legacy_configure_form_and_absolute_paths(self):
        helper = load_file('optional_pigz', ROOT / 'cmake/build_optional_pigz.py')
        with tempfile.TemporaryDirectory(prefix='legacy pigz ') as directory:
            build_dir = Path(directory) / 'build directory'
            relative_source = os.path.relpath(ROOT / 'src/pigz')
            relative_build = os.path.relpath(build_dir)
            with mock.patch.object(helper.subprocess, 'run') as run:
                result = helper.build_pigz(relative_source, relative_build)
            self.assertEqual(result, build_dir / ('pigz.exe' if os.name == 'nt' else 'pigz'))
            configure, build, verify = run.call_args_list
            self.assertEqual(configure.args[0][1], str(ROOT / 'src/pigz'))
            self.assertNotIn('-S', configure.args[0])
            self.assertNotIn('-B', configure.args[0])
            self.assertEqual(configure.kwargs['cwd'], str(build_dir))
            self.assertEqual(build.args[0][:3], ['cmake', '--build', '.'])
            self.assertEqual(build.kwargs['cwd'], str(build_dir))
            self.assertEqual(verify.args[0], [str(result), '--version'])

    @unittest.skipUnless(LEGACY_CMAKE, 'set MEGAHIT_TEST_LEGACY_CMAKE to test an older CMake')
    def test_real_pre_313_cmake_packages_pigz(self):
        with tempfile.TemporaryDirectory(prefix='legacy pigz ') as directory:
            root = Path(directory)
            generated = self.make_package(root, [
                '--pigz-source', str(ROOT / 'src/pigz'), '--pigz-build-dir', str(root / 'build'),
                '--pigz-cmake', LEGACY_CMAKE])
            self.assertTrue(generated.PIGZ_B64, generated.build_output)
            launcher, runtime_bin = self.prepare_launcher(generated, root)
            version = subprocess.check_output([str(runtime_bin / 'pigz'), '--version'],
                                              stderr=subprocess.STDOUT)
            self.assertIn(b'pigz', version)
            source, output = root / 'input.fastq', root / 'output.fa.gz'
            source.write_bytes(b'@read\nACGT\n+\nIIII\n')
            with mock.patch.dict(os.environ, {'PATH': ''}):
                launcher.fastq_to_fasta_gz(str(source), str(output), num_threads=2)
            self.assertEqual(gzip.decompress(output.read_bytes()), b'>read\nACGT\n')

    def test_failed_build_removes_stale_binary(self):
        helper = load_file('optional_pigz', ROOT / 'cmake/build_optional_pigz.py')
        with tempfile.TemporaryDirectory() as directory:
            binary = Path(directory) / ('pigz.exe' if os.name == 'nt' else 'pigz')
            binary.write_bytes(b'stale executable')
            with mock.patch.object(helper.subprocess, 'run', side_effect=OSError('compiler failed')), \
                    contextlib.redirect_stderr(io.StringIO()) as errors:
                result = helper.build_pigz(ROOT / 'src/pigz', directory)
            self.assertIsNone(result)
            self.assertFalse(binary.exists())
            self.assertIn('Python gzip fallback', errors.getvalue())

    def test_optional_unpack_failure_does_not_abort(self):
        with tempfile.TemporaryDirectory() as directory:
            generated = self.make_package(Path(directory))
            generated.PIGZ_B64 = 'present'
            with mock.patch.object(generated.tempfile, 'mkstemp', side_effect=OSError('read-only')), \
                    contextlib.redirect_stderr(io.StringIO()) as errors:
                generated._extract_optional_pigz(directory)
            self.assertIn('use Python gzip', errors.getvalue())

    @staticmethod
    def make_package(root, extra=()):
        core = root / 'core'
        core.write_bytes(b'fixture core, not executed')
        executable = root / 'megahit'
        result = subprocess.run([sys.executable, str(ROOT / 'cmake/embed_megahit.py'),
                                 '--launcher', str(ROOT / 'src/megahit'), '--core', str(core),
                                 '--output', str(executable)] + list(extra),
                                capture_output=True, text=True, check=True)
        generated = load_file('packaged_megahit', executable)
        generated.build_output = result.stdout + result.stderr
        return generated

    def prepare_launcher(self, generated, root):
        with mock.patch.dict(os.environ, {'MEGAHIT_EMBED_CACHE': str(root / 'cache')}), \
                mock.patch.object(generated.sys, 'argv', [str(root / 'megahit')]):
            runtime_bin = Path(generated._prepare_runtime())
        launcher = load_launcher()
        launcher.__file__ = str(runtime_bin / 'megahit')
        return launcher, runtime_bin

    def check_fallback(self, generated, root):
        self.assertFalse(generated.PIGZ_B64)
        launcher, runtime_bin = self.prepare_launcher(generated, root)
        self.assertFalse((runtime_bin / 'pigz').exists())
        source, output = root / 'input.fastq', root / 'output.fa.gz'
        source.write_bytes(b'@read\nACGT\n+\nIIII\n')
        with mock.patch.object(launcher.subprocess, 'Popen') as popen:
            launcher.fastq_to_fasta_gz(str(source), str(output), num_threads=4)
        popen.assert_not_called()
        self.assertEqual(gzip.decompress(output.read_bytes()), b'>read\nACGT\n')

    def test_disabled_bundle_packages_gzip_fallback(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            generated = self.make_package(root)
            self.check_fallback(generated, root)

    @unittest.skipUnless(shutil.which('cmake'), 'cmake is not available')
    def test_real_build_failure_still_packages_working_fallback(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            generated = self.make_package(root, [
                '--pigz-source', str(ROOT / 'src/pigz'), '--pigz-build-dir', str(root / 'pigz-build'),
                '--pigz-cc', str(root / 'missing-compiler')])
            self.assertIn('optional pigz build failed', generated.build_output)
            self.check_fallback(generated, root)

    @unittest.skipUnless(shutil.which('cmake') and shutil.which('cc'), 'C build tools unavailable')
    def test_real_bundle_compresses_without_system_path(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            generated = self.make_package(root, [
                '--pigz-source', str(ROOT / 'src/pigz'), '--pigz-build-dir', str(root / 'pigz-build')])
            self.assertTrue(generated.PIGZ_B64, generated.build_output)
            launcher, runtime_bin = self.prepare_launcher(generated, root)
            self.assertEqual(launcher.find_bundled_pigz(), str(runtime_bin / 'pigz'))
            source, output = root / 'input.fastq', root / 'output.fa.gz'
            source.write_bytes(b'@read\nACGT\n+\nIIII\n')
            with mock.patch.dict(os.environ, {'PATH': ''}):
                launcher.fastq_to_fasta_gz(str(source), str(output), num_threads=2)
            self.assertEqual(gzip.decompress(output.read_bytes()), b'>read\nACGT\n')
            # Repair the optional executable independently of the core cache marker.
            (runtime_bin / 'pigz').unlink()
            self.prepare_launcher(generated, root)
            self.assertTrue((runtime_bin / 'pigz').is_file())
