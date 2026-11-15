"""Packed flank extensions across the 32-base word boundary."""
import os
from pathlib import Path
import random
import struct
import subprocess
import tempfile
import unittest

from test_query_merge import canonical, reverse_complement

CORE = os.environ.get('MEGAHIT_TEST_CORE')
BASELINE = os.environ.get('MEGAHIT_TEST_BASELINE_CORE')


@unittest.skipUnless(CORE, 'set MEGAHIT_TEST_CORE to the built core')
class IterateStepTest(unittest.TestCase):
    def test_packed_extensions_and_reverse_complements(self):
        with tempfile.TemporaryDirectory(prefix='iterate-steps-') as directory:
            root = Path(directory)
            for step in (20, 28, 40):
                rng = random.Random(371 + step)
                dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
                k, gap = 101, 2 * step - 10
                sequence = dna(k + 1 + gap)
                contigs = [sequence[:k + step] + dna(70), dna(70) + sequence[gap - step + 1:]]
                contig = root / 'contigs.fa'
                contig.write_text(''.join('>c{} flag=0 multi=7\n{}\n'.format(i, seq)
                                          for i, seq in enumerate(contigs)))
                bubble = root / 'bubble.fa'
                bubble.write_text('')
                expected = {canonical(sequence[i:i + k + step + 1])
                            for i in range(gap - step + 1)}
                for reverse in (False, True):
                    read = root / 'read.fa'
                    read.write_text('>r\n' + (reverse_complement(sequence) if reverse else sequence) + '\n')
                    descriptor = root / 'library'
                    descriptor.write_text(str(read) + '\nse ' + str(read) + '\n')
                    subprocess.run([CORE, 'buildlib', str(descriptor), str(descriptor), '2'],
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, timeout=30)
                    output = root / 'iterated'
                    args = ['iterate', '-c', str(contig), '-b', str(bubble),
                            '-r', str(descriptor) + '.bin', '-k', str(k), '-s', str(step),
                            '-t', '2', '-o', str(output)]
                    subprocess.run([CORE] + args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                   check=True, timeout=30)
                    data = Path(str(output) + '.edges.0').read_bytes()
                    length = k + step + 1
                    words = (length * 2 + 16 + 31) // 32
                    self.assertEqual(len(data) % (words * 4), 0)
                    actual = set()
                    for record in struct.iter_unpack('=' + 'I' * words, data):
                        # Iterative edge records use the graph's reversed base order.
                        decoded = ''.join('ACGT'[(record[j // 16] >> (30 - 2 * (j % 16))) & 3]
                                          for j in range(length))
                        actual.add(canonical(decoded[::-1]))
                    self.assertEqual(actual, expected, (step, reverse))
                    if BASELINE and step <= 28:
                        subprocess.run([BASELINE] + args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                       check=True, timeout=30)
                        self.assertEqual(Path(str(output) + '.edges.0').read_bytes(), data)
