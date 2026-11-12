#!/usr/bin/env python3

#*******************************************************
# BSD license:
# Copyright 2026, Dr. Christoph Mayer, LIB Bonn, Germany
# Contributors: see README.md
# See License.txt for the full license text.
#*******************************************************

# Tests of the Python module: module attributes, and site pattern frequencies and effective lengths
# compared with a simple count in Python, for DNA and proteins.
# Run with:  python test_python_module.py   (or with pytest)
# Requires the module built with:  make python   or   pip install .

import itertools
import os
import random
import re
import tempfile

import numpy as np
import QuartetPatternCounterModule as QPC

HERE = os.path.dirname(os.path.abspath(__file__))


def test_module_attributes():
    assert QPC.DNA_PATTERN_ORDER == "ACGT"
    assert QPC.AA_PATTERN_ORDER == "ACDEFGHIKLMNPQRSTVWY"
    assert re.fullmatch(r"\d+\.\d+\.\d+", QPC.__version__)
    # The version is set in QuartetPattern-counter-core.hpp. Only checked in the source folder.
    header = os.path.join(HERE, "QuartetPattern-counter-core.hpp")
    if os.path.exists(header):
        with open(header) as f:
            assert f'#define QPC_VERSION "{QPC.__version__}"' in f.read()


def expected(sequences, alphabet, quartets):
    """Count the site patterns of each quartet in Python."""
    k = len(alphabet)
    index = {c: i for i, c in enumerate(alphabet)}
    frequencies = np.zeros((len(quartets), k ** 4), dtype=np.float32)
    lengths = np.zeros(len(quartets), dtype=np.uint32)
    for row, quartet in enumerate(quartets):
        for site in zip(*(sequences[i].upper() for i in quartet)):
            if all(c in index for c in site):
                column = 0
                for c in site:
                    column = column * k + index[c]
                frequencies[row, column] += 1
                lengths[row] += 1
        if lengths[row] > 0:
            frequencies[row] /= lengths[row]
    return frequencies, lengths


def check_against_python_count(data_type, alphabet, extra_symbols, num_seq=7, length=60, seed=1):
    rng = random.Random(seed)
    symbols = alphabet + alphabet.lower() + extra_symbols
    sequences = ["".join(rng.choice(symbols) for _ in range(length)) for _ in range(num_seq)]
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "alignment.fas")
        with open(path, "w") as f:
            f.write("".join(f">s{i}\n{s}\n" for i, s in enumerate(sequences)))
        counter = QPC.QuartetPatternCounter(path, data_type, 2, 0)
        assert counter.get_dataset_characteristics()[:2] == (num_seq, length)

        batch_size = 8  # Smaller than the number of quartets, so several batches are needed.
        freq = np.empty((batch_size, len(alphabet) ** 4), dtype=np.float32)
        lengths = np.empty(batch_size, dtype=np.uint32)
        quartets = np.empty((batch_size, 4), dtype=np.uint16)
        all_f, all_l, all_q = [], [], []
        while True:
            n = counter.next_quartet_pattern_frequencies_with_lengths_with_quartet_combinations(
                freq, lengths, quartets)
            if n == 0:
                break
            all_f.append(freq[:n].copy())
            all_l.append(lengths[:n].copy())
            all_q.extend(tuple(int(i) for i in q) for q in quartets[:n])

    assert all_q == list(itertools.combinations(range(num_seq), 4))
    exp_f, exp_l = expected(sequences, alphabet, all_q)
    assert np.array_equal(np.concatenate(all_l), exp_l)
    assert np.allclose(np.concatenate(all_f), exp_f, atol=1e-6)


def test_dna_against_python_count():
    check_against_python_count('D', QPC.DNA_PATTERN_ORDER, "-NR")


def test_protein_against_python_count():
    check_against_python_count('P', QPC.AA_PATTERN_ORDER, "-X")


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("passed:", name)
    print("All tests passed.")
