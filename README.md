# Quartet Pattern Counter

Author: Christoph Mayer, Leibniz Insitute or the analysis of biodiversity change
Version 1.8.0

The Quartet Pattern Counter computes site pattern frequencies of quartets of sequences in a multiple sequence alignment. For every quartet combination (four sequences of the alignment) it returns the relative frequencies of all 256 (DNA) or 160000 (protein) site patterns.

It is a C++17 library with Python bindings (pybind11). There is also a stand-alone program for alignments of exactly four sequences.

License: BSD 3-clause, see [License.txt](License.txt).

## What you can to with the Quartet Pattern Counter (QPC).

The QPC has the following different applicaitons.

### Stand alone program for counting all site pattern frequences of given four sequences in a fasta file


### C++ library and Python bindings for counting site pattern frequencies for all combinations of four sequences in a fasta file 

#### Python bindings

#### C++ library


### C++ library and Python bindings  counting site pattern frequencies for  selected combinations of four sequences in a fasta file 



- All C(n,4) quartet combinations of an alignment with n sequences, or only the four cluster
  quartet combinations: one sequence from each of four disjoint lists of sequences.
- Results are computed in batches into arrays you pass in. In Python these are numpy arrays filled
  in place, without copying.
- Multithreaded counting.
- Optional effective lengths: the number of sites used for each quartet combination.
- Optional random skipping of a fixed proportion of quartet combinations (`skip_proportion`).

## Counting rules

- The input is a FASTA file. All sequences must have the same length, and sequence names must be
  unique.
- Sites with a gap, ambiguity code or other invalid symbol in any of the four sequences are skipped.
  Lower case symbols are counted as upper case. RNA (U) is not supported.
- Frequencies are relative to the number of sites that were counted (the effective length). If no
  site can be used, all frequencies are 0 and the effective length is 0.

## Output layout

Each row of the frequency array belongs to one quartet combination (s1, s2, s3, s4) of sequence
indices, where indices refer to the order of the sequences in the FASTA file, starting at 0.
The site pattern (x1, x2, x3, x4), with xi the symbol of sequence si, is stored in column

- DNA: `64*x1 + 16*x2 + 4*x3 + x4`, with A=0, C=1, G=2, T=3, so the columns are AAAA, AAAC, ..., TTTT.
- Proteins: `8000*x1 + 400*x2 + 20*x3 + x4`, with the amino acids in alphabetical order of their
  one-letter codes, ACDEFGHIKLMNPQRSTVWY.

## Building

Requirements: a C++17 compiler. For the Python module also Python 3, pybind11 and numpy.

Install the Python module with pip (pip installs pybind11 and numpy if needed):

```bash
pip install .
```

Or use make:

```bash
make standalone      # Stand-alone program quartet-pattern-counter
make python          # Python module QuartetPatternCounterModule in this folder (calls setup.py)
make run_test_skip   # C++ tests, no test framework needed
make run_catch2_tests  # Catch2 tests
make all             # Python module and the Catch2 test program
make clean
```

`CXX`, `CXXFLAGS` and `LDFLAGS` can be set in the environment or on the make command line.
`-std=c++17` is always added.

Python tests (after installing the module): `pytest test_python_module.py test_skip_proportion.py`.

The version is set in one place, `QPC_VERSION` in `QuartetPattern-counter-core.hpp`. `setup.py` reads it
from there.

`catch.hpp` (Catch2 v2.13.10 single header) is included. To use another copy, pass its directory:
`make all CATCH2_DIR=/path/to/Catch2-single-header`.

## Python

```python
import numpy as np
import QuartetPatternCounterModule as QPC

counter = QPC.QuartetPatternCounter("alignment.fas", data_type='D', num_threads=4, verbosity=0)
num_sequences, alignment_length = counter.get_dataset_characteristics()

batch_size = 10000
frequencies          = np.empty((batch_size, 256), dtype=np.float32)   # 160000 columns for proteins
lengths              = np.empty(batch_size, dtype=np.uint32)
quartet_combinations = np.empty((batch_size, 4), dtype=np.uint16)

while True:
    n = counter.next_quartet_pattern_frequencies_with_lengths_with_quartet_combinations(
            frequencies, lengths, quartet_combinations)
    if n == 0:
        break
    # Rows 0..n-1 are valid: frequencies[i] belongs to the sequences quartet_combinations[i].
```

The arrays must have exactly these dtypes and be C-contiguous, otherwise an error is raised. Each
call fills as many rows as possible and returns the number of filled rows. It returns 0 when all
quartet combinations have been processed.

Module attributes: `QPC.__version__` (e.g. `"1.8.0"`), and the symbol orders of the columns,
`QPC.DNA_PATTERN_ORDER` (`"ACGT"`) and `QPC.AA_PATTERN_ORDER` (`"ACDEFGHIKLMNPQRSTVWY"`).

Functions for the next batch:

| Function | Fills |
|---|---|
| `next_quartet_pattern_frequencies(freq)` | frequencies |
| `next_quartet_pattern_frequencies_with_lengths(freq, lengths)` | frequencies, effective lengths |
| `next_quartet_pattern_frequencies_with_quartet_combinations(freq, quartet_combinations)` | frequencies, quartet combinations |
| `next_quartet_pattern_frequencies_with_lengths_with_quartet_combinations(freq, lengths, quartet_combinations)` | all three |

Four cluster quartet combinations: pass four disjoint integer arrays of sequence indices. Each
quartet combination contains one sequence of each array, and the four arrays must not overlap.

```python
counter = QPC.QuartetPatternCounter("alignment.fas",
                                    np.array([0, 1]), np.array([2]), np.array([3, 4]), np.array([5]),
                                    data_type='D', num_threads=4, verbosity=0)
```

Constructor arguments:

| Argument | Meaning |
|---|---|
| `filename` | FASTA file |
| `partition1` ... `partition4` | Only for four cluster quartet combinations: disjoint arrays of sequence indices |
| `data_type` | `'D'` for DNA, `'P'` for proteins |
| `num_threads` | Number of threads |
| `verbosity` | 0 for no output, larger values for more |
| `skip_proportion` | Proportion of quartet combinations skipped at random, default 0 |
| `seed` | Seed for the random skipping, default 0 |

## Skipping quartet combinations

With `skip_proportion = q > 0`, only a random subset of the quartet combinations is analysed. For
an array with N rows, each call draws the next round(N/(1-q)) quartet combinations and chooses
exactly N of them uniformly at random, keeping their order. So every full batch skips the same
number of quartet combinations. When fewer quartet combinations remain for the last batch, the
same proportion of them is kept.

Because the caller has to know which quartet combinations were kept, skipping requires the functions
that return the quartet combinations (`..._with_quartet_combinations`). The other functions raise an
error if `skip_proportion > 0`.

The selection is reproducible: the same `seed` gives the same quartet combinations on every
platform and for any number of threads.

```python
counter = QPC.QuartetPatternCounter("alignment.fas", 'D', 4, 0, skip_proportion=0.9, seed=42)
```

## C++

Include `QuartetPattern-counter-core.hpp` and compile together with `QuartetPattern-counter-core.cpp`
(with `-pthread`). The class `QuartetPatternCounter` has the same constructors as in Python. The
optional quartet combination and length arrays are passed as pointers:

```cpp
QuartetPatternCounter qpc("alignment.fas", 'D', 4, 0, /*skip_proportion=*/0.5, /*seed=*/1);
std::vector<float_type_freq> freqs(N*qpc.get_num_cols());
std::vector<uint32_t>        lengths(N);
std::vector<uint16_t>        quartet_combinations(N*4);
size_t n;
while ((n = qpc.next_quartet_pattern_frequencies_with_lengths(freqs.data(), lengths.data(), N,
                                                              quartet_combinations.data())) != 0)
{
  // Use rows 0..n-1.
}
```

Errors are reported as exceptions derived from `std::exception` (mostly `std::runtime_error`).

## Stand-alone program

```bash
quartet-pattern-counter [-p] infile.fas outfile.npy
```

The input file must contain exactly four sequences. The output is a 1-dimensional float32 numpy
array with the 256 (DNA) or 160000 (protein, `-p`) site pattern frequencies. The command line and
output are the same as in version 1.1.

## Files

| File | Content |
|---|---|
| `QuartetPattern-counter-core.hpp/.cpp` | Library |
| `quartet-combination-generators.h` | Generators of the quartet combinations |
| `QP_python_bindings.cpp`, `setup.py`, `pyproject.toml`, `MANIFEST.in` | Python module |
| `quartet-pattern-counter.cpp` | Stand-alone program |
| `test-skip-proportion.cpp`, `test_skip_proportion.py`, `test_python_module.py` | Tests |
| `test-QuartetPattern-counter_catch2.cpp`, `catch.hpp` | Catch2 tests |
| `.github/workflows/tests.yml` | Tests on GitHub Actions (Ubuntu and macOS) |
| `QuartetPattern_PythonTester.py` | Python usage examples |
| `simple_*.fas`, `random_12x200.fas` | Example and test alignments |
