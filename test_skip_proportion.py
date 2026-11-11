# Tests for the skip_proportion option of the QuartetPatternCounter (version 1.8).
# Run with:  python test_skip_proportion.py   (or with pytest)
# Requires the module built with:  make python

import itertools
import math
import os

import numpy as np
import QuartetPatternCounterModule as QPC

HERE = os.path.dirname(os.path.abspath(__file__))
FILENAME = os.path.join(HERE, "random_12x200.fas")
NUM_SEQ = 12
ALL_QUARTETS = list(itertools.combinations(range(NUM_SEQ), 4))


def make(skip=0.0, seed=0, threads=1, partitions=None):
    if partitions is None:
        return QPC.QuartetPatternCounter(FILENAME, 'D', threads, 0, skip_proportion=skip, seed=seed)
    return QPC.QuartetPatternCounter(FILENAME, *partitions, 'D', threads, 0, skip_proportion=skip, seed=seed)


def run_all(counter, batch_size, with_lengths=False):
    """Collect all batches. Returns quartets, frequencies, lengths and the batch sizes."""
    freq = np.empty((batch_size, 256), dtype=np.float32)
    lengths = np.empty(batch_size, dtype=np.uint32)
    quartets = np.empty((batch_size, 4), dtype=np.uint16)
    all_q, all_f, all_l, batches = [], [], [], []
    while True:
        if with_lengths:
            n = counter.next_quartet_pattern_frequencies_with_lengths_with_quartet_combinations(
                freq, lengths, quartet_combinations=quartets)
        else:
            n = counter.next_quartet_pattern_frequencies_with_quartet_combinations(freq, quartet_combinations=quartets)
        if n == 0:
            break
        batches.append(n)
        all_q.extend(tuple(int(x) for x in row) for row in quartets[:n])
        all_f.append(freq[:n].copy())
        all_l.append(lengths[:n].copy())
    f = np.concatenate(all_f) if all_f else np.empty((0, 256), np.float32)
    l = np.concatenate(all_l) if all_l else np.empty(0, np.uint32)
    return all_q, f, l, batches


def expected_batches(total, n, q):
    """Batch sizes that the counter must return for a skip proportion q > 0."""
    def round(x):  # Like C++ std::round: halves away from zero. Python's round() rounds halves to even.
        return math.floor(x + 0.5)
    full_draw = round(n / (1.0 - q))
    out = []
    while total > 0:
        if full_draw <= total:
            keep, total = n, total - full_draw
        else:
            keep, total = round(total * n / full_draw), 0
        if keep == 0:
            break
        out.append(keep)
    return out


def reference():
    q, f, l, _ = run_all(make(), 1000, with_lengths=True)
    return {quartet: (f[i], l[i]) for i, quartet in enumerate(q)}


def test_default_is_no_skipping():
    c = make()
    assert c.get_skip_proportion() == 0.0
    q, f, _, batches = run_all(c, 50)
    assert q == ALL_QUARTETS
    assert batches == [50] * 9 + [45]
    # The old functions without quartets still work and give the same frequencies.
    c2 = make()
    freq = np.empty((1000, 256), dtype=np.float32)
    n = c2.next_quartet_pattern_frequencies(freq)
    assert n == len(ALL_QUARTETS)
    assert np.array_equal(freq[:n], f)


def test_skip_batch_sizes_and_subset():
    ref = reference()
    for q_skip, batch_size in [(0.3, 50), (0.5, 7), (0.9, 10), (0.01, 50), (0.5, 1000)]:
        q, f, l, batches = run_all(make(q_skip, seed=1), batch_size, with_lengths=True)
        exp = expected_batches(len(ALL_QUARTETS), batch_size, q_skip)
        assert batches == exp, (q_skip, batch_size, batches, exp)
        # Kept quartets are a subset of all quartets, in generator order, without duplicates.
        idx = [ALL_QUARTETS.index(x) for x in q]
        assert idx == sorted(set(idx))
        # Each frequency row and length belongs to the returned quartet.
        for i, quartet in enumerate(q):
            assert np.array_equal(f[i], ref[quartet][0])
            assert l[i] == ref[quartet][1]
        # Overall proportion kept is close to 1-q.
        assert abs(len(q) / len(ALL_QUARTETS) - (1 - q_skip)) < 0.03, (q_skip, len(q))


def test_seed_reproducible_and_threads_independent():
    a = run_all(make(0.4, seed=7), 30)[0]
    b = run_all(make(0.4, seed=7, threads=4), 30)[0]
    c = run_all(make(0.4, seed=8), 30)[0]
    assert a == b
    assert a != c


def test_selection_is_uniform():
    # Every quartet should be kept with probability about 1-q.
    q_skip, runs = 0.5, 400
    counts = np.zeros(len(ALL_QUARTETS))
    pos = {x: i for i, x in enumerate(ALL_QUARTETS)}
    for seed in range(runs):
        for x in run_all(make(q_skip, seed=seed), 25)[0]:
            counts[pos[x]] += 1
    p = counts / runs
    sd = math.sqrt(0.25 / runs)
    assert abs(p.mean() - 0.5) < 0.01
    assert np.all(np.abs(p - 0.5) < 5 * sd), (p.min(), p.max())


def test_four_cluster_generator():
    parts = (np.array([0, 1, 2]), np.array([3, 4]), np.array([5, 6, 7]), np.array([8, 9, 10, 11]))
    full = run_all(make(partitions=parts), 1000)[0]
    assert len(full) == 3 * 2 * 3 * 4
    q, _, _, batches = run_all(make(0.25, seed=3, partitions=parts), 10)
    assert batches == expected_batches(len(full), 10, 0.25)
    idx = [full.index(x) for x in q]
    assert idx == sorted(set(idx))


def test_skip_requires_quartet_array():
    freq = np.empty((10, 256), dtype=np.float32)
    lengths = np.empty(10, dtype=np.uint32)
    for call in (lambda c: c.next_quartet_pattern_frequencies(freq),
                 lambda c: c.next_quartet_pattern_frequencies_with_lengths(freq, lengths)):
        try:
            call(make(0.2))
        except RuntimeError as e:
            assert "skip proportion" in str(e)
        else:
            raise AssertionError("Expected RuntimeError")


def test_four_cluster_lists_must_not_overlap():
    cases = [(([0, 1], [2], [3], [1]), "occurs in list 1 and in list 4"),
             (([0], [2, 5], [5], [3]), "occurs in list 2 and in list 3"),
             (([0], [2], [3, 3], [4]), "occurs more than once in list 3")]
    for lists, expected in cases:
        try:
            make(partitions=tuple(np.array(x) for x in lists))
        except RuntimeError as e:
            assert expected in str(e), str(e)
        else:
            raise AssertionError(f"Expected RuntimeError for {lists}")


def test_invalid_skip_proportion():
    for bad in (-0.1, 1.0, 1.5, float('nan')):
        try:
            make(bad)
        except RuntimeError as e:
            assert "skip proportion" in str(e)
        else:
            raise AssertionError(f"Expected RuntimeError for {bad}")


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("passed:", name)
    print("All tests passed.")
