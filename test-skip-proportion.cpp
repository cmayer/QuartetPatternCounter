// Tests for the skip_proportion option of the QuartetPatternCounter (version 1.8).
// Plain C++, no test framework required. Build and run with:  make run_test_skip

#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include "QuartetPattern-counter-core.hpp"

static int failures = 0;

#define CHECK(cond)                                                              \
  do {                                                                           \
    if (!(cond)) {                                                               \
      ++failures;                                                                \
      std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK failed: " #cond "\n"; \
    }                                                                            \
  } while (0)

static const std::string filename = "random_12x200.fas";
static const size_t total_quartets = 495; // C(12,4)

struct Result
{
  std::vector<std::array<uint16_t, 4>> quartets;
  std::vector<float_type_freq>         freqs;
  std::vector<size_t>                  batches;
};

static Result run_all(QuartetPatternCounter &qpc, size_t N)
{
  Result r;
  std::vector<float_type_freq> freq(N*256);
  std::vector<uint32_t>        lengths(N);
  std::vector<uint16_t>        quartets(N*4);
  size_t n;
  while ((n = qpc.next_quartet_pattern_frequencies_with_lengths(freq.data(), lengths.data(), N, quartets.data())) != 0)
  {
    r.batches.push_back(n);
    for (size_t i=0; i<n; ++i)
      r.quartets.push_back({quartets[4*i], quartets[4*i+1], quartets[4*i+2], quartets[4*i+3]});
    r.freqs.insert(r.freqs.end(), freq.begin(), freq.begin() + n*256);
  }
  return r;
}

static size_t index_of(const std::vector<std::array<uint16_t, 4>> &all, const std::array<uint16_t, 4> &q)
{
  for (size_t i=0; i<all.size(); ++i)
    if (all[i] == q)
      return i;
  return (size_t)-1;
}

int main()
{
  // Reference: no skipping.
  QuartetPatternCounter ref_counter(filename, 'D', 1, 0);
  Result ref = run_all(ref_counter, 50);
  CHECK(ref.quartets.size() == total_quartets);
  CHECK(ref_counter.get_skip_proportion() == 0.0);

  // q = 0.3, N = 50: full batches draw round(50/0.7) = 71 quartets and keep 50.
  // 6 full batches use 426 quartets. The last 69 keep round(69*50/71) = 49.
  {
    QuartetPatternCounter qpc(filename, 'D', 4, 0, 0.3, 11);
    Result r = run_all(qpc, 50);
    CHECK(r.batches == std::vector<size_t>({50, 50, 50, 50, 50, 50, 49}));

    size_t last = 0;
    bool first = true;
    for (size_t i=0; i<r.quartets.size(); ++i)
    {
      size_t idx = index_of(ref.quartets, r.quartets[i]);
      CHECK(idx != (size_t)-1);
      CHECK(first || idx > last); // Generator order, no duplicates.
      CHECK(std::memcmp(&r.freqs[i*256], &ref.freqs[idx*256], 256*sizeof(float_type_freq)) == 0);
      last = idx;
      first = false;
    }
  }

  // Same seed gives the same quartets, also with a different number of threads.
  {
    QuartetPatternCounter a(filename, 'D', 1, 0, 0.5, 3);
    QuartetPatternCounter b(filename, 'D', 4, 0, 0.5, 3);
    CHECK(run_all(a, 20).quartets == run_all(b, 20).quartets);
  }

  // A skip proportion > 0 requires the quartet array.
  {
    QuartetPatternCounter qpc(filename, 'D', 1, 0, 0.2);
    std::vector<float_type_freq> freq(10*256);
    bool thrown = false;
    try { qpc.next_quartet_pattern_frequencies(freq.data(), 10); }
    catch (const std::runtime_error &) { thrown = true; }
    CHECK(thrown);
  }

  // Invalid skip proportions.
  for (double bad : {-0.1, 1.0, 2.0, std::nan("")})
  {
    bool thrown = false;
    try { QuartetPatternCounter qpc(filename, 'D', 1, 0, bad); }
    catch (const std::runtime_error &) { thrown = true; }
    CHECK(thrown);
  }

  // Four cluster lists must not overlap and must not contain an index twice.
  {
    auto throws = [](std::vector<int> p1, std::vector<int> p2, std::vector<int> p3, std::vector<int> p4, const char *expected)
    {
      try
      {
        QuartetPatternCounter qpc(filename, p1.data(), p1.size(), p2.data(), p2.size(),
                                  p3.data(), p3.size(), p4.data(), p4.size(), 'D', 1, 0);
      }
      catch (const std::runtime_error &e)
      {
        return std::string(e.what()).find(expected) != std::string::npos;
      }
      return false;
    };
    CHECK(throws({0, 1}, {2}, {3}, {1},  "occurs in list 1 and in list 4"));
    CHECK(throws({0}, {2, 5}, {5}, {3},  "occurs in list 2 and in list 3"));
    CHECK(throws({0}, {2}, {3, 3}, {4},  "occurs more than once in list 3"));
    CHECK(!throws({0, 1}, {2}, {3}, {4}, ""));
  }

  if (failures == 0)
    std::cout << "All skip_proportion tests passed.\n";
  else
    std::cout << failures << " check(s) failed.\n";
  return failures == 0 ? 0 : 1;
}
