//*******************************************************
// BSD license:
// Copyright 2026, Dr. Christoph Mayer, Nikita Kulikov, Dr. Iker Irisarri LIB, Bonn/Hamburg, Germany
//
// Redistribution and use in source and binary forms, with or without modification, are
// permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this list of
// conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice, this list of 
// conditions and the following disclaimer in the documentation and/or other materials provided 
//with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors may be used to 
// endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED 
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY 
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//*******************************************************

// ============================================================================
// File: test_compute_catch2.cpp
// Unit tests using Catch2 framework (header-only version)
// ============================================================================

#include "QuartetPattern-counter-core.hpp"
#include <fstream>
#include <vector>
#include "symbol-cartesian-product_string.h"
#include <iomanip>

extern char order_of_AA_site_patterns; // 'a': alphabetical, 'p': paml

// Helper to create test file
void create_test_file(const std::string& filename, const std::vector<double>& data)
{
  std::ofstream file(filename);
  for (double val : data)
  {
    file << val << "\n";
  }
  file.close();
}

void print_pattern_frequencies_and_lengths_DEBUG(std::ostream &os, float_type_freq *freq, uint32_t *lengths, char dataType, size_t N_combinations, uint16_t *quartet_combination_array=nullptr)
{
  cartesian_product    cp;
  vector<string>       all_patterns;
  int                  num_cols;

  if (dataType == 'P')
  {
    if (order_of_AA_site_patterns == 'a')
      cp.operator()(4, "ACDEFGHIKLMNPQRSTVWY", all_patterns);
    else
      cp.operator()(4, "ARNDCQEGHILKMFPSTWYV", all_patterns);
    num_cols = 160000;
  }
  else
  {
    cp.operator()(4, "ACGT", all_patterns);
    num_cols = 256;
  }

  if (lengths)
    for (size_t r = 0; r < N_combinations; ++r)
    {
      os << "Quartet combination: " << r+1 << "\t" << " length: " << lengths[r] << endl;
      if (quartet_combination_array)
        os << "Combination indices: " << quartet_combination_array[0+r*4] << "," << quartet_combination_array[1+r*4] << "," << quartet_combination_array[2+r*4] << "," << quartet_combination_array[3+r*4] << "\n";
      for (int i=0; i<num_cols; ++i)
      {
        os << std::setw(8) << i << " " << all_patterns[i] << "\t" << freq[r*num_cols+i] << endl;
      }
    }
  else
    for (size_t r = 0; r < N_combinations; ++r)
    {
      os << "Combination: " << r+1 << endl;
      for (int i=0; i<num_cols; ++i)
      {
        os << all_patterns[i] << "\t" << freq[r*num_cols+i] << endl;
      }
    }
}


//#define CATCH_CONFIG_MAIN  // Tells Catch2 to provide main()
#include "catch.hpp"       // Catch2 single header version. Set CATCH2_DIR in the makefile.
                                                                             // OR if using catch_amalgamated.hpp:

#ifdef CATCH_CONFIG_MAIN
// #include "catch_amalgamated.hpp"

TEST_CASE("QuartetPatternCounter construction", "[QuartetPatternCounter][constructor]")
{
  SECTION("DNA data constructor")
  {
    std::string filename = "simple_4.fas";
    QuartetPatternCounter data_DNA(filename);
    CHECK(data_DNA.get_num_cols() == 256);
    CHECK(data_DNA.get_call_count() == 0);
    CHECK(data_DNA.get_dataset_characteristics()==std::pair<size_t,unsigned>(4,5));
  }
}

TEST_CASE("QuartetPatternCounter construct thread collection", "[QuartetPatternCounter][construct thread collection]")
{
  SECTION("DNA data constructor")
  {
    std::string filename = "simple_4.fas";
    QuartetPatternCounter data_DNA(filename);
    std::vector<float_type_freq> freqs(10*data_DNA.get_num_cols()); // Room for 10 quartets.
    data_DNA.next_quartet_pattern_frequencies(freqs.data(), 10);
  }


/*
*  SECTION("Constructor with custom num_cols")
*  {
*    DataProcessor processor(512);
*    CHECK(processor.get_num_cols() == 512);
*    CHECK_FALSE(processor.has_file_loaded());
*  }
*  */
}

// ============================================================================
// Tests with expected values computed by hand.
// ============================================================================

#include <algorithm>
#include <array>
#include <cstdio>
#include <filesystem>
#include <map>
#include <stdexcept>
#include <type_traits>

using Catch::Matchers::Contains;

// Writes a FASTA file to the temporary directory and returns its path.
static std::string write_temp_fasta(const std::string &name, const std::string &content)
{
  std::string path = (std::filesystem::temp_directory_path() / ("qpc_catch2_test_" + name + ".fas")).string();
  std::ofstream os(path);
  os << content;
  return path;
}

// Column of the DNA site pattern, e.g. "AAGG", in a frequency row. A=0, C=1, G=2, T=3.
static size_t DNA_column(const char *pattern)
{
  const std::string symbols = "ACGT";
  size_t col = 0;
  for (int i=0; i<4; ++i)
    col = 4*col + symbols.find(pattern[i]);
  return col;
}

// Column of the protein site pattern in a frequency row. Alphabetical order of the one-letter codes.
static size_t AA_column(const char *pattern)
{
  const std::string symbols = "ACDEFGHIKLMNPQRSTVWY";
  size_t col = 0;
  for (int i=0; i<4; ++i)
    col = 20*col + symbols.find(pattern[i]);
  return col;
}

// Checks that a frequency row has exactly the given non-zero values.
static void check_row(const float_type_freq *row, size_t num_cols, const std::map<size_t, float> &expected)
{
  for (size_t c=0; c<num_cols; ++c)
  {
    auto it = expected.find(c);
    float exp = (it == expected.end()) ? 0.0f : it->second;
    INFO("column " << c);
    CHECK(row[c] == Approx(exp));
  }
}

struct All_batches
{
  std::vector<std::array<uint16_t, 4>> quartet_combinations;
  std::vector<float_type_freq>         frequencies;
  std::vector<uint32_t>                lengths;
  std::vector<size_t>                  batch_sizes;
};

// Reads all batches with the given batch size.
static All_batches read_all_batches(QuartetPatternCounter &qpc, size_t N)
{
  All_batches r;
  const size_t cols = qpc.get_num_cols();
  std::vector<float_type_freq> freq(N*cols);
  std::vector<uint32_t>        lengths(N);
  std::vector<uint16_t>        quartet_combinations(N*4);
  size_t n;
  while ((n = qpc.next_quartet_pattern_frequencies_with_lengths(freq.data(), lengths.data(), N, quartet_combinations.data())) != 0)
  {
    r.batch_sizes.push_back(n);
    for (size_t i=0; i<n; ++i)
      r.quartet_combinations.push_back({quartet_combinations[4*i], quartet_combinations[4*i+1],
                                        quartet_combinations[4*i+2], quartet_combinations[4*i+3]});
    r.frequencies.insert(r.frequencies.end(), freq.begin(), freq.begin() + n*cols);
    r.lengths.insert(r.lengths.end(), lengths.begin(), lengths.begin() + n);
  }
  return r;
}

TEST_CASE("DNA frequencies of simple_4.fas", "[QuartetPatternCounter][DNA]")
{
  // s1 ACGTA, s2 ACGTA, s3 ACGTG, s4 ACGTG: the sites are AAAA, CCCC, GGGG, TTTT, AAGG.
  QuartetPatternCounter qpc("simple_4.fas", 'D', 1, 0);
  std::vector<float_type_freq> freq(10*256);
  std::vector<uint32_t>        lengths(10);
  std::vector<uint16_t>        quartet_combinations(10*4);

  REQUIRE(qpc.next_quartet_pattern_frequencies_with_lengths(freq.data(), lengths.data(), 10, quartet_combinations.data()) == 1);
  CHECK(lengths[0] == 5);
  CHECK(std::vector<uint16_t>(quartet_combinations.begin(), quartet_combinations.begin()+4) == std::vector<uint16_t>{0, 1, 2, 3});
  check_row(freq.data(), 256, {{DNA_column("AAAA"), 0.2f}, {DNA_column("CCCC"), 0.2f}, {DNA_column("GGGG"), 0.2f},
                               {DNA_column("TTTT"), 0.2f}, {DNA_column("AAGG"), 0.2f}});

  // All quartet combinations have been processed.
  CHECK(qpc.next_quartet_pattern_frequencies_with_lengths(freq.data(), lengths.data(), 10, quartet_combinations.data()) == 0);
  CHECK(qpc.get_call_count() == 2);
}

TEST_CASE("Protein frequencies of simple_4.fas", "[QuartetPatternCounter][protein]")
{
  // A, C, G, T are also amino acids (Ala, Cys, Gly, Thr).
  QuartetPatternCounter qpc("simple_4.fas", 'P', 1, 0);
  REQUIRE(qpc.get_num_cols() == 160000);
  std::vector<float_type_freq> freq(2*160000);
  std::vector<uint32_t>        lengths(2);

  REQUIRE(qpc.next_quartet_pattern_frequencies_with_lengths(freq.data(), lengths.data(), 2) == 1);
  CHECK(lengths[0] == 5);
  check_row(freq.data(), 160000, {{AA_column("AAAA"), 0.2f}, {AA_column("CCCC"), 0.2f}, {AA_column("GGGG"), 0.2f},
                                  {AA_column("TTTT"), 0.2f}, {AA_column("AAGG"), 0.2f}});
}

TEST_CASE("Sites with gaps are skipped", "[QuartetPatternCounter][gaps]")
{
  // simple_5_gap.fas: s1 TAAAAA, s2 AGAAAA, s3 AACAA-, s4 -AATAA, s5 TAAAGA
  QuartetPatternCounter qpc("simple_5_gap.fas", 'D', 2, 0);
  All_batches r = read_all_batches(qpc, 10);

  REQUIRE(r.quartet_combinations.size() == 5);

  // Quartet combination (0,1,2,3): sites 1 to 4 are used, site 0 and 5 have a gap.
  CHECK(r.quartet_combinations[0] == std::array<uint16_t, 4>{0, 1, 2, 3});
  CHECK(r.lengths[0] == 4);
  check_row(&r.frequencies[0], 256, {{DNA_column("AGAA"), 0.25f}, {DNA_column("AACA"), 0.25f},
                                     {DNA_column("AAAT"), 0.25f}, {DNA_column("AAAA"), 0.25f}});

  // Quartet combination (0,1,2,4): sites 0 to 4 are used, site 5 has a gap in s3.
  CHECK(r.quartet_combinations[1] == std::array<uint16_t, 4>{0, 1, 2, 4});
  CHECK(r.lengths[1] == 5);
  check_row(&r.frequencies[256], 256, {{DNA_column("TAAT"), 0.2f}, {DNA_column("AGAA"), 0.2f}, {DNA_column("AACA"), 0.2f},
                                       {DNA_column("AAAA"), 0.2f}, {DNA_column("AAAG"), 0.2f}});
}

TEST_CASE("Lower case symbols are counted as upper case", "[QuartetPatternCounter][input]")
{
  std::string path = write_temp_fasta("lower", ">s1\nacgta\n>s2\nACGTA\n>s3\nacgtg\n>s4\nACGTG\n");
  QuartetPatternCounter qpc(path, 'D', 1, 0);
  std::vector<float_type_freq> freq(256);
  REQUIRE(qpc.next_quartet_pattern_frequencies(freq.data(), 1) == 1);
  check_row(freq.data(), 256, {{DNA_column("AAAA"), 0.2f}, {DNA_column("CCCC"), 0.2f}, {DNA_column("GGGG"), 0.2f},
                               {DNA_column("TTTT"), 0.2f}, {DNA_column("AAGG"), 0.2f}});
  std::remove(path.c_str());
}

TEST_CASE("No valid site gives zero frequencies and length 0", "[QuartetPatternCounter][gaps]")
{
  std::string path = write_temp_fasta("nosite", ">s1\nA-N\n>s2\n-AA\n>s3\nAA-\n>s4\nAAA\n");
  QuartetPatternCounter qpc(path, 'D', 1, 0);
  std::vector<float_type_freq> freq(256, 1.0f);
  uint32_t length = 99;
  REQUIRE(qpc.next_quartet_pattern_frequencies_with_lengths(freq.data(), &length, 1) == 1);
  CHECK(length == 0);
  check_row(freq.data(), 256, {});
  std::remove(path.c_str());
}

TEST_CASE("Order of all quartet combinations", "[QuartetPatternCounter][generator]")
{
  QuartetPatternCounter qpc("simple_5.fas", 'D', 1, 0);
  All_batches r = read_all_batches(qpc, 10);
  CHECK(r.quartet_combinations == std::vector<std::array<uint16_t, 4>>{{0,1,2,3}, {0,1,2,4}, {0,1,3,4}, {0,2,3,4}, {1,2,3,4}});
}

TEST_CASE("Batch size and number of threads do not change the results", "[QuartetPatternCounter][batches]")
{
  QuartetPatternCounter reference("random_12x200.fas", 'D', 1, 0);
  All_batches ref = read_all_batches(reference, 1000);
  REQUIRE(ref.quartet_combinations.size() == 495); // C(12,4)
  CHECK(ref.batch_sizes == std::vector<size_t>{495});

  QuartetPatternCounter qpc("random_12x200.fas", 'D', 4, 0);
  All_batches r = read_all_batches(qpc, 7);
  CHECK(r.batch_sizes.size() == 71); // 70 batches of 7 and one of 5
  CHECK(r.batch_sizes.back() == 5);
  CHECK(r.quartet_combinations == ref.quartet_combinations);
  CHECK(r.lengths == ref.lengths);
  CHECK(r.frequencies == ref.frequencies);
}

TEST_CASE("Four cluster quartet combinations", "[QuartetPatternCounter][four cluster]")
{
  QuartetPatternCounter all("random_12x200.fas", 'D', 1, 0);
  All_batches ref = read_all_batches(all, 1000);

  int p1[] = {0, 1}, p2[] = {2}, p3[] = {3, 4}, p4[] = {5};
  QuartetPatternCounter qpc("random_12x200.fas", p1, 2, p2, 1, p3, 2, p4, 1, 'D', 1, 0);
  All_batches r = read_all_batches(qpc, 10);

  // One sequence of each list, the last list changes fastest.
  REQUIRE(r.quartet_combinations == std::vector<std::array<uint16_t, 4>>{{0,2,3,5}, {0,2,4,5}, {1,2,3,5}, {1,2,4,5}});

  // The frequencies are the same as for these quartet combinations in the run with all quartet combinations.
  for (size_t i=0; i<r.quartet_combinations.size(); ++i)
  {
    size_t j = std::find(ref.quartet_combinations.begin(), ref.quartet_combinations.end(), r.quartet_combinations[i])
               - ref.quartet_combinations.begin();
    REQUIRE(j < ref.quartet_combinations.size());
    CHECK(r.lengths[i] == ref.lengths[j]);
    CHECK(std::equal(&r.frequencies[i*256], &r.frequencies[(i+1)*256], &ref.frequencies[j*256]));
  }
}

TEST_CASE("Skipping quartet combinations", "[QuartetPatternCounter][skip]")
{
  QuartetPatternCounter all("random_12x200.fas", 'D', 1, 0);
  All_batches ref = read_all_batches(all, 1000);

  // q = 0.5, N = 20: each full batch draws 40 quartet combinations and keeps 20.
  // 12 full batches use 480. The last 15 keep round(15*20/40) = 8.
  QuartetPatternCounter qpc("random_12x200.fas", 'D', 4, 0, 0.5, 123);
  CHECK(qpc.get_skip_proportion() == 0.5);
  All_batches r = read_all_batches(qpc, 20);
  std::vector<size_t> expected(12, 20);
  expected.push_back(8);
  CHECK(r.batch_sizes == expected);

  // Kept quartet combinations are in generator order and have the right frequencies.
  size_t j = 0;
  for (size_t i=0; i<r.quartet_combinations.size(); ++i)
  {
    while (j < ref.quartet_combinations.size() && ref.quartet_combinations[j] != r.quartet_combinations[i])
      ++j;
    REQUIRE(j < ref.quartet_combinations.size());
    CHECK(std::equal(&r.frequencies[i*256], &r.frequencies[(i+1)*256], &ref.frequencies[j*256]));
    ++j;
  }

  // The same seed gives the same selection.
  QuartetPatternCounter again("random_12x200.fas", 'D', 1, 0, 0.5, 123);
  CHECK(read_all_batches(again, 20).quartet_combinations == r.quartet_combinations);
}

TEST_CASE("QuartetPatternCounter cannot be copied", "[QuartetPatternCounter][copy]")
{
  // A copy would point into the sequences of the original.
  STATIC_REQUIRE_FALSE(std::is_copy_constructible<QuartetPatternCounter>::value);
  STATIC_REQUIRE_FALSE(std::is_copy_assignable<QuartetPatternCounter>::value);
}

TEST_CASE("The filename may be a temporary string", "[QuartetPatternCounter][input]")
{
  QuartetPatternCounter qpc(std::string("simple_") + "4.fas", 'D', 1, 0);
  std::vector<float_type_freq> freq(256);
  CHECK(qpc.next_quartet_pattern_frequencies(freq.data(), 1) == 1);
}

TEST_CASE("Invalid input is rejected", "[QuartetPatternCounter][errors]")
{
  SECTION("Missing file")
  {
    CHECK_THROWS_WITH(QuartetPatternCounter("does_not_exist.fas", 'D', 1, 0), Contains("Could not open the input file"));
  }
  SECTION("Unknown data type")
  {
    CHECK_THROWS_WITH(QuartetPatternCounter("simple_4.fas", 'X', 1, 0), Contains("Unknown data type"));
  }
  SECTION("Sequences of unequal length")
  {
    std::string path = write_temp_fasta("unequal", ">s1\nACGT\n>s2\nACG\n>s3\nACGT\n>s4\nACGT\n");
    CHECK_THROWS_WITH(QuartetPatternCounter(path, 'D', 1, 0), Contains("do not have equal lengths"));
    std::remove(path.c_str());
  }
  SECTION("Duplicate sequence names")
  {
    std::string path = write_temp_fasta("duplicate", ">s1\nACGT\n>s2\nACGT\n>s1\nACGT\n>s4\nACGT\n");
    CHECK_THROWS_WITH(QuartetPatternCounter(path, 'D', 1, 0), Contains("occurs more than once in the input file"));
    std::remove(path.c_str());
  }
  SECTION("Fewer than four sequences")
  {
    std::string path = write_temp_fasta("three", ">s1\nACGT\n>s2\nACGT\n>s3\nACGT\n");
    CHECK_THROWS_AS(QuartetPatternCounter(path, 'D', 1, 0), std::invalid_argument);
    std::remove(path.c_str());
  }
  SECTION("Four cluster index out of range")
  {
    int p1[] = {0}, p2[] = {1}, p3[] = {2}, p4[] = {5};
    CHECK_THROWS_WITH(QuartetPatternCounter("simple_5.fas", p1, 1, p2, 1, p3, 1, p4, 1, 'D', 1, 0), Contains("out of range"));
  }
  SECTION("Four cluster lists overlap")
  {
    int p1[] = {0, 1}, p2[] = {1}, p3[] = {2}, p4[] = {3};
    CHECK_THROWS_WITH(QuartetPatternCounter("simple_5.fas", p1, 2, p2, 1, p3, 1, p4, 1, 'D', 1, 0), Contains("must not overlap"));
  }
  SECTION("Invalid skip proportion")
  {
    CHECK_THROWS_WITH(QuartetPatternCounter("simple_5.fas", 'D', 1, 0, 1.0), Contains("skip proportion"));
    CHECK_THROWS_WITH(QuartetPatternCounter("simple_5.fas", 'D', 1, 0, -0.5), Contains("skip proportion"));
  }
  SECTION("Skipping without quartet combination array")
  {
    QuartetPatternCounter qpc("simple_5.fas", 'D', 1, 0, 0.5);
    std::vector<float_type_freq> freq(10*256);
    CHECK_THROWS_WITH(qpc.next_quartet_pattern_frequencies(freq.data(), 10), Contains("skip proportion"));
  }
}

#endif


/*
 *TEST_CASE("DataProcessor construction", "[DataProcessor][constructor]")
 *{
 *  SECTION("Default constructor")
 *  {
 *    DataProcessor processor;
 *    CHECK(processor.get_num_cols() == 256);
 *    CHECK_FALSE(processor.has_file_loaded());
 *    CHECK(processor.get_call_count() == 0);
 *  }
 *
 *  SECTION("Constructor with custom num_cols")
 *  {
 *    DataProcessor processor(512);
 *    CHECK(processor.get_num_cols() == 512);
 *    CHECK_FALSE(processor.has_file_loaded());
 *  }
 *}
 *
 *TEST_CASE("DataProcessor file loading", "[DataProcessor][load]")
 *{
 *  DataProcessor processor;
 *
 *  SECTION("Load valid file")
 *  {
 *    std::vector<double> test_data = {1.0, 2.0, 3.0, 4.0, 5.0};
 *    create_test_file("test_load.txt", test_data);
 *
 *    processor.load_file("test_load.txt");
 *
 *    CHECK(processor.has_file_loaded());
 *    auto file_data = processor.get_file_data();
 *    CHECK(file_data.size() == 5);
 *    CHECK(file_data[0] == Catch::Detail::Approx(1.0));
 *    CHECK(file_data[4] == Catch::Detail::Approx(5.0));
 *  }
 *
 *  SECTION("Load empty file")
 *  {
 *    create_test_file("test_empty_load.txt", {});
 *
 *    processor.load_file("test_empty_load.txt");
 *
 *    CHECK(processor.has_file_loaded()); // File is loaded, just empty
 *    CHECK(processor.get_file_data().empty());
 *  }
 *
 *  SECTION("Load nonexistent file throws")
 *  {
 *    REQUIRE_THROWS_AS(
 *                      processor.load_file("nonexistent.txt"),
 *                      std::runtime_error
 *                      );
 *    CHECK_FALSE(processor.has_file_loaded());
 *  }
 *
 *  SECTION("Load multiple files in sequence")
 *  {
 *    std::vector<double> data1 = {1.0, 2.0};
 *    std::vector<double> data2 = {10.0, 20.0, 30.0};
 *    create_test_file("test_file1.txt", data1);
 *    create_test_file("test_file2.txt", data2);
 *
 *    processor.load_file("test_file1.txt");
 *    CHECK(processor.get_file_data().size() == 2);
 *
 *    processor.load_file("test_file2.txt");
 *    CHECK(processor.get_file_data().size() == 3);
 *    CHECK(processor.get_file_data()[0] == Catch::Detail::Approx(10.0));
 *  }
 *}
 *
 *TEST_CASE("DataProcessor statistics", "[DataProcessor][stats]")
 *{
 *  DataProcessor processor;
 *
 *  SECTION("Stats on loaded file")
 *  {
 *    std::vector<double> test_data = {1.0, 2.0, 3.0, 4.0, 5.0};
 *    create_test_file("test_stats_proc.txt", test_data);
 *
 *    processor.load_file("test_stats_proc.txt");
 *    auto [mean, min_val, max_val] = processor.compute_file_stats();
 *
 *    CHECK(mean == Catch::Detail::Approx(3.0));
 *    CHECK(min_val == Catch::Detail::Approx(1.0));
 *    CHECK(max_val == Catch::Detail::Approx(5.0));
 *  }
 *
 *  SECTION("Stats before loading file")
 *  {
 *    auto [mean, min_val, max_val] = processor.compute_file_stats();
 *    CHECK(mean == Catch::Detail::Approx(0.0));
 *  }
 *}
 *
 *TEST_CASE("DataProcessor process_data", "[DataProcessor][process]")
 *{
 *  DataProcessor processor(256);
 *
 *  std::vector<double> file_data = {10.0, 20.0, 30.0};
 *  create_test_file("test_process_class.txt", file_data);
 *  processor.load_file("test_process_class.txt");
 *
 *  double arr1[] = {1.0, 2.0, 3.0};
 *  double arr2[] = {4.0, 5.0, 6.0};
 *  double arr3[] = {7.0, 8.0, 9.0};
 *  double arr4[] = {0.1, 0.2, 0.3};
 *
 *  size_t N = 3;
 *  size_t num_cols = 256;
 *
 *  std::vector<double> output1(N * num_cols);
 *  std::vector<double> output2(N * num_cols);
 *  std::vector<double> output3(N * num_cols);
 *
 *  SECTION("Basic processing")
 *  {
 *    processor.process_data(
 *                           arr1, 3, arr2, 3, arr3, 3, arr4, 3,
 *                           output1.data(), output2.data(), output3.data(), N
 *                           );
 *
 *    // Check first element
 *    CHECK(output1[0] == Catch::Detail::Approx(1.0 + 4.0 + 10.0 + 0));
 *
 *    // Check call count incremented
 *    CHECK(processor.get_call_count() == 1);
 *
 *    // Check processing time tracked
 *    //    std::cout << "Value: " << processor.get_total_processing_time() << std::endl;
 *    CHECK(processor.get_total_processing_time() >= 0.0);
 *  }
 *
 *  SECTION("Multiple calls increment state")
 *  {
 *    processor.process_data(arr1, 3, arr2, 3, arr3, 3, arr4, 3,
 *                           output1.data(), output2.data(), output3.data(), N);
 *
 *    CHECK(processor.get_call_count() == 1);
 *    double time1 = processor.get_total_processing_time();
 *
 *    processor.process_data(arr1, 3, arr2, 3, arr3, 3, arr4, 3,
 *                           output1.data(), output2.data(), output3.data(), N);
 *
 *    CHECK(processor.get_call_count() == 2);
 *    CHECK(processor.get_total_processing_time() >= time1);
 *  }
 *
 *  SECTION("Processing without file loaded")
 *  {
 *    DataProcessor empty_processor;
 *
 *    empty_processor.process_data(arr1, 3, arr2, 3, arr3, 3, arr4, 3,
 *                                 output1.data(), output2.data(), output3.data(), N);
 *
 *    // Should work, file_val will be 0
 *    CHECK(output1[0] == Catch::Detail::Approx(1.0 + 4.0 + 0.0 + 0));
 *    CHECK(empty_processor.get_call_count() == 1);
 *  }
 *}
 *
 *
 *TEST_CASE("DataProcessor configuration", "[DataProcessor][config]")
 *{
 *  DataProcessor processor;
 *
 *  SECTION("Change num_cols")
 *  {
 *    CHECK(processor.get_num_cols() == 256);
 *
 *    processor.set_num_cols(512);
 *    CHECK(processor.get_num_cols() == 512);
 *
 *    processor.set_num_cols(128);
 *    CHECK(processor.get_num_cols() == 128);
 *  }
 *
 *  SECTION("num_cols affects output shape")
 *  {
 *    std::vector<double> file_data = {5.0};
 *    create_test_file("test_cols.txt", file_data);
 *    processor.load_file("test_cols.txt");
 *
 *    double arr[] = {1.0};
 *    std::vector<double> output1(128);
 *    std::vector<double> output2(128);
 *    std::vector<double> output3(128);
 *
 *    processor.set_num_cols(128);
 *    processor.process_data(arr, 1, arr, 1, arr, 1, arr, 1,
 *                           output1.data(), output2.data(), output3.data(), 1);
 *
 *    // Check that all 128 columns are filled
 *    CHECK(output1[0] == Catch::Detail::Approx(1.0 + 1.0 + 5.0 + 0));    // j=0
 *    CHECK(output1[127] == Catch::Detail::Approx(1.0 + 1.0 + 5.0 + 127)); // j=127
 *  }
 *}
 *
 *
 *TEST_CASE("DataProcessor reset", "[DataProcessor][reset]")
 *{
 *  DataProcessor processor;
 *
 *  std::vector<double> file_data = {10.0, 20.0};
 *  create_test_file("test_reset.txt", file_data);
 *  processor.load_file("test_reset.txt");
 *
 *  double arr[] = {1.0};
 *  std::vector<double> output1(256);
 *  std::vector<double> output2(256);
 *  std::vector<double> output3(256);
 *
 *  processor.process_data(arr, 1, arr, 1, arr, 1, arr, 1,
 *                         output1.data(), output2.data(), output3.data(), 1);
 *  processor.process_data(arr, 1, arr, 1, arr, 1, arr, 1,
 *                         output1.data(), output2.data(), output3.data(), 1);
 *
 *  SECTION("Reset clears all state")
 *  {
 *    CHECK(processor.has_file_loaded());
 *    CHECK(processor.get_call_count() == 2);
 *    CHECK(processor.get_total_processing_time() >= 0.0);
 *
 *    processor.reset();
 *
 *    CHECK_FALSE(processor.has_file_loaded());
 *    CHECK(processor.get_call_count() == 0);
 *    CHECK(processor.get_total_processing_time() == Catch::Detail::Approx(0.0));
 *    CHECK(processor.get_file_data().empty());
 *  }
 *
 *  SECTION("Can use processor after reset")
 *  {
 *    processor.reset();
 *
 *    processor.load_file("test_reset.txt");
 *    CHECK(processor.has_file_loaded());
 *
 *    processor.process_data(arr, 1, arr, 1, arr, 1, arr, 1,
 *                           output1.data(), output2.data(), output3.data(), 1);
 *    CHECK(processor.get_call_count() == 1);
 *  }
 *}
 *
 *
 *TEST_CASE("DataProcessor caching", "[DataProcessor][cache]")
 *{
 *  DataProcessor processor(10); // Small for testing
 *
 *  std::vector<double> file_data = {5.0};
 *  create_test_file("test_cache.txt", file_data);
 *  processor.load_file("test_cache.txt");
 *
 *  double arr[] = {1.0, 2.0};
 *  size_t N = 2;
 *  std::vector<double> output1(N * 10);
 *  std::vector<double> output2(N * 10);
 *  std::vector<double> output3(N * 10);
 *
 *  processor.process_data(arr, 2, arr, 2, arr, 2, arr, 2,
 *                         output1.data(), output2.data(), output3.data(), N);
 *
 *  SECTION("Last results are cached")
 *  {
 *    const auto& cached1 = processor.get_last_result1();
 *    const auto& cached2 = processor.get_last_result2();
 *    const auto& cached3 = processor.get_last_result3();
 *
 *    CHECK(cached1.size() == N * 10);
 *    CHECK(cached2.size() == N * 10);
 *    CHECK(cached3.size() == N * 10);
 *
 *    // Verify cached values match output
 *    CHECK(cached1[0] == Catch::Detail::Approx(output1[0]));
 *    CHECK(cached1[19] == Catch::Detail::Approx(output1[19]));
 *  }
 *
 *  SECTION("Cache updates on new call")
 *  {
 *    double arr_new[] = {10.0, 20.0};
 *    processor.process_data(arr_new, 2, arr_new, 2, arr_new, 2, arr_new, 2,
 *                           output1.data(), output2.data(), output3.data(), N);
 *
 *    const auto& cached = processor.get_last_result1();
 *    // New values should be different
 *    CHECK(cached[0] == Catch::Detail::Approx(10.0 + 10.0 + 5.0 + 0));
 *  }
 *
 *  SECTION("Reset clears cache")
 *  {
 *    CHECK_FALSE(processor.get_last_result1().empty());
 *
 *    processor.reset();
 **
 *    CHECK(processor.get_last_result1().empty());
 *    CHECK(processor.get_last_result2().empty());
 *    CHECK(processor.get_last_result3().empty());
 *  }
 *}
 *
 *
 *TEST_CASE("DataProcessor state isolation", "[DataProcessor][isolation]")
 *{
 *  SECTION("Multiple independent processors")
 *  {
 *    DataProcessor proc1(256);
 *    DataProcessor proc2(128);
 *
 *    std::vector<double> data1 = {10.0};
 *    std::vector<double> data2 = {20.0};
 *    create_test_file("test_proc1.txt", data1);
 *    create_test_file("test_proc2.txt", data2);
 *
 *    proc1.load_file("test_proc1.txt");
 *    proc2.load_file("test_proc2.txt");
 *
 *    CHECK(proc1.get_file_data()[0] == Catch::Detail::Approx(10.0));
 *    CHECK(proc2.get_file_data()[0] == Catch::Detail::Approx(20.0));
 *    CHECK(proc1.get_num_cols() == 256);
 *    CHECK(proc2.get_num_cols() == 128);
 *
 *    double arr[] = {1.0};
 *    std::vector<double> output(256);
 *
 *    proc1.process_data(arr, 1, arr, 1, arr, 1, arr, 1,
 *                       output.data(), output.data(), output.data(), 1);
 *
 *    CHECK(proc1.get_call_count() == 1);
 *    CHECK(proc2.get_call_count() == 0); // Unaffected
 *  }
 *}
 *
 *
 *TEST_CASE("DataProcessor edge cases", "[DataProcessor][edge]")
 *{
 *  DataProcessor processor;
 *
 *  SECTION("Process with N=0") {
 *    std::vector<double> output1, output2, output3;
 *    double arr[] = {1.0};
 *
 *    // Should not crash
 *    processor.process_data(arr, 1, arr, 1, arr, 1, arr, 1,
 *                           output1.data(), output2.data(), output3.data(), 0);
 *
 *    CHECK(processor.get_call_count() == 1);
 *  }
 *
 *  SECTION("Large N value")
 *  {
 *    std::vector<double> file_data = {5.0};
 *    create_test_file("test_large.txt", file_data);
 *    processor.load_file("test_large.txt");
 *
 *    size_t N = 10000;
 *    std::vector<double> arr(N, 1.0);
 *    std::vector<double> output1(N * 256);
 *    std::vector<double> output2(N * 256);
 *    std::vector<double> output3(N * 256);
 *
 *    processor.process_data(arr.data(), N, arr.data(), N, arr.data(), N, arr.data(), N,
 *                           output1.data(), output2.data(), output3.data(), N);
 *
 *    CHECK(processor.get_call_count() == 1);
 *    CHECK(output1[0] == Catch::Detail::Approx(1.0 + 1.0 + 5.0 + 0));
 *  }
 *
 *}
 */


#ifndef CATCH_CONFIG_MAIN

using namespace std;

void test_without_lengths_DNA()
{
  size_t count_res;

  cerr << "Welcome to the test main function.\n";
  cerr << "Testing DNA input sequence file.\n";

  //  std::string filename = "simple_5_gap.fas";
  std::string filename = "simple_4.fas";
  QuartetPatternCounter qpc(filename, 'D', 1, 10);
  float_type_freq *freqs = new float_type_freq[10*256];
  count_res = qpc.next_quartet_pattern_frequencies(freqs, 10);

  cout << "Number of quartet combinations for which frequencies have been stored: " << count_res << endl;
  print_pattern_frequencies_and_lengths_DEBUG(cout, freqs, nullptr, 'D', count_res);
}


void test_with_lengths_DNA()
{
  size_t count_res;

  //std::string filename = "simple_5_gap.fas";
  std::string filename = "simple_4.fas";
  cerr << "Testing DNA input sequence file:\n";
  cerr << filename << endl;

  QuartetPatternCounter qpc(filename, 'D', 1, 10);
  float_type_freq *freqs = new float_type_freq[10*256];
  uint32_t *lengths = new uint32_t[10];
  count_res = qpc.next_quartet_pattern_frequencies_with_lengths(freqs, lengths, 10);

  cout << "Number of quartet combinations for which frequencies have been stored: " << count_res << endl;
  print_pattern_frequencies_and_lengths_DEBUG(cout, freqs, lengths, 'D', count_res);
}

void test_without_lengths_Protein()
{
  size_t count_res;

  //  std::string filename = "simple_5_gap.fas";
  std::string filename = "simple_4.fas";
  cerr << "Testing protein input sequence file:\n";
  cerr << filename << endl;

  QuartetPatternCounter qpc(filename, 'P', 1, 10);
  float_type_freq *freqs = new float_type_freq[10*160000];
  count_res = qpc.next_quartet_pattern_frequencies(freqs, 10);

  cout << "Number of quartet combinations for which frequencies have been stored: " << count_res << endl;
  print_pattern_frequencies_and_lengths_DEBUG(cout, freqs, nullptr, 'P', count_res);
}

void test_with_lengths_Protein()
{
  size_t count_res;

  //  std::string filename = "simple_5_gap.fas";
  std::string filename = "simple_4.fas";
  cerr << "Testing protein input sequence file:\n";
  cerr << filename << endl;

  QuartetPatternCounter qpc(filename, 'P', 1, 10);
  float_type_freq *freqs = new float_type_freq[10*160000];
  uint32_t *lengths = new uint32_t[10];
  count_res = qpc.next_quartet_pattern_frequencies_with_lengths(freqs, lengths, 10);

  cout << "Number of quartet combinations for which frequencies have been stored: " << count_res << endl;
  print_pattern_frequencies_and_lengths_DEBUG(cout, freqs, lengths, 'P', count_res);
}

void test_four_cluster_quartet_freq_without_lengths_DNA()
{
  size_t count_res;

  // Five sequences are needed, since the lists below use the indices 0 to 4.
  std::string filename = "simple_5_gap.fas";
  cerr << "Testing four cluster quartet selection. Testing DNA input sequence file:\n";
  cerr << filename << endl;

/*
*  int part1[] = {0,4};
*  int part2[] = {1};
*  int part3[] = {2};
*  int part4[] = {3};
*  int size1 = 2, size2 = 1, size3 = 1, size4 = 1;
*/

  int part1[] = {0,3};
  int part2[] = {1};
  int part3[] = {2};
  int part4[] = {4};
  int size1 = 2, size2 = 1, size3 = 1, size4 = 1;

  QuartetPatternCounter qpc(filename, part1, size1, part2, size2, part3, size3, part4, size4, 'D', 1, 10);
  float_type_freq *freqs = new float_type_freq[10*256];
  uint32_t *lengths = new uint32_t[10];
  uint16_t *quartet_combination_array = new uint16_t [10*4];
  count_res = qpc.next_quartet_pattern_frequencies_with_lengths(freqs, lengths, 10, quartet_combination_array);

  cout << "Number of quartet combinations for which frequencies have been stored: " << count_res << endl;
  print_pattern_frequencies_and_lengths_DEBUG(cout, freqs, lengths, 'D', count_res, quartet_combination_array);
}

int main()
{
  std::cout << "Welcome to the main function in the test-QuartetPatternCounter-catch2.cpp file:\n";
  //  test_without_lengths_DNA();
  //  test_with_lengths_DNA();

  //  test_without_lengths_Protein();

  test_with_lengths_Protein();

  //  test_four_cluster_quartet_freq_without_lengths_DNA();
  //  test_four_cluster_quartet_freq_without_lengths_DNA();
}

#endif
