# Makefile for building both C++ tests and Python module

# CXX, CXXFLAGS and LDFLAGS can be set in the environment or on the command line.
CXX ?= g++
CXXFLAGS ?= -O3 -Wall -Wextra
# C++17 is always required, also if CXXFLAGS comes from the environment.
override CXXFLAGS := -std=c++17 $(CXXFLAGS)
#CXXFLAGS = -std=c++17 -g3 -ggdb -Wall -Wextra


# For C++ unit tests (simple macros - no dependencies)
#test: test_compute
#	./test_compute

#test_compute: test_compute.cpp QuartetPattern-counter-core.cpp QuartetPattern-counter-core.hpp
#	$(CXX) $(CXXFLAGS) -o test_compute test_compute.cpp QuartetPattern-counter-core.cpp

# Default target: the stand-alone program. Needs only a C++17 compiler (no pybind11, no Catch2).
# Same command line and output as quartet-pattern-counter v1.1: [-p] infile.fas outfile.npy
standalone: quartet-pattern-counter

quartet-pattern-counter: quartet-pattern-counter.cpp QuartetPattern-counter-core.cpp QuartetPattern-counter-core.hpp basic-DNA-RNA-AA-routines.h quartet-combination-generators.h
	$(CXX) $(CXXFLAGS) -pthread -o quartet-pattern-counter quartet-pattern-counter.cpp QuartetPattern-counter-core.cpp $(LDFLAGS)

# Run both tests and build Python module
all: python run_test-QP-counter_catch2

# For Catch2 tests (requires Catch2 installed)
run_test-QP-counter_catch2: test_QP-counter_catch2
	@echo "*****************************"
	@echo "** Run test_compute_catch2 **"
	@echo "*****************************"
#	./test_QP-counter_catch2 -s
	./test_QP-counter_catch2

# Directory containing the Catch2 single header catch.hpp, e.g.: make all CATCH2_DIR=/path/to/Catch2-single-header
CATCH2_DIR ?= .

test_QP-counter_catch2: test-QuartetPattern-counter_catch2.cpp QuartetPattern-counter-core.cpp QuartetPattern-counter-core.hpp
	$(CXX) $(CXXFLAGS) -pthread -I$(CATCH2_DIR) -o test_QP-counter_catch2 test-QuartetPattern-counter_catch2.cpp QuartetPattern-counter-core.cpp $(LDFLAGS)

# The Catch2 test cases in test-QuartetPattern-counter_catch2.cpp. They are compiled only if
# CATCH_CONFIG_MAIN is defined. Then Catch2 provides main() instead of the hand-written one.
run_catch2_tests: catch2_tests
	./catch2_tests

catch2_tests: test-QuartetPattern-counter_catch2.cpp QuartetPattern-counter-core.cpp QuartetPattern-counter-core.hpp quartet-combination-generators.h
	$(CXX) $(CXXFLAGS) -DCATCH_CONFIG_MAIN -pthread -I$(CATCH2_DIR) -o catch2_tests test-QuartetPattern-counter_catch2.cpp QuartetPattern-counter-core.cpp $(LDFLAGS)

# Tests of the skip_proportion option. Plain C++, no Catch2 needed.
run_test_skip: test_skip_proportion
	./test_skip_proportion

test_skip_proportion: test-skip-proportion.cpp QuartetPattern-counter-core.cpp QuartetPattern-counter-core.hpp quartet-combination-generators.h
	$(CXX) $(CXXFLAGS) -pthread -o test_skip_proportion test-skip-proportion.cpp QuartetPattern-counter-core.cpp $(LDFLAGS)

# For Python module (requires pybind11)
python: 
	python3 setup.py build_ext --inplace

# Clean build artifacts
clean:
	rm -f quartet-pattern-counter test_QP-counter_catch2 test_skip_proportion catch2_tests *.o *.so *.pyc
	rm -rf build __pycache__ *.egg-info
	rm -f test_*.txt



.PHONY: standalone test test-QP-counter_catch2 run_test_skip run_catch2_tests python clean all


