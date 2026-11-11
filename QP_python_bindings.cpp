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
// with the distribution.
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
// File: QP_python_bindings.cpp
// Python bindings for QuartetPatternCounter class
// ============================================================================
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "QuartetPattern-counter-core.hpp"

namespace py = pybind11;

// Python wrapper for QuartetPatternCounter class
// The output arrays are filled in place. They must have exactly this dtype and be C-contiguous.
// Together with py::arg().noconvert() this makes pybind11 reject other arrays instead of
// silently filling a temporary converted copy.
typedef py::array_t<float_type_freq, py::array::c_style> freq_array;
typedef py::array_t<uint32_t,        py::array::c_style> length_array;
typedef py::array_t<uint16_t,        py::array::c_style> quartet_combination_array;

// The partition arrays are only read. Other integer types (e.g. numpy's default int64) and
// non-contiguous views (e.g. arr[::2]) are converted by pybind11 into a contiguous int32 copy.
// Without c_style a view with the right dtype would be passed unchanged and read as if contiguous.
typedef py::array_t<int32_t, py::array::c_style | py::array::forcecast> partition_array;

class PyQuartetPatternCounter
{
  private:
    QuartetPatternCounter myQPcounter;

    void check_columns(const py::buffer_info &buf_freq) const
    {
        if ((size_t)buf_freq.shape[1] != myQPcounter.get_num_cols())
        {
            throw std::runtime_error("Frequency array must have shape (batch_size, " + std::to_string(myQPcounter.get_num_cols()) +
                                     ") for this data type, but has " + std::to_string(buf_freq.shape[1]) + " columns.");
        }
    }
    
    static QuartetPatternCounter create_QuartetPatternCounter(
        const std::string &filename,
        partition_array arr1, 
        partition_array arr2,
        partition_array arr3, 
        partition_array arr4,
        char dataType, int numThreads, int verbosity,
        double skip_proportion, uint64_t seed
        )
    {
        py::buffer_info buf1 = arr1.request();
        py::buffer_info buf2 = arr2.request();
        py::buffer_info buf3 = arr3.request();
        py::buffer_info buf4 = arr4.request();       

        // Error checking
        if (buf1.ndim != 1) throw std::runtime_error("partition1 must be a 1D array");
        if (buf2.ndim != 1) throw std::runtime_error("partition2 must be a 1D array");
        if (buf3.ndim != 1) throw std::runtime_error("partition3 must be a 1D array");
        if (buf4.ndim != 1) throw std::runtime_error("partition4 must be a 1D array");

        size_t size1 = buf1.size;
        size_t size2 = buf2.size;
        size_t size3 = buf3.size;
        size_t size4 = buf4.size;
        
        // Extract pointers and sizes
        const int32_t* ptr1 = static_cast<const int32_t*>(buf1.ptr);
        const int32_t* ptr2 = static_cast<const int32_t*>(buf2.ptr);
        const int32_t* ptr3 = static_cast<const int32_t*>(buf3.ptr);
        const int32_t* ptr4 = static_cast<const int32_t*>(buf4.ptr);
        
        // Construct and return C++ object
        return QuartetPatternCounter(filename, ptr1, size1, ptr2, size2, ptr3, size3, ptr4, size4, dataType, numThreads, verbosity,
                                     skip_proportion, seed);
    }
    
    
public:
    // Constructor for computing all quartet frequencies
    PyQuartetPatternCounter(const std::string &filename, char dataType, int numThreads, int verbosity,
                            double skip_proportion, uint64_t seed) : myQPcounter(filename, dataType, numThreads, verbosity, skip_proportion, seed)
    {}

    // Constructor for computing four cluster quartet frequencies
    PyQuartetPatternCounter(const std::string &filename,
                partition_array partition1, partition_array partition2,
                partition_array partition3, partition_array partition4,
                char dataType, int numThreads, int verbosity,
                double skip_proportion, uint64_t seed) : myQPcounter(create_QuartetPatternCounter(filename, partition1, partition2, partition3, partition4,
                                                                          dataType, numThreads, verbosity, skip_proportion, seed))
    {}
   
    
   
    std::pair<int32_t, int32_t> get_dataset_characteristics() const
    {
        return myQPcounter.get_dataset_characteristics();
    }

    // *************************************************
    // We can implement different variants of this function by implementing different parameter lists
    // or by using py::optional<py::array_t<double>> array_opt. 
    // Finally, I decided to implement this using different paramters lists, since this can improve code optimisation by the compiler. 
    // *************************************************

    int32_t next_quartet_pattern_frequencies(freq_array output)
    {
        // Get buffer info from the numpy array
        auto buf_freq = output.request(true);
        
        // Verify it's 2D
        if (buf_freq.ndim != 2)
        {
            throw std::runtime_error("In PyQuartetPatternCounter::next_quartet_pattern_frequencies: Output array must be 2D with shape (batch_size, 256 or 160000).");
        }
        
        size_t rows_freq = buf_freq.shape[0];
        check_columns(buf_freq);
//        size_t cols = buf_freq.shape[1];
        
        // Get pointer to the array data (no copy!)
        float_type_freq* ptr_freq = static_cast<float_type_freq*>(buf_freq.ptr);
        
        // Fill the array in-place
        return myQPcounter.next_quartet_pattern_frequencies(ptr_freq, rows_freq);
    }


    int32_t next_quartet_pattern_frequencies_with_quartet_combinations(freq_array output, quartet_combination_array quartet_combinations)
    {
        // Get buffer info from the numpy array
        auto buf_freq = output.request(true);
        auto buf_quart = quartet_combinations.request(true);
        
        // Verify it's 2D
        if (buf_freq.ndim != 2)
        {
            throw std::runtime_error("In PyQuartetPatternCounter::next_quartet_pattern_frequencies_with_quartet_combinations: Output array must be 2D with shape (batch_size, 256 or 160000).");
        }
        if (buf_quart.ndim != 2)
        {
            throw std::runtime_error("In PyQuartetPatternCounter::next_quartet_pattern_frequencies_with_quartet_combinations: quartet combination array must be 2D with shape (batch_size, 4)");
        }
        if (buf_quart.shape[1] != 4)
        {
            throw std::runtime_error("In PyQuartetPatternCounter::next_quartet_pattern_frequencies_with_quartet_combinations: quartet combination array must have shape (batch_size, 4)");
        }
          
        size_t rows_freq = buf_freq.shape[0];
        check_columns(buf_freq);

        if (rows_freq != (size_t)buf_quart.shape[0])
        {
           throw std::runtime_error("In PyQuartetPatternCounter::next_quartet_pattern_frequencies_with_quartet_combinations: "
                                    "Quartet combination array has to have the same size as the frequency array in the first dimension.");
        }

//        size_t cols = buf_freq.shape[1];

        // Get pointer to the array data (no copy!)
        float_type_freq* ptr_freq = static_cast<float_type_freq*>(buf_freq.ptr);
        uint16_t*        ptr_quart = static_cast<uint16_t*>(buf_quart.ptr);

        // Fill the array in-place
        return myQPcounter.next_quartet_pattern_frequencies(ptr_freq, rows_freq, ptr_quart);
    }


    int32_t next_quartet_pattern_frequencies_with_lengths(freq_array frequencies, length_array lengths)
    {
        // Get buffer info from the numpy array
        auto buf_freq = frequencies.request(true);
        auto buf_len = lengths.request(true);
        
        // Verify it's 2D
        if (buf_freq.ndim != 2)
        {
            throw std::runtime_error("In PyQuartetPatternCounter::next_quartet_pattern_frequencies_with_lengths: Output array must be 2D with shape (batch_size, 256 or 160000).");
        }
        if (buf_len.ndim != 1)
        {
            throw std::runtime_error("In PyQuartetPatternCounter::next_quartet_pattern_frequencies_with_lengths: length array must be 1D with shape (batch_size)");
        }
        
        size_t rows_freq = buf_freq.shape[0];
        check_columns(buf_freq);
//        size_t cols = buf_freq.shape[1];
        
        if (rows_freq != (size_t)buf_len.shape[0])
        {
           throw std::runtime_error("In PyQuartetPatternCounter::next_quartet_pattern_frequencies_with_lengths: "
                                    "Lengths array has to have the same size as the frequency "
                                    "array in the first dimension.");
        }
        
        // Get pointer to the array data (no copy!)
        float_type_freq*  ptr_freq = static_cast<float_type_freq*>(buf_freq.ptr);
        uint32_t*         ptr_len  = static_cast<uint32_t*>(buf_len.ptr);
        
        // Fill the array in-place
        return myQPcounter.next_quartet_pattern_frequencies_with_lengths(ptr_freq, ptr_len, rows_freq);
    }

    int32_t next_quartet_pattern_frequencies_with_lengths_with_quartet_combinations(freq_array frequencies,
                                                                        length_array lengths,
                                                                        quartet_combination_array quartet_combinations)
    {
        // Get buffer info from the numpy array
        auto buf_freq = frequencies.request(true);
        auto buf_len  = lengths.request(true);
        auto buf_quart = quartet_combinations.request(true);
        
        // Verify it's 2D
        if (buf_freq.ndim != 2)
        {
            throw std::runtime_error("In PyQuartetPatternCounter::next_quartet_pattern_frequencies_with_lengths_with_quartet_combinations: Output array must be 2D with shape (batch_size, 256 or 160000).");
        }
        if (buf_quart.ndim != 2)
        {
            throw std::runtime_error("In PyQuartetPatternCounter::next_quartet_pattern_frequencies_with_lengths_with_quartet_combinations: quartet combination array must be 2D with shape (batch_size, 4)");
        }
        if (buf_quart.shape[1] != 4)
        {
            throw std::runtime_error("In PyQuartetPatternCounter::next_quartet_pattern_frequencies_with_lengths_with_quartet_combinations: quartet combination array must have shape (batch_size, 4)");
        }
        if (buf_len.ndim != 1)
        {
            throw std::runtime_error("In PyQuartetPatternCounter::next_quartet_pattern_frequencies_with_lengths_with_quartet_combinations: length array must be 1D with shape (batch_size)");
        }
        
        size_t rows_freq = buf_freq.shape[0];
        check_columns(buf_freq);
//        size_t cols = buf_len.shape[1];
        
        if (rows_freq != (size_t)buf_len.shape[0])
        {
           throw std::runtime_error("In PyQuartetPatternCounter::next_quartet_pattern_frequencies_with_lengths_with_quartet_combinations: "
                                    "Lengths array has to have the same size as the frequency array in the first dimension.");
        }

        if (rows_freq != (size_t)buf_quart.shape[0])
        {
           throw std::runtime_error("In PyQuartetPatternCounter::next_quartet_pattern_frequencies_with_lengths_with_quartet_combinations: "
                                    "Quartet combination array has to have the same size as the frequency array in the first dimension.");
        }
        
        // Get pointer to the array data (no copy!)
        float_type_freq*  ptr_freq = static_cast<float_type_freq*>(buf_freq.ptr);
        uint32_t*         ptr_len  = static_cast<uint32_t*>(buf_len.ptr);
        uint16_t*         ptr_quart = static_cast<uint16_t*>(buf_quart.ptr);
        
        // Fill the array in-place
        return myQPcounter.next_quartet_pattern_frequencies_with_lengths(ptr_freq, ptr_len, rows_freq, ptr_quart);
    }

    
    size_t get_num_cols() const {
        return myQPcounter.get_num_cols();
    }

    double get_skip_proportion() const {
        return myQPcounter.get_skip_proportion();
    }
    
    size_t get_call_count() const {
        return myQPcounter.get_call_count();
    }
    
    float_type_freq get_total_processing_time() const {
        return myQPcounter.get_total_processing_time();
    }
    
};


PYBIND11_MODULE(QuartetPatternCounterModule, m) {
    m.doc() = "High-performance computation module with stateful processing";
    
    // Expose the QuartetPatternCounter class
    py::class_<PyQuartetPatternCounter>(m, "QuartetPatternCounter")
        // Bind constructor for computing all quartet frequencies
        .def(py::init<std::string, char, int, int, double, uint64_t>(), 
             py::arg("filename"),
             py::arg("data_type"),
             py::arg("num_threads"),
             py::arg("verbosity"),
             py::arg("skip_proportion") = 0.0,
             py::arg("seed") = QuartetPatternCounter::default_seed,
             "Create a new QuartetPatternCounter Object for computing all quartet frequencies\n\n"
             "data_type:       'D' for DNA (256 site patterns) or 'P' for proteins (160000 site patterns).\n"
             "num_threads:     Number of threads used to count site patterns.\n"
             "verbosity:       0 for no output, larger values for more output.\n"
             "skip_proportion: Proportion of quartet combinations skipped at random, 0 <= skip_proportion < 1.\n"
             "                 For an array with N rows, round(N/(1-skip_proportion)) quartet combinations are drawn\n"
             "                 and exactly N of them are chosen at random, keeping their order. If > 0, only the\n"
             "                 *_with_quartet_combinations functions can be used.\n"
             "seed:            Seed of the random number generator used for skipping.")

        // Bind constructor for computing four cluster quartet frequencies
        .def(py::init<std::string,
                     partition_array, partition_array,
                     partition_array, partition_array, char, int, int, double, uint64_t>(), 
             py::arg("filename"),
             py::arg("partition1"), py::arg("partition2"), py::arg("partition3"), py::arg("partition4"),
             py::arg("data_type"),
             py::arg("num_threads"),
             py::arg("verbosity"),
             py::arg("skip_proportion") = 0.0,
             py::arg("seed") = QuartetPatternCounter::default_seed,
             "Create a new QuartetPatternCounter Object for four cluster quartet frequencies\n\n"
             "partition1..4:   Disjoint 1D integer arrays of sequence indices. Each quartet combination contains\n"
             "                 one sequence of each partition.\n"
             "data_type:       'D' for DNA (256 site patterns) or 'P' for proteins (160000 site patterns).\n"
             "num_threads:     Number of threads used to count site patterns.\n"
             "verbosity:       0 for no output, larger values for more output.\n"
             "skip_proportion: Proportion of quartet combinations skipped at random, 0 <= skip_proportion < 1.\n"
             "                 For an array with N rows, round(N/(1-skip_proportion)) quartet combinations are drawn\n"
             "                 and exactly N of them are chosen at random, keeping their order. If > 0, only the\n"
             "                 *_with_quartet_combinations functions can be used.\n"
             "seed:            Seed of the random number generator used for skipping.")
  
        .def("next_quartet_pattern_frequencies", &PyQuartetPatternCounter::next_quartet_pattern_frequencies,
             "Fill a pre-allocated numpy array in-place with pattern frequencies (no copy)",
             py::arg("freq").noconvert())            
        .def("next_quartet_pattern_frequencies_with_quartet_combinations", &PyQuartetPatternCounter::next_quartet_pattern_frequencies_with_quartet_combinations,
             "Fill pre-allocated numpy arrays in-place with pattern frequencies and the quartet combinations (sequence indices) they belong to (no copy)",
             py::arg("freq").noconvert(), py::arg("quartet_combinations").noconvert())

        .def("next_quartet_pattern_frequencies_with_lengths", &PyQuartetPatternCounter::next_quartet_pattern_frequencies_with_lengths,
             "Fill a pre-allocated numpy array in-place with pattern frequencies and effective lengths(no copy)",
             py::arg("freq").noconvert(), py::arg("lengths").noconvert())
        .def("next_quartet_pattern_frequencies_with_lengths_with_quartet_combinations", &PyQuartetPatternCounter::next_quartet_pattern_frequencies_with_lengths_with_quartet_combinations,
             "Fill pre-allocated numpy arrays in-place with pattern frequencies, effective lengths and the quartet combinations (sequence indices) they belong to (no copy)",
             py::arg("freq").noconvert(), py::arg("lengths").noconvert(), py::arg("quartet_combinations").noconvert())
        
        .def("get_dataset_characteristics", &PyQuartetPatternCounter::get_dataset_characteristics,
             "Get number of sequences and sequnce lengths")
        .def("get_skip_proportion", &PyQuartetPatternCounter::get_skip_proportion,
             "Get the proportion of quartets that are skipped at random")
    
     
     
     
     
         
/* // OLD VERSION: Version using only functions, no classes.
        .def("load_file", &PyQuartetPatternCounter::load_file,
             "Load data from file",
             py::arg("filename"))
        .def("compute_file_stats", &PyQuartetPatternCounter::compute_file_stats,
             "Compute mean, min, max of loaded file")
        .def("process_data", &PyQuartetPatternCounter::process_data,
             "Process arrays with loaded file data",
             py::arg("arr1"),
             py::arg("arr2"),
             py::arg("arr3"),
             py::arg("arr4"))
        .def("set_num_cols", &PyQuartetPatternCounter::set_num_cols,
             "Set number of output columns",
             py::arg("num_cols"))
        .def("get_num_cols", &PyQuartetPatternCounter::get_num_cols,
             "Get number of output columns")
        .def("get_call_count", &PyQuartetPatternCounter::get_call_count,
             "Get number of process_data calls")
        .def("get_total_processing_time", &PyQuartetPatternCounter::get_total_processing_time,
             "Get total processing time in seconds")
        .def("reset", &PyQuartetPatternCounter::reset,
             "Reset all state")
        .def("process_data", &PyQuartetPatternCounter::process_data,
             "Process arrays with loaded file data",
             py::arg("arr1"),
             py::arg("arr2"),
             py::arg("arr3"),
             py::arg("arr4"))
        .def("assign_data", &PyQuartetPatternCounter::assign_data,
             "Fill a pre-allocated numpy array in-place (no copy)",
             py::arg("output"))
*/
         ;
             

}
