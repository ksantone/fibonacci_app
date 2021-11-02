cmake_minimum_required(VERSION 3.4)
project(pybind_tut)

add_subdirectory(pybind11)
add_subdirectory(DataStructures)

set(ALGORITHMS_executables)
include(algorithms.cmake)

foreach (algo ${ALGORITHMS_executables})
	pybind11_add_module(denovo_sequencing ${algo}.cpp)
	target_include_directories(denovo_sequencing PRIVATE DataStructures/)
	target_link_libraries(denovo_sequencing PRIVATE data_structures)
endforeach(algo)