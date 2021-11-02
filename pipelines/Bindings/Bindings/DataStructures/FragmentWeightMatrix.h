#include "SpectrumGraph.h"
#include <array>

#ifndef FragmentWeightMatrix_HEADER
#define FragmentWeightMatrix_HEADER

class FragmentWeightMatrix {
public:
	int spectrumPeaks;
	int num_of_threads = 8;
	int smallest_available_thread = 8;
	std::vector<bool> row_complete;
	std::vector<int> available_indices;
	std::vector<double> fragmentWeightMatrix;
	void construct_triangle(int);
	void compute_row(int, int);
	void traverse_fragment_matrix(int);
	std::map<std::pair<int, int>, std::string> edges;
	FragmentWeightMatrix(SpectrumGraph, int);
};

#endif