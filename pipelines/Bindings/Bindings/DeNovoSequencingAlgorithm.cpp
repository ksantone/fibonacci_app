#include <pybind11/pybind11.h>
#include <iostream>
#include <fstream>
#include <thread>
#include <chrono>
#include <list>
#include "DataStructures/FragmentWeightMatrix.h"
#include "DataStructures/PeptideSequence.h"

namespace py = pybind11;

std::list<PeptideSequence> performDeNovoSequencingOnGraphs(std::list<SpectrumGraph> spectrumGraphs)
{
	return std::list<PeptideSequence>();
}

void run_denovo()
{
	int argc = 2;
	std::string argv = "C:\\Users\\Kassim Santone\\Desktop\\DjangoProjects\\MyOmicsFlows\\src\\spectrum_list_29.txt";//"C:\\Users\\Kassim Santone\\source\\repos\\BiologicalAlgorithmLibrary\\out\\install\\out\\build\\spectrum_list.txt";
	if (argc != 2) {
		std::cout << argc << " Incorrect number of arguments!" << std::endl;
	}
	else if (argc == 2) {
		std::string spectraFileName = argv;
		std::ifstream spectraFile(spectraFileName);
		int first = 0;
		std::list<SpectrumGraph> spectrumGraphs, filteredSpectrumGraphs;
		std::string spectralData;
		while (getline(spectraFile, spectralData) && first <= 1) {
			if (first == 1) {
				spectrumGraphs.push_back(SpectrumGraph(spectralData));
			}
			first += 1;
		}
		std::list<SpectrumGraph>::iterator it;
		std::list<FragmentWeightMatrix> fragmentWeightMatrices;
		for (it = spectrumGraphs.begin(); it != spectrumGraphs.end(); it++) {
			fragmentWeightMatrices.push_back(FragmentWeightMatrix(*it, it->spectralPeaks.size()));
		}
	}
}

PYBIND11_MODULE(denovo_sequencing, handle) {
	handle.doc() = "This is the DeNovo sequencing algorithm.";
	handle.def("denovo_sequencing_algorithm", &run_denovo);

	py::class_<FragmentWeightMatrix>(
		handle, "FragmentWeightMatrix"
		);
}