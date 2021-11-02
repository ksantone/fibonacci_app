#include "FragmentWeightMatrix.h"
#include <thread>
#include <condition_variable>
#include <deque>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <mutex>

std::string desktop_loc = "../../../../../Desktop/fragment_2.txt";
std::string logging_loc = "./logging/logging.txt";
std::ofstream fragment_file;
std::ofstream logging_file;

std::mutex mtx;

double findMaximumSubfragment(std::vector<double> fragmentWeightMatrix, int i, int j, std::map<std::pair<int, int>, std::string> edges, int spectrumPeaks) {
	double currentMaxfragment = -1;
	std::vector<std::pair<int, int>> edges_without_amino_acids;
	transform(edges.begin(), edges.end(), std::back_inserter(edges_without_amino_acids), [&](std::pair<std::pair<int, int>, std::string> edge) { return edge.first; });
	std::pair<int, int> final_pair;
	std::string amino_acid;
	for (int k = j - 1; k >= 0; k--) {
		double current_entry = fragmentWeightMatrix.at(i * spectrumPeaks + k);
		if (current_entry > currentMaxfragment && std::binary_search(edges_without_amino_acids.begin(), edges_without_amino_acids.end(), std::pair<int, int>(k, j))) {
			currentMaxfragment = current_entry + 1;
			final_pair = std::pair<int, int>(i, k);
			amino_acid = edges[std::pair<int, int>(k, j)];
		}
	}
	for (int k = i - 1; k >= 0; k--) {
		double current_entry = fragmentWeightMatrix.at(k * spectrumPeaks + j);
		if (current_entry > currentMaxfragment && std::binary_search(edges_without_amino_acids.begin(), edges_without_amino_acids.end(), std::pair<int, int>(k, i))) {
			currentMaxfragment = current_entry + 1;
			final_pair = std::pair<int, int>(k, j);
			amino_acid = edges[std::pair<int, int>(k, i)];
		}
	}
	if (currentMaxfragment > -1) {
		mtx.lock();
		fragment_file << "Computed fragment at position: (" << i << ", " << j << ") to have " << currentMaxfragment << " peaks." << std::endl;
		mtx.unlock();
	}
	return currentMaxfragment;
}

std::vector<std::string> glue_edge_to_all_sequences(std::string amino_acid, std::vector<std::string> subpeptides, bool left) {
	std::vector<std::string> new_subpeptides;
	for (int i = 0; i < subpeptides.size(); i++) {
		if (left) {
			new_subpeptides.push_back(amino_acid + subpeptides.at(i));
		}
		else {
			new_subpeptides.push_back(subpeptides.at(i) + amino_acid);
		}
	}
	return new_subpeptides;
}

std::vector<std::string> find_sequence_from_fragment_position(std::vector<double> fragmentWeightMatrix, std::map<std::pair<int, int>, std::string> edges, std::vector<std::pair<int, int>> edges_without_amino_acids, double max_num_of_peaks, std::pair<int, int> optimal_fragment_position, int spectrumPeaks, std::string prefix_suffix_connector) {
	std::vector<std::string> empty_string_vector;
	if (max_num_of_peaks == 0) {
		empty_string_vector.push_back(prefix_suffix_connector);
		return empty_string_vector;
	}
	int i = optimal_fragment_position.first, j = optimal_fragment_position.second;
	double max_prev_frag = -std::numeric_limits<double>::infinity();
	std::pair<int, int> prev_frag_pair;
	bool found = false;
	std::vector<std::string> vector_1, vector_2;
	for (int k = i - 1; k >= 0; k--) {
		if (fragmentWeightMatrix.at(k * spectrumPeaks + j) > max_prev_frag && std::binary_search(edges_without_amino_acids.begin(), edges_without_amino_acids.end(), std::pair<int, int>(k, i))) {
			max_prev_frag = fragmentWeightMatrix.at(k * spectrumPeaks + j);
			prev_frag_pair = std::pair<int, int>(k, i);
			found = true;
			vector_1 = glue_edge_to_all_sequences(edges[prev_frag_pair], find_sequence_from_fragment_position(fragmentWeightMatrix, edges, edges_without_amino_acids, max_num_of_peaks - 1, std::pair<int, int>(k, j), spectrumPeaks, prefix_suffix_connector), true);
		}
	}
	for (int k = j - 1; k >= 0; k--) {
		if (fragmentWeightMatrix.at(i * spectrumPeaks + k) > max_prev_frag && std::binary_search(edges_without_amino_acids.begin(), edges_without_amino_acids.end(), std::pair<int, int>(k, j))) {
			max_prev_frag = fragmentWeightMatrix.at(i * spectrumPeaks + k);
			prev_frag_pair = std::pair<int, int>(k, j);
			found = true;
			vector_2 = glue_edge_to_all_sequences(edges[prev_frag_pair], find_sequence_from_fragment_position(fragmentWeightMatrix, edges, edges_without_amino_acids, max_num_of_peaks - 1, std::pair<int, int>(i, k), spectrumPeaks, prefix_suffix_connector), false);
		}
	}
	if (!found) {
		return empty_string_vector;
	}
	vector_1.insert(vector_1.end(), vector_2.begin(), vector_2.end());
	return vector_1;
}

std::condition_variable current_cv;
std::condition_variable threads_completing_cv;
std::vector<bool> triangle_complete;
std::mutex m1, m2, m3;
int threads_executing = 0;
int rows_completed = 0;

void FragmentWeightMatrix::construct_triangle(int section) {
	std::cout << "Computing triangle." << std::endl;
	std::cout << "Num of threads: " << std::to_string(num_of_threads) << std::endl;
	// Eventually change from hardcoded value of 8 to something which works for all edge cases
	int current_row = section * 8;
	for (int i = current_row; i < current_row + num_of_threads; i++) {
		for (int j = i + 1; j < current_row + num_of_threads; j++) {
			if (i == 369 || j == 370) {
				std::cout << "But how?" << std::endl;
				std::cout << "(" + std::to_string(i) + ", " + std::to_string(j) + ")." << std::endl;
			}
			if (num_of_threads < 8) {
				std::cout << "Computing triangle index: (" + std::to_string(i) + ", " + std::to_string(j) + ")." << std::endl;
			}
			// Computing subrow on top of triangle
			fragmentWeightMatrix.at(i * spectrumPeaks + j) = findMaximumSubfragment(fragmentWeightMatrix, i, j, edges, spectrumPeaks);
			// Compute triangle underneath diagonal by symmetry
			fragmentWeightMatrix.at(j * spectrumPeaks + i) = fragmentWeightMatrix.at(i * spectrumPeaks + j);
		}
	}
	// After triangle has been constructed
	if (num_of_threads == 8) {
		triangle_complete.at(section) = true;
	}
	std::cout << "Triangle computation completed." << std::endl;
}

void FragmentWeightMatrix::compute_row(int section, int row) {
	std::cout << "Computing row "+std::to_string(row) << std::endl;
	int actual_row = section * num_of_threads + row; //Compute actual row from section*num_of_threads or something
	for (int j = actual_row + 1; j < spectrumPeaks; j++) {
		// Compute on row to right of diagonal
		double temp = findMaximumSubfragment(fragmentWeightMatrix, actual_row, j, edges, spectrumPeaks);
		fragmentWeightMatrix.at(actual_row * spectrumPeaks + j) = temp;
		// Follows that column below diagonal is the same since matrix is symmetric
		fragmentWeightMatrix.at(j * spectrumPeaks + actual_row) = fragmentWeightMatrix.at(actual_row * spectrumPeaks + j);
	}
	threads_executing -= 1;
	logging_file << "Row " << actual_row << " complete" << std::endl;
	std::cout << "Row "+std::to_string(row)+" computation complete." << std::endl;
	std::cout << "Remaining threads "+std::to_string(threads_executing)<< "." << std::endl;
}

void FragmentWeightMatrix::traverse_fragment_matrix(int id) {
	std::unique_lock<std::mutex> traversal_lock(m3);
	if (num_of_threads == 8) {
		int current_section = (int)(id / num_of_threads);
		int row_index = id % num_of_threads;
		if (row_index == 0) {
			construct_triangle(current_section);
			current_cv.notify_all();
		}
		current_cv.wait(traversal_lock, [&]() { return triangle_complete.at(current_section); });
		compute_row(current_section, row_index);
		threads_completing_cv.notify_all();
	} else {
		int current_section = (int)(id / 8);
		construct_triangle(current_section);
	}
	/*if (row_index==0) {
		if (current_section > 0) {
			current_cv.wait(lk, [&]() {
					for (int i = 0; i < num_of_threads; i++) {
						if (row_complete.at(id - num_of_threads + i) == false) {
							return false;
						}
					}
					return true;
				});
		}
		construct_triangle(current_section);
		current_cv.notify_all();
		smallest_available_thread += 1;
		traverse_fragment_matrix(smallest_available_thread - 1);
	}
	else {*/
	//}
}

FragmentWeightMatrix::FragmentWeightMatrix(SpectrumGraph spectrumGraph, int spectrumPeaks) {
	FragmentWeightMatrix::spectrumPeaks = spectrumPeaks;
	edges = spectrumGraph.spectralEdges;
	fragment_file.open(desktop_loc);
	logging_file.open(logging_loc);
	logging_file << "Total rows: " << spectrumPeaks << std::endl;
	int matrix_sections = (int)(spectrumPeaks / num_of_threads);
	for (int i = 0; i < spectrumPeaks; i++) {
		row_complete.push_back(false);
	}
	for (int i = 0; i < matrix_sections; i++) {
		available_indices.push_back(0);
		triangle_complete.push_back(false);
	}
	for (int i = 0; i < spectrumPeaks * spectrumPeaks; i++) {
		fragmentWeightMatrix.push_back(0);
	}
	for (int i = 1; i < spectrumPeaks; i++) {
		fragmentWeightMatrix.at(i + i * spectrumPeaks) = -std::numeric_limits<double>::infinity();
	}
	std::unique_lock<std::mutex> execution_lock(m1);
	std::unique_lock<std::mutex> completion_lock(m2);
	std::condition_variable threads_executing_cv;
	std::cout << "Just before thread execution!" << std::endl;
	for (int i = 0; i <= (int)(spectrumPeaks/8); i++) {
		std::vector<std::thread> threads;
		int threads_to_execute = 8;
		if (i == (int)(spectrumPeaks / 8)) {
			num_of_threads = (spectrumPeaks - 1) % 8;
			threads.push_back(std::thread(&FragmentWeightMatrix::traverse_fragment_matrix, this, i * 8));
			threads.at(0).join();
		}
		else {
			for (int j = 0; j < num_of_threads; j++) {
				threads_executing_cv.wait(execution_lock, [&]() { return threads_executing < 8; });
				std::cout << "Thread " + std::to_string(j + i * 8) + " has been created!" << std::endl;
				threads.push_back(std::thread(&FragmentWeightMatrix::traverse_fragment_matrix, this, j + i * 8));
				threads_executing += 1;
			}
			for (int j = 0; j < 8; j++) {
				threads.at(j).join();
			}
			threads_completing_cv.wait(completion_lock, [&]() { return threads_executing == 0; });
		}
		std::cout << "Threads completed successfully!" << std::endl;
		try {
			threads.resize(0);
		} catch(const std::exception& e) {
			std::cout << e.what();
		}
		/*for (int j = 0; j < 8; j++) {
			try {
				delete threads[j];
			} catch (const std::exception& e) {
				std::cout << e.what();
			}
		}*/
		/*for (int j = 0; j < 8; j++) {
			threads.erase(threads.begin());
		}*/
		std::cout << "Threads erased successfully!" << std::endl;
	}
	/*for (int i = 0; i < spectrumPeaks - 1; i++) {
		for (int j = i + 1; j < spectrumPeaks; j++) {
			// Compute on row to right of diagonal
			fragmentWeightMatrix.at(i * spectrumPeaks + j) = findMaximumSubfragment(fragmentWeightMatrix, i, j, edges, spectrumPeaks);
			// Follows that column below diagonal is the same since matrix is symmetric
			fragmentWeightMatrix.at(j * spectrumPeaks + i) = fragmentWeightMatrix.at(i * spectrumPeaks + j);
		}
		logging_file << i << std::endl;
	}*/
	double max_num_of_peaks = 0;
	std::pair<int, int> optimal_fragment_position;
	std::vector<std::pair<int, int>> edges_without_amino_acids;
	transform(edges.begin(), edges.end(), std::back_inserter(edges_without_amino_acids), [&](std::pair<std::pair<int, int>, std::string> edge) { return edge.first; });
	std::vector<std::string> prefix_suffix_connectors;
	std::vector<std::pair<int, int>> optimal_fragment_positions;
	std::map<std::pair<int, int>, double> fragment_to_num_of_peaks;
	for (int i = 0; i < spectrumPeaks; i++) {
		for (int j = 0; j < spectrumPeaks; j++) {
			if (fragmentWeightMatrix.at(i * spectrumPeaks + j) == -std::numeric_limits<double>::infinity()) {
				std::cout << fragmentWeightMatrix.at(i * spectrumPeaks + j) << " ";
			}
			else {
				std::cout << (int)fragmentWeightMatrix.at(i * spectrumPeaks + j) << " ";
			}
			if (fragmentWeightMatrix[i * spectrumPeaks + j] > max_num_of_peaks) {
				max_num_of_peaks = fragmentWeightMatrix[i * spectrumPeaks + j];
			} 
			optimal_fragment_position = std::pair<int, int>(i, j);
			optimal_fragment_positions.push_back(optimal_fragment_position);
			fragment_to_num_of_peaks[optimal_fragment_position] = fragmentWeightMatrix[i * spectrumPeaks + j];
		}
		std::cout << std::endl;
	}
	std::vector<std::pair<double, std::pair<int, int>>> actual_optimal_fragment_positions;
	for (int i = 0; i < optimal_fragment_positions.size(); i++) {
		if (fragment_to_num_of_peaks[optimal_fragment_positions.at(i)] > max_num_of_peaks - 3 && std::binary_search(edges_without_amino_acids.begin(), edges_without_amino_acids.end(), optimal_fragment_positions.at(i))) {
			actual_optimal_fragment_positions.push_back(std::pair<double, std::pair<int, int>>(fragment_to_num_of_peaks[optimal_fragment_positions.at(i)], optimal_fragment_positions.at(i)));
			prefix_suffix_connectors.push_back(edges[optimal_fragment_positions.at(i)]);
		}
	}
	//optimal_fragment_positions.erase(std::remove_if(optimal_fragment_positions.begin(), optimal_fragment_positions.end(), [&](std::pair<int, int> fragment_position) { return fragment_to_num_of_peaks[fragment_position] != max_num_of_peaks; }));
	std::cout << "The maximum number of peaks in a potential peptide from the spectrum is: " + std::to_string(max_num_of_peaks);
	std::vector<std::pair<int, int>> optimal_fragment_positions_without_peaks;
	transform(actual_optimal_fragment_positions.begin(), actual_optimal_fragment_positions.end(), std::back_inserter(optimal_fragment_positions_without_peaks), [&](std::pair<double, std::pair<int, int>> actual_optimal_fragment_position) { return actual_optimal_fragment_position.second; });
	std::vector<std::vector<std::string>> best_fragments;
	for (int i = 0; i < optimal_fragment_positions_without_peaks.size(); i++) {
		best_fragments.push_back(find_sequence_from_fragment_position(fragmentWeightMatrix, edges, edges_without_amino_acids, actual_optimal_fragment_positions.at(i).first, optimal_fragment_positions_without_peaks.at(i), spectrumPeaks, prefix_suffix_connectors.at(i)));
	}
	for (int i = 0; i < optimal_fragment_positions_without_peaks.size(); i++) {
		std::cout << "The optimal peptide is with prefix-suffix connector " << prefix_suffix_connectors.at(i) << " is " << best_fragments.at(i).at(0) << "." << std::endl;
	}
	std::string spectraFileName = "C:\\Users\\Kassim Santone\\source\\repos\\BiologicalAlgorithmLibrary\\out\\install\\out\\build\\";
	const size_t last_slash_idx = spectraFileName.rfind('\\');
	std::string directoryName;
	std::string peptideFile = spectraFileName + "temporary_files\\peptide_list.txt";
	std::ofstream file(peptideFile);
	std::cout << peptideFile;
	file << best_fragments.at(0).at(0);
	file.close();
}