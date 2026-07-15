#include "include/t0/dssd.h"
#include "include/event/t0/dssd_match.h"

#include <TH1D.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

namespace brill {

namespace {

struct Hit {
	int strip = 0;
	double energy = 0.0;
	double time = 0.0;
};

double NormalizeEnergy(double raw_energy, double p0, double p1, double p2) {
	return p0 + p1 * raw_energy + p2 * raw_energy * raw_energy;
}

int ReadOneSide(
	const std::string &path,
	int strips,
	double *p0,
	double *p1,
	double *p2
) {
	std::ifstream fin(path);
	if (!fin.good()) {
		std::cerr << "Error: Open normalize parameter file " << path << " failed.\n";
		return -1;
	}

	std::string line;
	if (!std::getline(fin, line)) {
		std::cerr << "Error: Read header from " << path << " failed.\n";
		return -1;
	}

	while (std::getline(fin, line)) {
		if (line.empty()) continue;
		std::istringstream iss(line);
		int index = -1;
		double value0 = 0.0;
		double value1 = 0.0;
		double value2 = 0.0;
		if (!(iss >> index >> value0 >> value1 >> value2)) continue;
		if (index < 0 || index >= strips) {
			std::cerr << "Error: Strip index " << index
				<< " out of range in " << path << ".\n";
			return -1;
		}
		p0[index] = value0;
		p1[index] = value1;
		p2[index] = value2;
	}
	return 0;
}

int WriteOneSide(
	const std::string &path,
	int strips,
	const double *p0,
	const double *p1,
	const double *p2
) {
	std::ofstream fout(path);
	if (!fout.good()) {
		std::cerr << "Error: Open output normalize parameter file " << path << " failed.\n";
		return -1;
	}
	fout << "strip p0 p1 p2\n";
	for (int i = 0; i < strips; ++i) {
		fout << i << " " << p0[i] << " " << p1[i] << " " << p2[i] << "\n";
	}
	return 0;
}

void SortHits(int &num, int *strip, double *energy, double *time) {
	std::vector<Hit> hits;
	hits.reserve(num);
	for (int i = 0; i < num; ++i) {
		hits.push_back(Hit{strip[i], energy[i], time[i]});
	}
	std::sort(
		hits.begin(),
		hits.end(),
		[](const Hit &left, const Hit &right) {
			return left.energy > right.energy;
		}
	);
	for (int i = 0; i < num; ++i) {
		strip[i] = hits[i].strip;
		energy[i] = hits[i].energy;
		time[i] = hits[i].time;
	}
}

} // namespace

int WriteDssdNormalizeParameters(
	const std::string &front_path,
	const std::string &back_path,
	const DssdNormalizeParameters &parameters
) {
	if (WriteOneSide(
		front_path,
		parameters.front_strips,
		parameters.front_p0,
		parameters.front_p1,
		parameters.front_p2
	)) {
		return -1;
	}
	if (WriteOneSide(
		back_path,
		parameters.back_strips,
		parameters.back_p0,
		parameters.back_p1,
		parameters.back_p2
	)) {
		return -1;
	}
	return 0;
}

int ReadDssdNormalizeParameters(
	const std::string &front_path,
	const std::string &back_path,
	DssdNormalizeParameters &parameters
) {
	if (ReadOneSide(
		front_path,
		parameters.front_strips,
		parameters.front_p0,
		parameters.front_p1,
		parameters.front_p2
	)) {
		return -1;
	}
	if (ReadOneSide(
		back_path,
		parameters.back_strips,
		parameters.back_p0,
		parameters.back_p1,
		parameters.back_p2
	)) {
		return -1;
	}
	return 0;
}

void ApplyDssdNormalize(
	const DssdEvent &input,
	const DssdNormalizeParameters &parameters,
	DssdEvent &output
) {
	output.front_num = input.front_num;
	output.back_num = input.back_num;
	for (int i = 0; i < input.front_num; ++i) {
		int strip = input.front_strip[i];
		output.front_strip[i] = strip;
		output.front_energy[i] = NormalizeEnergy(
			input.front_energy[i],
			parameters.front_p0[strip],
			parameters.front_p1[strip],
			parameters.front_p2[strip]
		);
		output.front_time[i] = input.front_time[i];
	}
	for (int i = 0; i < input.back_num; ++i) {
		int strip = input.back_strip[i];
		output.back_strip[i] = strip;
		output.back_energy[i] = NormalizeEnergy(
			input.back_energy[i],
			parameters.back_p0[strip],
			parameters.back_p1[strip],
			parameters.back_p2[strip]
		);
		output.back_time[i] = input.back_time[i];
	}
	SortHits(output.front_num, output.front_strip, output.front_energy, output.front_time);
	SortHits(output.back_num, output.back_strip, output.back_energy, output.back_time);
}

namespace {

bool MatchSimple(
	const DssdEvent &input,
	const SquareDetectorConfig &detector,
	DssdMatchEvent &output,
	TH1D *h_energy_diff
) {
	int n_front = input.front_num;
	int n_back = input.back_num;
	const int *front_strip = input.front_strip;
	const int *back_strip = input.back_strip;
	const double *front_e = input.front_energy;
	const double *back_e = input.back_energy;
	const double *front_t = input.front_time;
	const double *back_t = input.back_time;
	double match_tolerance = detector.match_tolerance;

	// sort front hits by strip number (ascending)
	std::vector<int> sorted_front_strip(front_strip, front_strip + n_front);
	std::vector<double> sorted_front_e(front_e, front_e + n_front);
	std::vector<double> sorted_front_t(front_t, front_t + n_front);
	for (int i = 0; i < n_front - 1; ++i) {
		for (int j = i + 1; j < n_front; ++j) {
			if (sorted_front_strip[j] < sorted_front_strip[i]) {
				std::swap(sorted_front_strip[i], sorted_front_strip[j]);
				std::swap(sorted_front_e[i], sorted_front_e[j]);
				std::swap(sorted_front_t[i], sorted_front_t[j]);
			}
		}
	}
	front_strip = sorted_front_strip.data();
	front_e = sorted_front_e.data();
	front_t = sorted_front_t.data();

	// sort back hits by strip number (ascending)
	std::vector<int> sorted_back_strip(back_strip, back_strip + n_back);
	std::vector<double> sorted_back_e(back_e, back_e + n_back);
	std::vector<double> sorted_back_t(back_t, back_t + n_back);
	for (int i = 0; i < n_back - 1; ++i) {
		for (int j = i + 1; j < n_back; ++j) {
			if (sorted_back_strip[j] < sorted_back_strip[i]) {
				std::swap(sorted_back_strip[i], sorted_back_strip[j]);
				std::swap(sorted_back_e[i], sorted_back_e[j]);
				std::swap(sorted_back_t[i], sorted_back_t[j]);
			}
		}
	}
	back_strip = sorted_back_strip.data();
	back_e = sorted_back_e.data();
	back_t = sorted_back_t.data();

	// 1f-1b matching
	if (n_front == 1 && n_back == 1) {
		if (std::fabs(front_e[0] - back_e[0]) < match_tolerance) {
			double diff = std::fabs(front_e[0] - back_e[0]);
			if (h_energy_diff) h_energy_diff->Fill(diff);
			AppendMatchResult(output, detector,
				double(front_strip[0]), double(back_strip[0]),
				front_e[0], front_t[0], 0, diff);
			return true;
		}
		return false;
	}

	// 2f-1b matching
	if (n_front == 2 && n_back == 1) {
		int best_mode = 0;
		double min_diff = std::fabs(front_e[0] - back_e[0]);
		if (std::fabs(front_e[1] - back_e[0]) < min_diff) {
			min_diff = std::fabs(front_e[1] - back_e[0]);
			best_mode = 1;
		}
		if (std::fabs((front_e[0] + front_e[1]) - back_e[0]) < min_diff) {
			min_diff = std::fabs((front_e[0] + front_e[1]) - back_e[0]);
			best_mode = 2;
		}

		if (min_diff < match_tolerance) {
			if (h_energy_diff) h_energy_diff->Fill(min_diff);
			if (best_mode == 0 || best_mode == 1) {
				AppendMatchResult(output, detector,
					double(front_strip[best_mode]), double(back_strip[0]),
					front_e[best_mode], front_t[best_mode], 0, min_diff);
			} else if (best_mode == 2) {
				bool adjacent = std::abs(front_strip[0] - front_strip[1]) == 1;
				if (adjacent) {
					double time_val = front_t[0];
					if (front_e[1] > front_e[0]) time_val = front_t[1];
					AppendMatchResult(output, detector,
						WeightedStrip(front_strip[0], front_e[0], front_strip[1], front_e[1]),
						double(back_strip[0]),
						front_e[0] + front_e[1], time_val, 1, min_diff);
				} else {
				double frac0 = front_e[0] / (front_e[0] + front_e[1]);
				double frac1 = front_e[1] / (front_e[0] + front_e[1]);
				double e0 = back_e[0] * frac0;
				double e1 = back_e[0] * frac1;

				AppendMatchResult(output, detector,
					double(front_strip[0]), double(back_strip[0]),
					e0, front_t[0], 3, std::fabs(e0 - front_e[0]));
				AppendMatchResult(output, detector,
					double(front_strip[1]), double(back_strip[0]),
					e1, front_t[1], 3, std::fabs(e1 - front_e[1]));
			}
			}
			return true;
		}
		return false;
	}

	// 1f-2b matching
	if (n_front == 1 && n_back == 2) {
		int best_mode = 0;
		double min_diff = std::fabs(back_e[0] - front_e[0]);
		if (std::fabs(back_e[1] - front_e[0]) < min_diff) {
			min_diff = std::fabs(back_e[1] - front_e[0]);
			best_mode = 1;
		}
		if (std::fabs((back_e[0] + back_e[1]) - front_e[0]) < min_diff) {
			min_diff = std::fabs((back_e[0] + back_e[1]) - front_e[0]);
			best_mode = 2;
		}

		if (min_diff < match_tolerance) {
			if (h_energy_diff) h_energy_diff->Fill(min_diff);
			if (best_mode == 0 || best_mode == 1) {
				AppendMatchResult(output, detector,
					double(front_strip[0]), double(back_strip[best_mode]),
					front_e[0], front_t[0], 0, min_diff);
			} else if (best_mode == 2) {
				bool adjacent = std::abs(back_strip[0] - back_strip[1]) == 1;
				if (adjacent) {
					AppendMatchResult(output, detector,
						double(front_strip[0]),
						WeightedStrip(back_strip[0], back_e[0], back_strip[1], back_e[1]),
						front_e[0], front_t[0], 1, min_diff);
				} else {
					double frac_0 = back_e[0] / (back_e[0] + back_e[1]);
					double frac_1 = back_e[1] / (back_e[0] + back_e[1]);
					AppendMatchResult(output, detector,
						double(front_strip[0]), double(back_strip[0]),
						front_e[0] * frac_0, front_t[0], 3,
						std::fabs(front_e[0] * frac_0 - back_e[0]));
					AppendMatchResult(output, detector,
						double(front_strip[0]), double(back_strip[1]),
						front_e[0] * frac_1, front_t[0], 3,
						std::fabs(front_e[0] * frac_1 - back_e[1]));
				}
			}
			return true;
		}
		return false;
	}

	return false;
}

void MatchComplex(
	const DssdEvent &input,
	const SquareDetectorConfig &detector,
	DssdMatchEvent &output,
	TH1D *h_energy_diff
) {
	int n_front = input.front_num;
	int n_back = input.back_num;
	const int *front_strip = input.front_strip;
	const int *back_strip = input.back_strip;
	const double *front_e = input.front_energy;
	const double *back_e = input.back_energy;
	const double *front_t = input.front_time;
	const double *back_t = input.back_time;
	double match_tolerance = detector.match_tolerance;

	// sort front hits by strip number (ascending)
	std::vector<int> sorted_front_strip(front_strip, front_strip + n_front);
	std::vector<double> sorted_front_e(front_e, front_e + n_front);
	std::vector<double> sorted_front_t(front_t, front_t + n_front);
	for (int i = 0; i < n_front - 1; ++i) {
		for (int j = i + 1; j < n_front; ++j) {
			if (sorted_front_strip[j] < sorted_front_strip[i]) {
				std::swap(sorted_front_strip[i], sorted_front_strip[j]);
				std::swap(sorted_front_e[i], sorted_front_e[j]);
				std::swap(sorted_front_t[i], sorted_front_t[j]);
			}
		}
	}
	front_strip = sorted_front_strip.data();
	front_e = sorted_front_e.data();
	front_t = sorted_front_t.data();

	std::map<int, int> available_back;
	for (int i = 0; i < n_back; ++i) {
		available_back[back_strip[i]] = i;
	}

	int skip_next_front = 0;
	int front_matched = 0;

	for (int i = 0; i < n_front; ++i) {
		if (available_back.empty()) break;
		if (skip_next_front == 1) {
			skip_next_front = 0;
			continue;
		}

		std::map<int, double> candidate_energies;
		candidate_energies[0] = front_e[i];
		if (i > 0 && std::abs(front_strip[i] - front_strip[i - 1]) == 1 && front_matched == 0) {
			candidate_energies[1] = front_e[i] + front_e[i - 1];
		} else {
			candidate_energies[1] = 999999.0;
		}
		if (i < (n_front - 1) && std::abs(front_strip[i] - front_strip[i + 1]) == 1) {
			candidate_energies[2] = front_e[i] + front_e[i + 1];
		} else {
			candidate_energies[2] = 999999.0;
		}

		double min_diff = 9999.0;
		int front_merge_mode = 0;
		int back_merge_mode = 0;
		int back_strip_1 = 0;
		int back_strip_2 = 0;

		for (int ie = 0; ie < 3; ++ie) {
			if (candidate_energies[ie] >= 999998.0) continue;

			for (auto it = available_back.begin(); it != available_back.end(); ++it) {
				double diff_single = std::fabs(candidate_energies[ie] - back_e[it->second]);
				if (diff_single < min_diff) {
					min_diff = diff_single;
					front_merge_mode = ie;
					back_merge_mode = 0;
					back_strip_1 = it->first;
				}

				auto prev_it = it;
				if (prev_it != available_back.begin()) {
					--prev_it;
					if (std::abs(it->first - prev_it->first) == 1) {
						double diff_merged = std::fabs(
							candidate_energies[ie] - back_e[it->second] - back_e[prev_it->second]);
						if (diff_merged < min_diff) {
							min_diff = diff_merged;
							front_merge_mode = ie;
							back_merge_mode = 1;
							back_strip_2 = prev_it->first;
							back_strip_1 = back_strip_2 + 1;
						}
					}
				}

				auto next_it = it;
				++next_it;
				if (next_it != available_back.end()) {
					if (std::abs(next_it->first - it->first) == 1) {
						double diff_merged = std::fabs(
							candidate_energies[ie] - back_e[it->second] - back_e[next_it->second]);
						if (diff_merged < min_diff) {
							min_diff = diff_merged;
							front_merge_mode = ie;
							back_merge_mode = 2;
							back_strip_2 = next_it->first;
							back_strip_1 = back_strip_2 - 1;
						}
					}
				}
			}
		}

		front_matched = 0;
		if (min_diff < match_tolerance) {
			if (h_energy_diff) h_energy_diff->Fill(min_diff);
			front_matched = 1;

			if (front_merge_mode == 0 && back_merge_mode == 0) {
				auto it = available_back.find(back_strip_1);
				AppendMatchResult(output, detector,
					double(front_strip[i]), double(back_strip[it->second]),
					front_e[i], front_t[i], 0, min_diff);
				available_back.erase(back_strip_1);
			} else if (front_merge_mode == 0 && (back_merge_mode == 1 || back_merge_mode == 2)) {
				auto it1 = available_back.find(back_strip_1);
				auto it2 = available_back.find(back_strip_2);
				double eb1 = back_e[it1->second];
				double eb2 = back_e[it2->second];

				AppendMatchResult(output, detector,
					double(front_strip[i]),
					WeightedStrip(back_strip[it1->second], eb1, back_strip[it2->second], eb2),
					front_e[i], front_t[i], 1, min_diff);
				available_back.erase(back_strip_1);
				available_back.erase(back_strip_2);
			} else if ((front_merge_mode == 1 || front_merge_mode == 2) && back_merge_mode == 0) {
				int partner_idx = i - 1;
				if (front_merge_mode == 2) {
					partner_idx = i + 1;
					skip_next_front = 1;
				}
				auto it = available_back.find(back_strip_1);
				double ef1 = front_e[i];
				double ef2 = front_e[partner_idx];
				double time_val = front_t[i];
				if (ef2 > ef1) time_val = front_t[partner_idx];

				AppendMatchResult(output, detector,
					WeightedStrip(front_strip[i], ef1, front_strip[partner_idx], ef2),
					double(back_strip[it->second]),
					ef1 + ef2, time_val, 2, min_diff);
				available_back.erase(back_strip_1);
			} else if ((front_merge_mode == 1 || front_merge_mode == 2) && (back_merge_mode == 1 || back_merge_mode == 2)) {
				int partner_idx = i - 1;
				if (front_merge_mode == 2) {
					partner_idx = i + 1;
					skip_next_front = 1;
				}
				auto it1 = available_back.find(back_strip_1);
				auto it2 = available_back.find(back_strip_2);
				double eb1 = back_e[it1->second];
				double eb2 = back_e[it2->second];
				double ef1 = front_e[i];
				double ef2 = front_e[partner_idx];
				double time_val = front_t[i];
				if (ef2 > ef1) time_val = front_t[partner_idx];
				if (eb2 > eb1) time_val = back_t[it2->second];

				AppendMatchResult(output, detector,
					WeightedStrip(front_strip[i], ef1, front_strip[partner_idx], ef2),
					WeightedStrip(back_strip[it1->second], eb1, back_strip[it2->second], eb2),
					ef1 + ef2, time_val, 3, min_diff);
				available_back.erase(back_strip_1);
				available_back.erase(back_strip_2);
			}
		}

		candidate_energies.clear();
	}
	available_back.clear();
}

} // namespace

void MatchDssdEvent(
	const DssdEvent &input,
	const SquareDetectorConfig &detector,
	DssdMatchEvent &output,
	TH1D *h_energy_diff
) {
	Reset(output);

	int n_front = input.front_num;
	int n_back = input.back_num;

	if (n_front == 0 || n_back == 0) return;

	MatchSimple(input, detector, output, h_energy_diff);

	if (output.num == 0) {
		MatchComplex_v1(input, detector, output, h_energy_diff);
	}

	for (int i = 0; i < output.num - 1; ++i) {
		for (int j = i + 1; j < output.num; ++j) {
			if (output.energy[j] > output.energy[i]) {
				std::swap(output.front_strip[i], output.front_strip[j]);
				std::swap(output.back_strip[i], output.back_strip[j]);
				std::swap(output.energy[i], output.energy[j]);
				std::swap(output.time[i], output.time[j]);
				std::swap(output.merge_tag[i], output.merge_tag[j]);
				std::swap(output.energy_diff[i], output.energy_diff[j]);
				std::swap(output.x[i], output.x[j]);
				std::swap(output.y[i], output.y[j]);
				std::swap(output.z[i], output.z[j]);
			}
		}
	}
}

} // namespace brill