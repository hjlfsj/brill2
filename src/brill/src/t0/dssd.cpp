#include "include/t0/dssd.h"

#include <algorithm>
#include <fstream>
#include <iostream>
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
	fout << "index p0 p1 p2\n";
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

} // namespace brill
