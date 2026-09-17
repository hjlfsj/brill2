#include "include/t0/calibrate_t0_utils.h"

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TCutG.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>

#include "include/config.h"
#include "include/energy_calculator/range_energy_calculator.h"
#include "include/utils.h"

namespace brill {

const char *kT0LayerNames[kT0LayerCount] = {
	"t0d1", "t0d2", "t0d3", "t0d4", "t0s"
};

const T0LayerPairInfo kT0LayerPairs[4] = {
	{"d1d2", 0, "g_d1d2"},
	{"d2d3", 1, "g_d2d3"},
	{"d3d4", 2, "g_d3d4"},
	{"d4s",  3, "g_d4s"},
};

const std::vector<T0ParticlePidInfo> kT0PidInfo = {
	{0, 2,  4,  2250.0, 12000.0,       0},
	{0, 3,  6,  4355.0, 12000.0,   13000},

	{1, 2,  4,  3890.0,  8250.0,   26000},
	{1, 4,  7, 10930.0, 17320.0,   36000},
	{1, 6, 12, 22050.0, 45800.0,   55000},

	{2, 2,  4,  3400.0,  7100.0,  102000},
	{2, 4,  7,  9600.0, 18200.0,  110000},
	{2, 6, 12, 22050.0, 45800.0,  130000},

	{3, 1,  1,   960.0,  2300.0,  180000},
	{3, 2,  4,  3890.0,  9430.0,  184000},
	{3, 3,  6,  7550.0, 16600.0,  200000},
};

const std::vector<ParticleIdentity> kCommonParticles = {
	{1, 1}, {1, 2}, {1, 3},
	{2, 3}, {2, 4}, {2, 6},
	{3, 6}, {3, 7}, {3, 8}, {3, 9},
	{4, 7}, {4, 9}, {4, 10},
	{5, 10}, {5, 11}, {5, 12},
	{6, 10}, {6, 11}, {6, 12}, {6, 13},
	{7, 12}, {7, 13}, {7, 14},
	{8, 13}, {8, 14},
};

ParticleIdentity ParseParticleName(const std::string &name) {
	static const std::map<std::string, ParticleIdentity> kParticleMap = {
		{"1H",  {1, 1}},
		{"4He", {2, 4}},
		{"6Li", {3, 6}},
		{"7Be", {4, 7}},
		{"12C", {6, 12}},
	};
	auto it = kParticleMap.find(name);
	if (it != kParticleMap.end()) return it->second;
	return {0, 0};
}

void WriteCalibrationParameters(
	const char *path,
	const double *parameters,
	int layer_count
) {
	std::filesystem::path file_path(path);
	if (!file_path.parent_path().empty()) {
		std::filesystem::create_directories(file_path.parent_path());
	}
	std::ofstream fout(path);
	if (!fout.good()) {
		std::cerr << "Error: Open output parameter file "
			<< path << " failed.\n";
		return;
	}

	fout << "# layer p0 p1\n";
	for (int layer = 0; layer < layer_count; ++layer) {
		fout << layer << " "
			<< parameters[layer * 2] << " "
			<< parameters[layer * 2 + 1] << "\n";
	}
}

std::string FindPreCalibrationFile(
	const std::string &estimate_dir,
	const std::string &prefix,
	int run
) {
	for (const auto &entry :
		std::filesystem::directory_iterator(estimate_dir)) {
		if (!entry.is_regular_file()) continue;
		std::string name = entry.path().filename().string();
		char expected[256];
		std::snprintf(
			expected, sizeof(expected),
			"pre_calibration_%s%04d_",
			prefix.c_str(), run
		);
		if (name.find(expected) == 0
			&& name.find(".root") != std::string::npos) {
			return entry.path().string();
		}
	}
	return "";
}

int LoadCuts(
	const std::string &cut_dir,
	const std::string &cut_suffix,
	std::vector<LoadedCut> &cuts
) {
	for (const auto &pair : kT0LayerPairs) {
		for (const auto &entry :
			std::filesystem::directory_iterator(cut_dir)) {
			if (!entry.is_regular_file()) continue;
			std::string name = entry.path().filename().string();
			if (name.size() <= 2
				|| name.substr(name.size() - 2) != ".C") continue;
			std::string basename = name.substr(0, name.size() - 2);
			std::string expected_prefix =
				std::string(pair.name) + "_";
			if (basename.find(expected_prefix) != 0) continue;
			if (basename.find(cut_suffix) == std::string::npos)
				continue;

			std::string particle_name = basename.substr(
				expected_prefix.size(),
				basename.size() - expected_prefix.size()
					- cut_suffix.size()
			);
			ParticleIdentity pid =
				ParseParticleName(particle_name);
			if (pid.charge == 0 && pid.mass == 0) {
				std::cerr << "Warning: Unknown particle '"
					<< particle_name
					<< "' in cut file " << name
					<< ", skipping.\n";
				continue;
			}

			TString cut_object_name = TString::Format(
				"%s_%s%s",
				pair.name,
				particle_name.c_str(),
				cut_suffix.c_str()
			);
			std::string macro_path = entry.path().string();
			gROOT->ProcessLine(TString::Format(
				".x %s", macro_path.c_str()));
			TCutG *cutg = dynamic_cast<TCutG*>(
				gROOT->FindObject(cut_object_name));
			if (!cutg) {
				std::cerr << "Error: Load TCutG "
					<< cut_object_name
					<< " from " << macro_path
					<< " failed.\n";
				return -1;
			}

			LoadedCut loaded;
			loaded.pair = &pair;
			loaded.particle = pid;
			loaded.cut.reset(
				dynamic_cast<TCutG*>(cutg->Clone()));
			cuts.push_back(std::move(loaded));
			std::cout << "Loaded cut: " << name
				<< " (Z=" << pid.charge
				<< ", A=" << pid.mass << ")\n";
		}
	}
	return cuts.empty() ? -1 : 0;
}

const T0ParticlePidInfo *GetPidInfo(
	int layer,
	int charge,
	int mass
) {
	for (const auto &info : kT0PidInfo) {
		if (info.layer == layer
			&& info.charge == charge
			&& info.mass == mass) {
			return &info;
		}
	}
	return nullptr;
}

TGraph *GenerateTheoryCurve(
	const RangeEnergyCalculator &calculator,
	double t_shallow,
	double t_deep
) {
	auto *curve = new TGraph();
	curve->SetLineColor(kRed);
	curve->SetLineWidth(1);

	double E_min = calculator.Energy(t_shallow);
	double E_max = calculator.Energy(t_shallow + t_deep);

	constexpr double kStep = 0.1;
	for (double E = E_min; E <= E_max; E += kStep) {
		double r_residual = calculator.Range(E) - t_shallow;
		if (r_residual <= 0.0) continue;
		if (r_residual > t_deep) break;

		double E_residual = calculator.Energy(r_residual);
		double delta_E = E - E_residual;
		if (delta_E < 0.0 || E_residual < 0.0) continue;

		curve->SetPoint(curve->GetN(), E_residual, delta_E);
	}

	return curve;
}

} // namespace brill