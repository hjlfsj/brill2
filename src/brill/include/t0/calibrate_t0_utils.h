#pragma once

#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include <TCutG.h>
#include <TGraph.h>
#include <TH2F.h>
#include <TString.h>

#include "include/config.h"
#include "include/energy_calculator/delta_energy_calculator.h"
#include "include/energy_calculator/range_energy_calculator.h"

namespace brill {

constexpr int kT0LayerCount = 5;

extern const char *kT0LayerNames[kT0LayerCount];

struct T0LayerPairInfo {
	const char *name;
	int layer;
	const char *graph_name;
};

extern const T0LayerPairInfo kT0LayerPairs[4];

struct T0ParticlePidInfo {
	int layer;
	int charge;
	int mass;
	double left;
	double right;
	double offset;
	double weight = 1.0;
};

extern const std::vector<T0ParticlePidInfo> kT0PidInfo;

struct ParticleIdentity {
	int charge = 0;
	int mass = 0;
};

ParticleIdentity ParseParticleName(const std::string &name);

void WriteCalibrationParameters(
	const char *path,
	const double *parameters,
	int layer_count
);

struct LoadedCut {
	const T0LayerPairInfo *pair = nullptr;
	ParticleIdentity particle;
	std::unique_ptr<TCutG> cut;
};

std::string FindPreCalibrationFile(
	const std::string &estimate_dir,
	const std::string &prefix,
	int run
);

int LoadCuts(
	const std::string &cut_dir,
	const std::string &cut_suffix,
	std::vector<LoadedCut> &cuts
);

const T0ParticlePidInfo *GetPidInfo(int layer, int charge, int mass);

int GetT0CalibrationRun(int run);

extern const std::vector<ParticleIdentity> kCommonParticles;

TGraph *GenerateTheoryCurve(
	const RangeEnergyCalculator &calculator,
	double t_shallow,
	double t_deep
);

} // namespace brill