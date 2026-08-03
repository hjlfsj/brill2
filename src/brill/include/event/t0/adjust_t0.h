#pragma once

#include <string>
#include <vector>

#include "include/config.h"

namespace brill {

struct T0AdjustResult {
	std::string detector_name;
	double dx = 0.0;
	double dy = 0.0;
	double dsigma_x = 0.0;
	double dsigma_y = 0.0;
	double residual_x_mean = 0.0;
	double residual_y_mean = 0.0;
	double residual_x_sigma = 0.0;
	double residual_y_sigma = 0.0;
	int num_events = 0;
};

struct T0EnergyCut {
	double min = 0.0;
	double max = 0.0;
};

struct T0AdjustConfig {
	std::string track_path;
	std::string match_dir;
	std::string trigger;
	int run = 0;
	std::string output_dir;
	std::vector<std::string> detector_names;
	std::vector<T0EnergyCut> energy_cuts;
};

int AdjustT0Step1(
	const T0AdjustConfig &adjust_config,
	const std::vector<SquareDetectorConfig> &detectors,
	std::vector<T0AdjustResult> &results
);

int AdjustT0Step2(
	const T0AdjustConfig &adjust_config,
	const std::vector<SquareDetectorConfig> &detectors,
	std::vector<T0AdjustResult> &results
);

} // namespace brill