#pragma once

#include <string>

#include "include/config.h"
#include "include/event/ppac/ppac_track_event.h"

namespace brill {

struct PpacPositionAdjust {
	double dx2 = 0.0;
	double dx3 = 0.0;
	double dy2 = 0.0;
	double dy3 = 0.0;
	double dsigma_x2 = 0.0;
	double dsigma_x3 = 0.0;
	double dsigma_y2 = 0.0;
	double dsigma_y3 = 0.0;
	double baseline_residual_x = 0.0;
	double best_residual_x = 0.0;
	double baseline_residual_y = 0.0;
	double best_residual_y = 0.0;
	double improvement_x = 0.0;
	double improvement_y = 0.0;
};

int AdjustPpacPosition(
	const std::string &track_path,
	const PpacConfig &ppac_config,
	int n,
	PpacPositionAdjust &result
);

} // namespace brill