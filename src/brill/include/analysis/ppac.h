#pragma once

#include <string>

#include "include/analysis/config.h"
#include "include/event/analysis/ppac_hit_event.h"
#include "include/event/analysis/ppac_track_event.h"
#include "include/event/forge/ppac_event.h"

namespace glimmer {

struct PpacNormalizeParameters {
	double x_offset[kMaxPpac] = {0.0};
	double y_offset[kMaxPpac] = {0.0};
	double x_scale[kMaxPpac] = {1.0, 1.0, 1.0};
	double y_scale[kMaxPpac] = {1.0, 1.0, 1.0};
};

bool ComputePpacRawPosition(
	const PpacEvent &event,
	int index,
	double &x,
	double &y
);

void BuildPpacHitEvent(
	const PpacEvent &event,
	const PpacNormalizeParameters &parameters,
	PpacHitEvent &hit
);

bool FitPpacTrace(
	const double *z,
	const double *value,
	const bool *valid,
	double &intercept,
	double &slope
);

int RunPpacNormalize(const AppConfig &config, const std::string &trigger, int run, int end_run);
int RunPpacTrack(const AppConfig &config, const std::string &trigger, int run);
int ReadPpacNormalizeParameters(const std::string &path, PpacNormalizeParameters &parameters);

} // namespace glimmer
