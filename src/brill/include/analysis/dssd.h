#pragma once

#include <string>

#include "include/analysis/config.h"
#include "include/event/analysis/dssd_merge_event.h"
#include "include/event/analysis/dssd_normalized_event.h"
#include "include/event/forge/dssd_event.h"

class TH2F;

namespace glimmer {

struct DssdNormalizeParameters {
	int front_strips = 0;
	int back_strips = 0;
	double front_offset[kMaxStrips] = {0.0};
	double front_scale[kMaxStrips] = {0.0};
	double back_offset[kMaxStrips] = {0.0};
	double back_scale[kMaxStrips] = {0.0};
};

int RunDssdNormalize(
	const AppConfig &config,
	const std::string &detector,
	const std::string &trigger,
	int run,
	int end_run
);

int ReadDssdNormalizeParameters(
	const std::string &path,
	DssdNormalizeParameters &parameters
);

void ApplyDssdNormalize(
	const DssdEvent &input,
	const DssdNormalizeParameters &parameters,
	DssdNormalizedEvent &output
);

void MergeDssdEvent(
	const DssdNormalizedEvent &input,
	const SquareDetectorConfig &detector,
	DssdMergeEvent &output
);

int RunDssdMerge(
	const AppConfig &config,
	const std::string &detector,
	const std::string &trigger,
	int run
);

void FillDssdPidHistogram(const DssdMergeEvent &left, const DssdMergeEvent &right, TH2F &hist);

} // namespace glimmer
