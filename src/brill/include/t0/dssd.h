#pragma once

#include <string>

#include "include/config.h"
#include "include/event/forge/dssd_event.h"

namespace brill {

struct DssdNormalizeParameters {
	int front_strips = 0;
	int back_strips = 0;
	double front_p0[kMaxStrips] = {0.0};
	double front_p1[kMaxStrips] = {0.0};
	double front_p2[kMaxStrips] = {0.0};
	double back_p0[kMaxStrips] = {0.0};
	double back_p1[kMaxStrips] = {0.0};
	double back_p2[kMaxStrips] = {0.0};
};

// int RunDssdNormalize(
// 	const AppConfig &config,
// 	const std::string &detector,
// 	const std::string &trigger,
// 	int run,
// 	int end_run
// );

int WriteDssdNormalizeParameters(
	const std::string &front_path,
	const std::string &back_path,
	const DssdNormalizeParameters &parameters
);

int ReadDssdNormalizeParameters(
	const std::string &front_path,
	const std::string &back_path,
	DssdNormalizeParameters &parameters
);

void ApplyDssdNormalize(
	const DssdEvent &input,
	const DssdNormalizeParameters &parameters,
	DssdEvent &output
);

// void MatchDssdEvent(
// 	const DssdNormalizedEvent &input,
// 	const SquareDetectorConfig &detector,
// 	DssdMergeEvent &output
// );

// int RunDssdMerge(
// 	const AppConfig &config,
// 	const std::string &detector,
// 	const std::string &trigger,
// 	int run
// );

// void FillDssdPidHistogram(const DssdMergeEvent &left, const DssdMergeEvent &right, TH2F &hist);

} // namespace brill
