#pragma once

#include <string>

#include <TTree.h>

#include "include/analysis/config.h"

namespace glimmer {

struct DssdNormalizedEvent {
	int front_num = 0;
	int front_strip[kMaxDssdHit] = {0};
	double front_energy[kMaxDssdHit] = {0.0};
	double front_time[kMaxDssdHit] = {0.0};
	int back_num = 0;
	int back_strip[kMaxDssdHit] = {0};
	double back_energy[kMaxDssdHit] = {0.0};
	double back_time[kMaxDssdHit] = {0.0};
};

void SetupInput(TTree *tree, DssdNormalizedEvent &event, const std::string &prefix = "");
void SetupOutput(TTree *tree, DssdNormalizedEvent &event);

} // namespace glimmer
