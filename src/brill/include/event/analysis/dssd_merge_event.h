#pragma once

#include <string>

#include <TTree.h>

#include "include/analysis/config.h"

namespace glimmer {

struct DssdMergeEvent {
	int num = 0;
	int front_strip[kMaxDssdHit] = {0};
	int back_strip[kMaxDssdHit] = {0};
	double front_energy[kMaxDssdHit] = {0.0};
	double back_energy[kMaxDssdHit] = {0.0};
	double energy[kMaxDssdHit] = {0.0};
	double time[kMaxDssdHit] = {0.0};
	double x[kMaxDssdHit] = {0.0};
	double y[kMaxDssdHit] = {0.0};
	double z[kMaxDssdHit] = {0.0};
	int merge_tag[kMaxDssdHit] = {0};
};

void SetupInput(TTree *tree, DssdMergeEvent &event, const std::string &prefix = "");
void SetupOutput(TTree *tree, DssdMergeEvent &event);

} // namespace glimmer
