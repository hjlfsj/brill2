#pragma once

#include <string>

#include <TTree.h>

namespace glimmer {

struct DssdEvent {
	int front_num;
	int front_strip[8];
	double front_energy[8];
	double front_time[8];
	int back_num;
	int back_strip[8];
	double back_energy[8];
	double back_time[8];
};

void SetupInput(TTree *tree, DssdEvent &event, const std::string &prefix);

void SetupOutput(TTree *tree, DssdEvent &event);

} // namespace glimmer
