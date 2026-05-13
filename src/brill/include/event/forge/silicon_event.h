#ifndef __SILICON_EVENT_H__
#define __SILICON_EVENT_H__

#include <string>

#include <TTree.h>

namespace glimmer {

struct SiliconEvent {
	bool valid;
	double time;
	int energy;
};

void SetupInput(TTree *tree, SiliconEvent &event, const std::string &prefix = "");

void SetupOutput(TTree *tree, SiliconEvent &event);

}

#endif
