#ifndef __TRIGGER_EVENT_H__
#define __TRIGGER_EVENT_H__

#include <string>

#include <TTree.h>

namespace glimmer {

struct TriggerEvent {
	int flag;
	bool valid[6];
	double time[6];
	long long external_time[6];
};

void SetupInput(TTree *tree, TriggerEvent &event, const std::string &prefix = "");

void SetupOutput(TTree *tree, TriggerEvent &event);

}

#endif
