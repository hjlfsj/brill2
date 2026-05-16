#ifndef __ALIGN_EVENT_H__
#define __ALIGN_EVENT_H__

#include <TTree.h>
#include <string>

namespace glimmer {

struct AlignEvent {
	// entry after forged
	long long xia_entry;
	// entry after converted
	long long vme_entry;
	// in nanoseconds
	long long xia_time;
	long long vme_time;
};

void SetupInput(TTree *tree, AlignEvent &event, const std::string &prefix = "");

void SetupOutput(TTree *tree, AlignEvent &event);


}

#endif