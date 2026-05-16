#include "include/event/forge/trigger_event.h"

namespace glimmer {

void SetupInput(TTree *tree, TriggerEvent &event, const std::string &prefix) {
	tree->SetBranchAddress((prefix+"flag").c_str(), &event.flag);
	tree->SetBranchAddress((prefix+"valid").c_str(), event.valid);
	tree->SetBranchAddress((prefix+"time").c_str(), event.time);
	tree->SetBranchAddress((prefix+"external_time").c_str(), event.external_time);
}

void SetupOutput(TTree *tree, TriggerEvent &event) {
	tree->Branch("flag", &event.flag, "flag/I");
	tree->Branch("valid", event.valid, "valid[6]/O");
	tree->Branch("time", event.time, "time[6]/D");
	tree->Branch("external_time", event.external_time, "ets[6]/L");
}

}
