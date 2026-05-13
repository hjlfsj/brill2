#include "include/event/forge/silicon_event.h"

namespace glimmer {

void SetupInput(TTree *tree, SiliconEvent &event, const std::string &prefix) {
	tree->SetBranchAddress((prefix+"valid").c_str(), &event.valid);
	tree->SetBranchAddress((prefix+"time").c_str(), &event.time);
	tree->SetBranchAddress((prefix+"energy").c_str(), &event.energy);
}

void SetupOutput(TTree *tree, SiliconEvent &event) {
	tree->Branch("valid", &event.valid, "valid/O");
	tree->Branch("time", &event.time, "time/D");
	tree->Branch("energy", &event.energy, "energy/I");
}

}
