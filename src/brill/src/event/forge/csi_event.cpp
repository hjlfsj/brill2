#include "include/event/forge/csi_event.h"

namespace glimmer {

void SetupInput(TTree *tree, CsiEvent &event, const std::string &prefix) {
	tree->SetBranchAddress((prefix+"flag").c_str(), &event.flag);
	tree->SetBranchAddress((prefix+"valid").c_str(), event.valid);
	tree->SetBranchAddress((prefix+"time").c_str(), event.time);
	tree->SetBranchAddress((prefix+"energy").c_str(), event.energy);
}

void SetupOutput(TTree *tree, CsiEvent &event) {
	tree->Branch("flag", &event.flag, "flag/l");
	tree->Branch("valid", event.valid, "valid[36]/O");
	tree->Branch("time", event.time, "time[36]/D");
	tree->Branch("energy", event.energy, "energy[36]/I");
}

}
