#include "include/event/analysis/t0_track_event.h"

namespace glimmer {

void SetupInput(TTree *tree, T0TrackEvent &event, const std::string &prefix) {
	tree->SetBranchAddress((prefix + "num").c_str(), &event.num);
	tree->SetBranchAddress((prefix + "flag").c_str(), event.flag);
	tree->SetBranchAddress((prefix + "energy").c_str(), event.energy);
	tree->SetBranchAddress((prefix + "x").c_str(), event.x);
	tree->SetBranchAddress((prefix + "y").c_str(), event.y);
	tree->SetBranchAddress((prefix + "z").c_str(), event.z);
}

void SetupOutput(TTree *tree, T0TrackEvent &event) {
	tree->Branch("num", &event.num, "num/I");
	tree->Branch("flag", event.flag, "flag[num]/i");
	tree->Branch("energy", event.energy, "energy[num][4]/D");
	tree->Branch("x", event.x, "x[num][4]/D");
	tree->Branch("y", event.y, "y[num][4]/D");
	tree->Branch("z", event.z, "z[num][4]/D");
}

} // namespace glimmer
