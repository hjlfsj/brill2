#include "include/event/analysis/ppac_track_event.h"

namespace glimmer {

void SetupInput(TTree *tree, PpacTrackEvent &event, const std::string &prefix) {
	tree->SetBranchAddress((prefix + "valid").c_str(), &event.valid);
	tree->SetBranchAddress((prefix + "used_ppac").c_str(), &event.used_ppac);
	tree->SetBranchAddress((prefix + "target_x").c_str(), &event.target_x);
	tree->SetBranchAddress((prefix + "target_y").c_str(), &event.target_y);
	tree->SetBranchAddress((prefix + "dir_x").c_str(), &event.dir_x);
	tree->SetBranchAddress((prefix + "dir_y").c_str(), &event.dir_y);
	tree->SetBranchAddress((prefix + "ppac_x").c_str(), event.ppac_x);
	tree->SetBranchAddress((prefix + "ppac_y").c_str(), event.ppac_y);
	tree->SetBranchAddress((prefix + "ppac_valid").c_str(), event.ppac_valid);
}

void SetupOutput(TTree *tree, PpacTrackEvent &event) {
	tree->Branch("valid", &event.valid, "valid/I");
	tree->Branch("used_ppac", &event.used_ppac, "used_ppac/i");
	tree->Branch("target_x", &event.target_x, "target_x/D");
	tree->Branch("target_y", &event.target_y, "target_y/D");
	tree->Branch("dir_x", &event.dir_x, "dir_x/D");
	tree->Branch("dir_y", &event.dir_y, "dir_y/D");
	tree->Branch("ppac_x", event.ppac_x, "ppac_x[3]/D");
	tree->Branch("ppac_y", event.ppac_y, "ppac_y[3]/D");
	tree->Branch("ppac_valid", event.ppac_valid, "ppac_valid[3]/I");
}

} // namespace glimmer
