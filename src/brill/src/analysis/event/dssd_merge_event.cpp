#include "include/event/analysis/dssd_merge_event.h"

namespace glimmer {

void SetupInput(TTree *tree, DssdMergeEvent &event, const std::string &prefix) {
	tree->SetBranchAddress((prefix + "num").c_str(), &event.num);
	tree->SetBranchAddress((prefix + "front_strip").c_str(), event.front_strip);
	tree->SetBranchAddress((prefix + "back_strip").c_str(), event.back_strip);
	tree->SetBranchAddress((prefix + "front_energy").c_str(), event.front_energy);
	tree->SetBranchAddress((prefix + "back_energy").c_str(), event.back_energy);
	tree->SetBranchAddress((prefix + "energy").c_str(), event.energy);
	tree->SetBranchAddress((prefix + "time").c_str(), event.time);
	tree->SetBranchAddress((prefix + "x").c_str(), event.x);
	tree->SetBranchAddress((prefix + "y").c_str(), event.y);
	tree->SetBranchAddress((prefix + "z").c_str(), event.z);
	tree->SetBranchAddress((prefix + "merge_tag").c_str(), event.merge_tag);
}

void SetupOutput(TTree *tree, DssdMergeEvent &event) {
	tree->Branch("num", &event.num, "num/I");
	tree->Branch("front_strip", event.front_strip, "front_strip[num]/I");
	tree->Branch("back_strip", event.back_strip, "back_strip[num]/I");
	tree->Branch("front_energy", event.front_energy, "front_energy[num]/D");
	tree->Branch("back_energy", event.back_energy, "back_energy[num]/D");
	tree->Branch("energy", event.energy, "energy[num]/D");
	tree->Branch("time", event.time, "time[num]/D");
	tree->Branch("x", event.x, "x[num]/D");
	tree->Branch("y", event.y, "y[num]/D");
	tree->Branch("z", event.z, "z[num]/D");
	tree->Branch("merge_tag", event.merge_tag, "merge_tag[num]/I");
}

} // namespace glimmer
