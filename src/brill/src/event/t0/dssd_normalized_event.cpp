#include "include/event/analysis/dssd_normalized_event.h"

namespace glimmer {

void SetupInput(TTree *tree, DssdNormalizedEvent &event, const std::string &prefix) {
	tree->SetBranchAddress((prefix + "front_num").c_str(), &event.front_num);
	tree->SetBranchAddress((prefix + "front_strip").c_str(), event.front_strip);
	tree->SetBranchAddress((prefix + "front_energy").c_str(), event.front_energy);
	tree->SetBranchAddress((prefix + "front_time").c_str(), event.front_time);
	tree->SetBranchAddress((prefix + "back_num").c_str(), &event.back_num);
	tree->SetBranchAddress((prefix + "back_strip").c_str(), event.back_strip);
	tree->SetBranchAddress((prefix + "back_energy").c_str(), event.back_energy);
	tree->SetBranchAddress((prefix + "back_time").c_str(), event.back_time);
}

void SetupOutput(TTree *tree, DssdNormalizedEvent &event) {
	tree->Branch("front_num", &event.front_num, "front_num/I");
	tree->Branch("front_strip", event.front_strip, "front_strip[front_num]/I");
	tree->Branch("front_energy", event.front_energy, "front_energy[front_num]/D");
	tree->Branch("front_time", event.front_time, "front_time[front_num]/D");
	tree->Branch("back_num", &event.back_num, "back_num/I");
	tree->Branch("back_strip", event.back_strip, "back_strip[back_num]/I");
	tree->Branch("back_energy", event.back_energy, "back_energy[back_num]/D");
	tree->Branch("back_time", event.back_time, "back_time[back_num]/D");
}

} // namespace glimmer
