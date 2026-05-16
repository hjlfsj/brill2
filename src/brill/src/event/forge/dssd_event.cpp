#include "include/event/forge/dssd_event.h"

namespace brill {

void SetupInput(TTree *tree, DssdEvent &event, const std::string &prefix) {
	tree->SetBranchAddress((prefix+"front_num").c_str(), &event.front_num);
	tree->SetBranchAddress((prefix+"front_strip").c_str(), event.front_strip);
	tree->SetBranchAddress((prefix+"front_energy").c_str(), event.front_energy);
	tree->SetBranchAddress((prefix+"front_time").c_str(), event.front_time);
	tree->SetBranchAddress((prefix+"back_num").c_str(), &event.back_num);
	tree->SetBranchAddress((prefix+"back_strip").c_str(), event.back_strip);
	tree->SetBranchAddress((prefix+"back_energy").c_str(), event.back_energy);
	tree->SetBranchAddress((prefix+"back_time").c_str(), event.back_time);
}

void SetupOutput(TTree *tree, DssdEvent &event) {
	tree->Branch("front_num", &event.front_num, "fn/I");
	tree->Branch("front_strip", event.front_strip, "fs[fn]/I");
	tree->Branch("front_energy", event.front_energy, "fe[fn]/D");
	tree->Branch("front_time", event.front_time, "ft[fn]/D");
	tree->Branch("back_num", &event.back_num, "bn/I");
	tree->Branch("back_strip", event.back_strip, "bs[bn]/I");
	tree->Branch("back_energy", event.back_energy, "be[bn]/D");
	tree->Branch("back_time", event.back_time, "bt[bn]/D");
}

void Reset(DssdEvent &event) {
	event.front_num = 0;
	event.back_num = 0;
}

void Update(
	DssdEvent &event,
	bool front,
	int strip,
	double energy,
	double time
) {
	if (front) {
		event.front_strip[event.front_num] = strip;
		event.front_energy[event.front_num] = energy;
		event.front_time[event.front_num] = time;
		++event.front_num;
	} {
		event.back_strip[event.back_num] = strip;
		event.back_energy[event.back_num] = energy;
		event.back_time[event.back_num] = time;
		++event.back_num;
	}
}

}
