#include "include/10C+4He/10C_4He_event.h"

namespace brill {

void SetupInputC10He4(TTree *tree, C10He4Event &event) {
	tree->SetBranchAddress("run_number", &event.run_number);
	tree->SetBranchAddress("entry", &event.entry);
	tree->SetBranchAddress("e1_10C", &event.e1_10C);
	tree->SetBranchAddress("e2_10C", &event.e2_10C);
	tree->SetBranchAddress("e3_10C", &event.e3_10C);
	tree->SetBranchAddress("e1_4He", &event.e1_4He);
	tree->SetBranchAddress("e2_4He", &event.e2_4He);
	tree->SetBranchAddress("e3_4He", &event.e3_4He);
	tree->SetBranchAddress("e4_4He", &event.e4_4He);
	tree->SetBranchAddress("e5_4He", &event.e5_4He);
	tree->SetBranchAddress("is_14O", &event.is_14O);
	tree->SetBranchAddress("is_13N", &event.is_13N);
	tree->SetBranchAddress("is_12C", &event.is_12C);
	tree->SetBranchAddress("ppac_valid", &event.ppac_valid);
	tree->SetBranchAddress("target_x", &event.target_x);
	tree->SetBranchAddress("target_y", &event.target_y);
	tree->SetBranchAddress("dir_x", &event.dir_x);
	tree->SetBranchAddress("dir_y", &event.dir_y);
	tree->SetBranchAddress("t0d2_10C_x", &event.t0d2_10C_x);
	tree->SetBranchAddress("t0d2_10C_y", &event.t0d2_10C_y);
	tree->SetBranchAddress("t0d2_10C_z", &event.t0d2_10C_z);
	tree->SetBranchAddress("t0d2_4He_x", &event.t0d2_4He_x);
	tree->SetBranchAddress("t0d2_4He_y", &event.t0d2_4He_y);
	tree->SetBranchAddress("t0d2_4He_z", &event.t0d2_4He_z);
	tree->SetBranchAddress("theta_beam", &event.theta_beam);
	tree->SetBranchAddress("phi_beam", &event.phi_beam);
	tree->SetBranchAddress("theta_4He", &event.theta_4He);
	tree->SetBranchAddress("theta_10C", &event.theta_10C);
	tree->SetBranchAddress("opening_angle", &event.opening_angle);
}

void SetupOutputC10He4(TTree *tree, C10He4Event &event) {
	tree->Branch("run_number", &event.run_number, "run_number/I");
	tree->Branch("entry", &event.entry, "entry/L");
	tree->Branch("e1_10C", &event.e1_10C, "e1_10C/D");
	tree->Branch("e2_10C", &event.e2_10C, "e2_10C/D");
	tree->Branch("e3_10C", &event.e3_10C, "e3_10C/D");
	tree->Branch("e1_4He", &event.e1_4He, "e1_4He/D");
	tree->Branch("e2_4He", &event.e2_4He, "e2_4He/D");
	tree->Branch("e3_4He", &event.e3_4He, "e3_4He/D");
	tree->Branch("e4_4He", &event.e4_4He, "e4_4He/D");
	tree->Branch("e5_4He", &event.e5_4He, "e5_4He/D");
	tree->Branch("is_14O", &event.is_14O, "is_14O/O");
	tree->Branch("is_13N", &event.is_13N, "is_13N/O");
	tree->Branch("is_12C", &event.is_12C, "is_12C/O");
	tree->Branch("ppac_valid", &event.ppac_valid, "ppac_valid/O");
	tree->Branch("target_x", &event.target_x, "target_x/D");
	tree->Branch("target_y", &event.target_y, "target_y/D");
	tree->Branch("dir_x", &event.dir_x, "dir_x/D");
	tree->Branch("dir_y", &event.dir_y, "dir_y/D");
	tree->Branch("t0d2_10C_x", &event.t0d2_10C_x, "t0d2_10C_x/D");
	tree->Branch("t0d2_10C_y", &event.t0d2_10C_y, "t0d2_10C_y/D");
	tree->Branch("t0d2_10C_z", &event.t0d2_10C_z, "t0d2_10C_z/D");
	tree->Branch("t0d2_4He_x", &event.t0d2_4He_x, "t0d2_4He_x/D");
	tree->Branch("t0d2_4He_y", &event.t0d2_4He_y, "t0d2_4He_y/D");
	tree->Branch("t0d2_4He_z", &event.t0d2_4He_z, "t0d2_4He_z/D");
	tree->Branch("theta_beam", &event.theta_beam, "theta_beam/D");
	tree->Branch("phi_beam", &event.phi_beam, "phi_beam/D");
	tree->Branch("theta_4He", &event.theta_4He, "theta_4He/D");
	tree->Branch("theta_10C", &event.theta_10C, "theta_10C/D");
	tree->Branch("opening_angle", &event.opening_angle, "opening_angle/D");
}

} // namespace brill