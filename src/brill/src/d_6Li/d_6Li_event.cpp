#include "include/d_6Li/d_6Li_event.h"

namespace brill {

void SetupInputD6Li(TTree *tree, D6LiEvent &event) {
	tree->SetBranchAddress("run_number", &event.run_number);
	tree->SetBranchAddress("entry", &event.entry);
	tree->SetBranchAddress("e1", &event.e1);
	tree->SetBranchAddress("e2", &event.e2);
	tree->SetBranchAddress("e3", &event.e3);
	tree->SetBranchAddress("e4", &event.e4);
	tree->SetBranchAddress("e1_10C", &event.e1_10C);
	tree->SetBranchAddress("e2_10C", &event.e2_10C);
	tree->SetBranchAddress("e3_10C", &event.e3_10C);
	tree->SetBranchAddress("e4_10C", &event.e4_10C);
	tree->SetBranchAddress("e1_6Li", &event.e1_6Li);
	tree->SetBranchAddress("e2_6Li", &event.e2_6Li);
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
	tree->SetBranchAddress("t0d2_6Li_x", &event.t0d2_6Li_x);
	tree->SetBranchAddress("t0d2_6Li_y", &event.t0d2_6Li_y);
	tree->SetBranchAddress("t0d2_6Li_z", &event.t0d2_6Li_z);
	tree->SetBranchAddress("theta_beam", &event.theta_beam);
	tree->SetBranchAddress("phi_beam", &event.phi_beam);
	tree->SetBranchAddress("theta_6Li", &event.theta_6Li);
	tree->SetBranchAddress("theta_10C", &event.theta_10C);
	tree->SetBranchAddress("opening_angle", &event.opening_angle);
}

void SetupOutputD6Li(TTree *tree, D6LiEvent &event) {
	tree->Branch("run_number", &event.run_number, "run_number/I");
	tree->Branch("entry", &event.entry, "entry/L");
	tree->Branch("e1", &event.e1, "e1/D");
	tree->Branch("e2", &event.e2, "e2/D");
	tree->Branch("e3", &event.e3, "e3/D");
	tree->Branch("e4", &event.e4, "e4/D");
	tree->Branch("e1_10C", &event.e1_10C, "e1_10C/D");
	tree->Branch("e2_10C", &event.e2_10C, "e2_10C/D");
	tree->Branch("e3_10C", &event.e3_10C, "e3_10C/D");
	tree->Branch("e4_10C", &event.e4_10C, "e4_10C/D");
	tree->Branch("e1_6Li", &event.e1_6Li, "e1_6Li/D");
	tree->Branch("e2_6Li", &event.e2_6Li, "e2_6Li/D");
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
	tree->Branch("t0d2_6Li_x", &event.t0d2_6Li_x, "t0d2_6Li_x/D");
	tree->Branch("t0d2_6Li_y", &event.t0d2_6Li_y, "t0d2_6Li_y/D");
	tree->Branch("t0d2_6Li_z", &event.t0d2_6Li_z, "t0d2_6Li_z/D");
	tree->Branch("theta_beam", &event.theta_beam, "theta_beam/D");
	tree->Branch("phi_beam", &event.phi_beam, "phi_beam/D");
	tree->Branch("theta_6Li", &event.theta_6Li, "theta_6Li/D");
	tree->Branch("theta_10C", &event.theta_10C, "theta_10C/D");
	tree->Branch("opening_angle", &event.opening_angle, "opening_angle/D");
}

} // namespace brill