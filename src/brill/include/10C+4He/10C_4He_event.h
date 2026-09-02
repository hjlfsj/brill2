#pragma once

#include <TTree.h>

namespace brill {

struct C10He4Event {
	int run_number = 0;
	Long64_t entry = 0;

	double e1_10C = 0.0;
	double e2_10C = 0.0;
	double e3_10C = 0.0;

	double e1_4He = 0.0;
	double e2_4He = 0.0;
	double e3_4He = 0.0;
	double e4_4He = 0.0;
	double e5_4He = 0.0;

	bool is_14O = false;
	bool is_13N = false;
	bool is_12C = false;

	bool ppac_valid = false;
	double target_x = 0.0;
	double target_y = 0.0;
	double dir_x = 0.0;
	double dir_y = 0.0;

	double t0d2_10C_x = 0.0;
	double t0d2_10C_y = 0.0;
	double t0d2_10C_z = 0.0;
	double t0d2_4He_x = 0.0;
	double t0d2_4He_y = 0.0;
	double t0d2_4He_z = 0.0;

	double theta_beam = 0.0;
	double phi_beam = 0.0;
	double theta_4He = 0.0;
	double theta_10C = 0.0;
	double opening_angle = 0.0;
};

void SetupInputC10He4(TTree *tree, C10He4Event &event);
void SetupOutputC10He4(TTree *tree, C10He4Event &event);

} // namespace brill