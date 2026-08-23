#pragma once

#include <TTree.h>

namespace brill {

struct D6LiEvent {
	int run_number = 0;
	Long64_t entry = 0;

	double e1 = 0.0;
	double e2 = 0.0;
	double e3 = 0.0;
	double e4 = 0.0;

	double e1_10C = 0.0;
	double e2_10C = 0.0;
	double e3_10C = 0.0;
	double e4_10C = 0.0;

	double e1_6Li = 0.0;
	double e2_6Li = 0.0;

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
	double t0d2_6Li_x = 0.0;
	double t0d2_6Li_y = 0.0;
	double t0d2_6Li_z = 0.0;

	double theta_beam = 0.0;
	double phi_beam = 0.0;
	double theta_6Li = 0.0;
	double theta_10C = 0.0;
	double opening_angle = 0.0;
};

void SetupInputD6Li(TTree *tree, D6LiEvent &event);
void SetupOutputD6Li(TTree *tree, D6LiEvent &event);

} // namespace brill