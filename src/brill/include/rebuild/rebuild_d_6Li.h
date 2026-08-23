#pragma once

#include "include/d_6Li/d_6Li_event.h"

#include <string>

class TCutG;

namespace brill {

struct D6LiKinematics {
	double E_6Li = 0.0;
	double E_10C = 0.0;
};

D6LiKinematics ComputeKinematics(const D6LiEvent &event);

TCutG* LoadCutGFromFile(const std::string &filepath);

} // namespace brill