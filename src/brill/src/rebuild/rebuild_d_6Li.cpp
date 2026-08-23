#include "include/rebuild/rebuild_d_6Li.h"

#include <TCutG.h>
#include <TROOT.h>

#include <stdexcept>

namespace brill {

D6LiKinematics ComputeKinematics(const D6LiEvent &event) {
	D6LiKinematics kin;
	kin.E_6Li = event.e1_6Li + event.e2_6Li;
	kin.E_10C = event.e1_10C + event.e2_10C + event.e3_10C + event.e4_10C;
	return kin;
}

TCutG* LoadCutGFromFile(const std::string &filepath) {
	gROOT->ProcessLine(TString::Format(".x %s", filepath.c_str()));

	std::string filename = filepath;
	size_t pos = filename.find_last_of("/\\");
	if (pos != std::string::npos) {
		filename = filename.substr(pos + 1);
	}
	pos = filename.find_last_of('.');
	if (pos != std::string::npos) {
		filename = filename.substr(0, pos);
	}

	TCutG *cut = dynamic_cast<TCutG*>(gROOT->FindObject(filename.c_str()));
	if (!cut) {
		cut = dynamic_cast<TCutG*>(gROOT->FindObject("cut"));
	}
	if (!cut) {
		throw std::runtime_error("Cannot find TCutG in " + filepath);
	}
	return cut;
}

} // namespace brill