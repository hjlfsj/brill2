#include "include/event/ppac/ppac_track.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <TF1.h>
#include <TGraph.h>

namespace brill {

int ReadPpacOffsetParams(const std::string &path, PpacOffsetParams params[3]) {
	std::ifstream file(path);
	if (!file.is_open()) {
		std::cerr << "Error: Cannot open offset file: " << path << "\n";
		return 1;
	}

	std::string line;
	std::getline(file, line);

	for (int i = 0; i < 3; ++i) {
		if (!std::getline(file, line)) {
			std::cerr << "Error: Missing line for PPAC" << (i + 1) << " in offset file.\n";
			return 1;
		}
		std::istringstream iss(line);
		int ppac_id;
		double chi2_ndf_x, chi2_ndf_y, anode_eff, x_cathode_eff, y_cathode_eff;
		iss >> ppac_id
		    >> params[i].delay_x_peak >> params[i].delay_x_sigma
		    >> params[i].delay_y_peak >> params[i].delay_y_sigma
		    >> params[i].position_x_p0 >> params[i].position_x_p1 >> chi2_ndf_x
		    >> params[i].position_y_p0 >> params[i].position_y_p1 >> chi2_ndf_y
		    >> anode_eff >> x_cathode_eff >> y_cathode_eff;
	}

	file.close();
	return 0;
}

static bool FitLinear(
	const std::vector<double> &z,
	const std::vector<double> &pos,
	double &intercept,
	double &slope
) {
	int n = z.size();
	if (n < 2) return false;

	TGraph gr(n, z.data(), pos.data());
	gr.Fit("pol1", "Q");
	TF1 *fit = gr.GetFunction("pol1");
	if (!fit) return false;

	intercept = fit->GetParameter(0);
	slope = fit->GetParameter(1);
	return true;
}

void TrackPpac(
	const PpacEvent &event,
	const PpacConfig &ppac_config,
	const PpacOffsetParams offset[3],
	PpacTrackEvent &track
) {
	track.valid = 0;
	track.x_used_ppac = 0;
	track.y_used_ppac = 0;
	track.target_x = 0.0;
	track.target_y = 0.0;
	track.dir_x = 0.0;
	track.dir_y = 0.0;
	for (int i = 0; i < kMaxPpac; ++i) {
		track.ppac_x[i] = NAN;
		track.ppac_y[i] = NAN;
	}

	for (int ppac_idx = 0; ppac_idx < 3; ++ppac_idx) {
		int ch_x1 = PpacChannelMap::ch_x1[ppac_idx];
		int ch_x2 = PpacChannelMap::ch_x2[ppac_idx];
		int ch_y1 = PpacChannelMap::ch_y1[ppac_idx];
		int ch_y2 = PpacChannelMap::ch_y2[ppac_idx];
		int ch_anode = PpacChannelMap::ch_anode[ppac_idx];

		bool anode_valid = event.valid[ch_anode] && event.time[ch_anode] > 0;
		bool x_valid = anode_valid && event.valid[ch_x1] && event.time[ch_x1] > 0
		               && event.valid[ch_x2] && event.time[ch_x2] > 0;
		bool y_valid = anode_valid && event.valid[ch_y1] && event.time[ch_y1] > 0
		               && event.valid[ch_y2] && event.time[ch_y2] > 0;

		if (x_valid) {
			double delay_x = event.time[ch_x1] + event.time[ch_x2] - 2 * event.time[ch_anode];
			if (std::fabs(delay_x - offset[ppac_idx].delay_x_peak) < 3 * offset[ppac_idx].delay_x_sigma) {
				double dtx = event.time[ch_x1] - event.time[ch_x2];
				track.ppac_x[ppac_idx] = offset[ppac_idx].position_x_p0 + offset[ppac_idx].position_x_p1 * dtx;
				track.x_used_ppac |= (1u << ppac_idx);
			}
		}

		if (y_valid) {
			double delay_y = event.time[ch_y1] + event.time[ch_y2] - 2 * event.time[ch_anode];
			if (std::fabs(delay_y - offset[ppac_idx].delay_y_peak) < 3 * offset[ppac_idx].delay_y_sigma) {
				double dty = event.time[ch_y1] - event.time[ch_y2];
				track.ppac_y[ppac_idx] = offset[ppac_idx].position_y_p0 + offset[ppac_idx].position_y_p1 * dty;
				track.y_used_ppac |= (1u << ppac_idx);
			}
		}
	}

	std::vector<double> z_x, pos_x, z_y, pos_y;
	for (int i = 0; i < 3; ++i) {
		if (track.x_used_ppac & (1u << i)) {
			z_x.push_back(ppac_config.z_mm[i]);
			pos_x.push_back(track.ppac_x[i]);
		}
		if (track.y_used_ppac & (1u << i)) {
			z_y.push_back(ppac_config.z_mm[i]);
			pos_y.push_back(track.ppac_y[i]);
		}
	}

	bool x_valid = false;
	if (z_x.size() >= 2) {
		x_valid = FitLinear(z_x, pos_x, track.target_x, track.dir_x);
	} else if (z_x.size() == 1 && (track.x_used_ppac & 1u)) {
		track.target_x = track.ppac_x[0];
		track.dir_x = 0.0;
		x_valid = true;
	}

	bool y_valid = false;
	if (z_y.size() >= 2) {
		y_valid = FitLinear(z_y, pos_y, track.target_y, track.dir_y);
	} else if (z_y.size() == 1 && (track.y_used_ppac & 1u)) {
		track.target_y = track.ppac_y[0];
		track.dir_y = 0.0;
		y_valid = true;
	}

	if (x_valid && y_valid) {
		track.valid = 1;
	}
}

} // namespace brill