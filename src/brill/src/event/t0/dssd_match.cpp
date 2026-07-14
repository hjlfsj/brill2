#include "include/event/t0/dssd_match.h"

#include <TH1D.h>

#include <algorithm>
#include <cmath>
#include <set>
#include <vector>

namespace brill {

double WeightedStrip(int strip0, double energy0, int strip1, double energy1) {
	if (energy0 * 0.1 > energy1) return double(strip0);
	return 0.5 * double(strip0 + strip1);
}

double StripPosition(double center, double size, int strips, double strip) {
	return center + size * ((strip + 0.5) / double(strips) - 0.5);
}

void FillPhysicalPosition(
	const SquareDetectorConfig &detector,
	double front_strip,
	double back_strip,
	double &x,
	double &y,
	double &z
) {
	if (detector.name == "t0d2") {
		x = StripPosition(
			detector.center_x_mm,
			detector.size_x_mm,
			detector.back_strips,
			back_strip
		);
		y = StripPosition(
			detector.center_y_mm,
			detector.size_y_mm,
			detector.front_strips,
			63-front_strip
		);
	} else {
		x = StripPosition(
			detector.center_x_mm,
			detector.size_x_mm,
			detector.front_strips,
			front_strip
		);
		y = StripPosition(
			detector.center_y_mm,
			detector.size_y_mm,
			detector.back_strips,
			back_strip
		);
	}
	z = detector.z_mm;
}

void AppendMatchResult(
	DssdMatchEvent &output,
	const SquareDetectorConfig &detector,
	double front_strip_val,
	double back_strip_val,
	double energy_val,
	double time_val,
	int merge_tag_val,
	double energy_diff_val
) {
	if (output.num >= 8) return;
	int index = output.num++;
	output.front_strip[index] = front_strip_val;
	output.back_strip[index] = back_strip_val;
	output.energy[index] = energy_val;
	output.time[index] = time_val;
	output.merge_tag[index] = merge_tag_val;
	output.energy_diff[index] = energy_diff_val;
	FillPhysicalPosition(
		detector,
		front_strip_val,
		back_strip_val,
		output.x[index],
		output.y[index],
		output.z[index]
	);
}

namespace {

struct SideCombo {
	int type = 0;
	int strip1 = -1;
	int strip2 = -1;
	double energy1 = 0.0;
	double energy2 = 0.0;
	double total_energy = 0.0;
	double time1 = 0.0;
	double time2 = 0.0;
	double time = 0.0;
	int idx1 = -1;
	int idx2 = -1;
};

struct DeEntry {
	double diff;
	int fi;
	int bi;
};

} // namespace

void MatchComplex_v1(
	const DssdEvent &input,
	const SquareDetectorConfig &detector,
	DssdMatchEvent &output,
	TH1D *h_energy_diff
) {
	int n_front = input.front_num;
	int n_back = input.back_num;
	const int *front_strip = input.front_strip;
	const int *back_strip = input.back_strip;
	const double *front_e = input.front_energy;
	const double *back_e = input.back_energy;
	const double *front_t = input.front_time;
	const double *back_t = input.back_time;
	double match_tolerance = detector.match_tolerance;

	if (n_front > 8) n_front = 8;
	if (n_back > 8) n_back = 8;

	std::vector<SideCombo> front_combos;
	std::vector<SideCombo> back_combos;

	auto build_combos = [](
		int n, const int *strip, const double *energy, const double *time,
		std::vector<SideCombo> &combos
	) {
		for (int i = 0; i < n; ++i) {
			SideCombo c;
			c.type = 0;
			c.strip1 = strip[i];
			c.energy1 = energy[i];
			c.total_energy = energy[i];
			c.time1 = time[i];
			c.time = time[i];
			c.idx1 = i;
			combos.push_back(c);
		}
		for (int i = 0; i < n; ++i) {
			for (int j = i + 1; j < n; ++j) {
				SideCombo c;
				c.type = (std::abs(strip[i] - strip[j]) == 1) ? 1 : 2;
				c.strip1 = strip[i];
				c.strip2 = strip[j];
				c.energy1 = energy[i];
				c.energy2 = energy[j];
				c.total_energy = energy[i] + energy[j];
				c.time1 = time[i];
				c.time2 = time[j];
				c.time = (energy[i] >= energy[j]) ? time[i] : time[j];
				c.idx1 = i;
				c.idx2 = j;
				combos.push_back(c);
			}
		}
	};

	build_combos(n_front, front_strip, front_e, front_t, front_combos);
	build_combos(n_back, back_strip, back_e, back_t, back_combos);

	std::vector<DeEntry> de_entries;
	for (size_t fi = 0; fi < front_combos.size(); ++fi) {
		for (size_t bi = 0; bi < back_combos.size(); ++bi) {
			double diff = std::fabs(
				front_combos[fi].total_energy - back_combos[bi].total_energy
			);
			if (diff < match_tolerance) {
				de_entries.push_back({diff, int(fi), int(bi)});
			}
		}
	}

	std::sort(de_entries.begin(), de_entries.end(),
		[](const DeEntry &a, const DeEntry &b) {
			return a.diff < b.diff;
		}
	);

	std::set<int> used_front_idx;
	std::set<int> used_back_idx;

	for (const DeEntry &entry : de_entries) {
		if (output.num >= 8) break;

		const SideCombo &fc = front_combos[entry.fi];
		const SideCombo &bc = back_combos[entry.bi];

		if (used_front_idx.count(fc.idx1)) continue;
		if (fc.idx2 >= 0 && used_front_idx.count(fc.idx2)) continue;
		if (used_back_idx.count(bc.idx1)) continue;
		if (bc.idx2 >= 0 && used_back_idx.count(bc.idx2)) continue;

		double fe = fc.total_energy;
		double be = bc.total_energy;
		double diff = entry.diff;

		if (!(fe > diff && be > diff)) continue;

		int match_type = -1;
		if (fc.type == 0 && bc.type == 0) {
			match_type = 0;
		} else if ((fc.type == 0 && bc.type == 1) || (fc.type == 1 && bc.type == 0)) {
			match_type = 1;
		} else if (fc.type == 1 && bc.type == 1) {
			match_type = 2;
		} else if ((fc.type == 0 && bc.type == 2) || (fc.type == 2 && bc.type == 0)) {
			match_type = 3;
		} else if ((fc.type == 1 && bc.type == 2) || (fc.type == 2 && bc.type == 1)) {
			match_type = 4;
		} else {
			continue;
		}

		if (h_energy_diff) h_energy_diff->Fill(diff);

		if (match_type == 0) {
			AppendMatchResult(output, detector,
				double(fc.strip1), double(bc.strip1),
				fe, fc.time, 0, diff);
		} else if (match_type == 1) {
			if (fc.type == 1) {
				double wstrip = WeightedStrip(fc.strip1, fc.energy1, fc.strip2, fc.energy2);
				AppendMatchResult(output, detector,
					wstrip, double(bc.strip1),
					fe, fc.time, 1, diff);
			} else {
				double wstrip = WeightedStrip(bc.strip1, bc.energy1, bc.strip2, bc.energy2);
				AppendMatchResult(output, detector,
					double(fc.strip1), wstrip,
					fe, fc.time, 1, diff);
			}
		} else if (match_type == 2) {
			double wfront = WeightedStrip(fc.strip1, fc.energy1, fc.strip2, fc.energy2);
			double wback = WeightedStrip(bc.strip1, bc.energy1, bc.strip2, bc.energy2);
			double time_val = fc.time;
			if (bc.energy1 > fc.energy1) time_val = bc.time;
			AppendMatchResult(output, detector,
				wfront, wback, fe, time_val, 2, diff);
		} else if (match_type == 3) {
			const SideCombo &nonadj = (fc.type == 2) ? fc : bc;
			const SideCombo &single = (fc.type == 0) ? fc : bc;
			bool is_front_nonadj = (fc.type == 2);

			if (output.num + 1 >= 8) {
				used_front_idx.insert(fc.idx1);
				if (fc.idx2 >= 0) used_front_idx.insert(fc.idx2);
				used_back_idx.insert(bc.idx1);
				if (bc.idx2 >= 0) used_back_idx.insert(bc.idx2);
				break;
			}

			double frac1 = nonadj.energy1 / nonadj.total_energy;
			double frac2 = nonadj.energy2 / nonadj.total_energy;

			double e1 = single.total_energy * frac1;
			double e2 = single.total_energy * frac2;
			double diff1 = std::fabs(e1 - nonadj.energy1);
			double diff2 = std::fabs(e2 - nonadj.energy2);

			if (is_front_nonadj) {
				AppendMatchResult(output, detector,
					double(nonadj.strip1), double(single.strip1),
					e1, nonadj.time1, 3, diff1);
				AppendMatchResult(output, detector,
					double(nonadj.strip2), double(single.strip1),
					e2, nonadj.time2, 3, diff2);
			} else {
				AppendMatchResult(output, detector,
					double(single.strip1), double(nonadj.strip1),
					e1, nonadj.time1, 3, diff1);
				AppendMatchResult(output, detector,
					double(single.strip1), double(nonadj.strip2),
					e2, nonadj.time2, 3, diff2);
			}
		} else if (match_type == 4) {
			const SideCombo &nonadj = (fc.type == 2) ? fc : bc;
			const SideCombo &adj = (fc.type == 1) ? fc : bc;
			bool is_front_nonadj = (fc.type == 2);

			if (output.num + 1 >= 8) {
				used_front_idx.insert(fc.idx1);
				if (fc.idx2 >= 0) used_front_idx.insert(fc.idx2);
				used_back_idx.insert(bc.idx1);
				if (bc.idx2 >= 0) used_back_idx.insert(bc.idx2);
				break;
			}

			double frac1 = nonadj.energy1 / nonadj.total_energy;
			double frac2 = nonadj.energy2 / nonadj.total_energy;

			double e1 = adj.total_energy * frac1;
			double e2 = adj.total_energy * frac2;
			double diff1 = std::fabs(e1 - nonadj.energy1);
			double diff2 = std::fabs(e2 - nonadj.energy2);

			double wadj = WeightedStrip(adj.strip1, adj.energy1, adj.strip2, adj.energy2);

			if (is_front_nonadj) {
				AppendMatchResult(output, detector,
					double(nonadj.strip1), wadj,
					e1, nonadj.time1, 4, diff1);
				AppendMatchResult(output, detector,
					double(nonadj.strip2), wadj,
					e2, nonadj.time2, 4, diff2);
			} else {
				AppendMatchResult(output, detector,
					wadj, double(nonadj.strip1),
					e1, nonadj.time1, 4, diff1);
				AppendMatchResult(output, detector,
					wadj, double(nonadj.strip2),
					e2, nonadj.time2, 4, diff2);
			}
		}

		used_front_idx.insert(fc.idx1);
		if (fc.idx2 >= 0) used_front_idx.insert(fc.idx2);
		used_back_idx.insert(bc.idx1);
		if (bc.idx2 >= 0) used_back_idx.insert(bc.idx2);
	}
}

} // namespace brill