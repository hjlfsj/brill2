#pragma once

#include "include/config.h"
#include "include/event/ingot/dssd_event.h"
#include "include/event/t0/dssd_match_event.h"

class TH1D;

namespace brill {

double WeightedStrip(int strip0, double energy0, int strip1, double energy1);

double StripPosition(double center, double size, int strips, double strip);

void FillPhysicalPosition(
	const SquareDetectorConfig &detector,
	double front_strip,
	double back_strip,
	double &x,
	double &y,
	double &z
);

void AppendMatchResult(
	DssdMatchEvent &output,
	const SquareDetectorConfig &detector,
	double front_strip_val,
	double back_strip_val,
	double energy_val,
	double time_val,
	int merge_tag_val,
	double energy_diff_val
);

void MatchComplex_v1(
	const DssdEvent &input,
	const SquareDetectorConfig &detector,
	DssdMatchEvent &output,
	TH1D *h_energy_diff = nullptr
);

} // namespace brill