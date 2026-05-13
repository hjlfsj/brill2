#pragma once

#include "include/analysis/config.h"

namespace glimmer {

double StripCenter(double size_mm, int strips, double strip);

void DssdStripToPosition(
	const SquareDetectorConfig &detector,
	double front_strip,
	double back_strip,
	double &x_mm,
	double &y_mm,
	double &z_mm
);

} // namespace glimmer
