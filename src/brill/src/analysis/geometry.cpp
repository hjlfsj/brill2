#include "include/analysis/geometry.h"

namespace glimmer {

double StripCenter(double size_mm, int strips, double strip) {
	return size_mm * ((strip + 0.5) / static_cast<double>(strips) - 0.5);
}

void DssdStripToPosition(
	const SquareDetectorConfig &detector,
	double front_strip,
	double back_strip,
	double &x_mm,
	double &y_mm,
	double &z_mm
) {
	x_mm = detector.center_x_mm + StripCenter(detector.size_x_mm, detector.back_strips, back_strip);
	y_mm = detector.center_y_mm + StripCenter(detector.size_y_mm, detector.front_strips, front_strip);
	z_mm = detector.z_mm;
}

} // namespace glimmer
