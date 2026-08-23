#include "include/physics/kinematics.h"

#include <cmath>

namespace brill {

double AngleBetween(const double v1[3], const double v2[3]) {
	double dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
	double mag1 = std::sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
	double mag2 = std::sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);
	if (mag1 == 0.0 || mag2 == 0.0) return 0.0;
	double cos_theta = dot / (mag1 * mag2);
	if (cos_theta > 1.0) cos_theta = 1.0;
	if (cos_theta < -1.0) cos_theta = -1.0;
	return std::acos(cos_theta) * 180.0 / M_PI;
}

double AngleWithZ(const double v[3]) {
	double mag = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
	if (mag == 0.0) return 0.0;
	double cos_theta = v[2] / mag;
	if (cos_theta > 1.0) cos_theta = 1.0;
	if (cos_theta < -1.0) cos_theta = -1.0;
	return std::acos(cos_theta) * 180.0 / M_PI;
}

} // namespace brill