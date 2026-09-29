#include "include/Lise++/target_energy_loss.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>

#include <catima/catima.h>
#include <TString.h>

#include "include/energy_calculator/range_energy_calculator.h"
#include "include/utils.h"

namespace brill {

static catima::Material CreateTargetMaterial(TargetType target) {
	catima::Material material;
	if (target == TargetType::kCD2) {
		material.add_element(12.0, 6.0, 1.0);
		material.add_element(2.0, 1.0, 2.0);
	} else {
		material.add_element(12.0, 6.0, 1.0);
		material.add_element(1.0, 1.0, 2.0);
	}
	material.density(0.79867);
	return material;
}

static double HalfThicknessUm(TargetType target) {
	if (target == TargetType::kCD2) {
		return 122.0 / 2.0;
	} else {
		return 180.0 / 2.0;
	}
}

static const char *TargetName(TargetType target) {
	return (target == TargetType::kCD2) ? "cd2" : "ch2";
}

TargetEnergyLoss::TargetEnergyLoss(
	const AppConfig &config,
	int charge,
	int mass
)
: charge_(charge)
, mass_(mass)
, cache_dir_(JoinPath(config.workspace, config.paths.energy_calculator))
{}

double TargetEnergyLoss::IncidentEnergy(
	double E,
	double theta,
	TargetType target,
	double &effective_half_thickness_um
) const {
	if (theta < 0.0 || theta > 70.0) {
		throw std::runtime_error("theta must be in [0, 70] degrees.");
	}

	double half_thickness_um = HalfThicknessUm(target);
	double cos_theta = std::cos(theta * M_PI / 180.0);
	effective_half_thickness_um = half_thickness_um / cos_theta;

	if (E <= 0.0) {
		return 0.0;
	}

	catima::Material material = CreateTargetMaterial(target);
	TString cache_path = TString::Format(
		"%s/target_z%d_a%d_%s.root",
		cache_dir_.c_str(),
		charge_,
		mass_,
		TargetName(target)
	);

	RangeEnergyCalculator calculator(charge_, mass_, material, cache_path.Data());

	double current_range = calculator.Range(E);
	double total_range = current_range + effective_half_thickness_um;
	double E0 = calculator.Energy(total_range);

	return E0;
}

} // namespace brill