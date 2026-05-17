#include "include/energy_calculator/lost_energy_calculator.h"

#include <algorithm>

namespace brill {

LostEnergyCalculator::LostEnergyCalculator(
	int charge,
	int mass,
	const catima::Material &material,
	double thickness_um,
	const std::string &cache_path
)
: calculator_(charge, mass, material, cache_path)
, thickness_um_(thickness_um)
{}

double LostEnergyCalculator::ResidualEnergy(double incident_energy) const {
	if (incident_energy <= 0.0) return 0.0;
	double incident_range = calculator_.Range(incident_energy);
	if (incident_range <= thickness_um_) return 0.0;
	return std::max(0.0, calculator_.Energy(incident_range - thickness_um_));
}

double LostEnergyCalculator::EnergyLost(double incident_energy) const {
	if (incident_energy <= 0.0) return 0.0;
	return std::max(0.0, incident_energy - ResidualEnergy(incident_energy));
}

double LostEnergyCalculator::IncidentEnergy(double residual_energy) const {
	if (residual_energy <= 0.0) {
		return std::max(0.0, calculator_.Energy(thickness_um_));
	}
	return std::max(0.0, calculator_.Energy(calculator_.Range(residual_energy) + thickness_um_));
}

double LostEnergyCalculator::EnergyLossFromResidual(double residual_energy) const {
	if (residual_energy <= 0.0) return 0.0;
	return std::max(0.0, IncidentEnergy(residual_energy) - residual_energy);
}

} // namespace brill
