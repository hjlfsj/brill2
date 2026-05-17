#ifndef __BRILL_LOST_ENERGY_CALCULATOR_H__
#define __BRILL_LOST_ENERGY_CALCULATOR_H__

#include <string>

#include <catima/catima.h>

#include "include/energy_calculator/range_energy_calculator.h"

namespace brill {

class LostEnergyCalculator {
public:
	LostEnergyCalculator(
		int charge,
		int mass,
		const catima::Material &material,
		double thickness_um,
		const std::string &cache_path
	);

	double ResidualEnergy(double incident_energy) const;
	double EnergyLost(double incident_energy) const;
	double IncidentEnergy(double residual_energy) const;
	double EnergyLossFromResidual(double residual_energy) const;

private:
	RangeEnergyCalculator calculator_;
	double thickness_um_ = 0.0;
};

} // namespace brill

#endif
