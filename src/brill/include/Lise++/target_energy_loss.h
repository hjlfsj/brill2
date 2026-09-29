#ifndef __BRILL_TARGET_ENERGY_LOSS_H__
#define __BRILL_TARGET_ENERGY_LOSS_H__

#include <memory>
#include <string>

#include "include/config.h"

namespace brill {

class RangeEnergyCalculator;

enum class TargetType {
	kCD2,
	kCH2
};

class TargetEnergyLoss {
public:
	TargetEnergyLoss(const AppConfig &config, int charge, int mass);

	double IncidentEnergy(
		double E, double theta, TargetType target,
		double &effective_half_thickness_um
	) const;

private:
	int charge_;
	int mass_;
	std::string cache_dir_;
};

} // namespace brill

#endif