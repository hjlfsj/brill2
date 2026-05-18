#ifndef __BRILL_DELTA_ENERGY_CALCULATOR_H__
#define __BRILL_DELTA_ENERGY_CALCULATOR_H__

#include <memory>
#include <string>
#include <vector>

#include <TSpline.h>

#include "include/config.h"

namespace brill {

class DeltaEnergyCalculator {
public:
	DeltaEnergyCalculator(
		const AppConfig &config,
		int charge,
		int mass
	);

	DeltaEnergyCalculator(
		const std::string &config_path,
		int charge,
		int mass
	);

	double Energy(unsigned short layer, double delta_energy) const;
	double DeltaEnergy(unsigned short layer, double energy) const;

	static int Initialize(
		const AppConfig &config,
		int charge,
		int mass
	);

	static int Initialize(
		const std::string &config_path,
		int charge,
		int mass
	);

private:
	AppConfig config_;
	int charge_ = 0;
	int mass_ = 0;
	std::vector<std::unique_ptr<TSpline3>> de_e_funcs_;
	std::vector<std::unique_ptr<TSpline3>> e_de_funcs_;

	int Load();
	std::string CachePath() const;
};

} // namespace brill

#endif
