#ifndef __BRILL_RANGE_ENERGY_CALCULATOR_H__
#define __BRILL_RANGE_ENERGY_CALCULATOR_H__

#include <memory>
#include <string>

#include <TSpline.h>
#include <catima/catima.h>

namespace brill {

class RangeEnergyCalculator {
public:
	RangeEnergyCalculator(
		int charge,
		int mass,
		const catima::Material &material,
		const std::string &cache_path
	);

	double Range(double energy) const;
	double Energy(double range) const;

private:
	int charge_ = 0;
	int mass_ = 0;
	catima::Material material_;
	std::string cache_path_;
	double max_energy_ = 0.0;
	double max_range_ = 0.0;
	std::unique_ptr<TSpline3> range_energy_func_;
	std::unique_ptr<TSpline3> energy_range_func_;

	bool LoadSplines();
	void BuildSplines();
};

catima::Material SiliconMaterial();

} // namespace brill

#endif
