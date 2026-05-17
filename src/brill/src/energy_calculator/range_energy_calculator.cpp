#include "include/energy_calculator/range_energy_calculator.h"

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <TFile.h>

namespace {

constexpr int kSplinePoints = 1200;
constexpr double kRangeEpsilon = 1e-6;
constexpr double kSiliconDensity = 2.329;

double RangeGcm2ToUm(double range_g_cm2, const catima::Material &material) {
	double density = material.density();
	if (density <= 0.0) {
		throw std::runtime_error("Invalid material density.");
	}
	return range_g_cm2 / density * 1e4;
}

} // namespace

namespace brill {

RangeEnergyCalculator::RangeEnergyCalculator(
	int charge,
	int mass,
	const catima::Material &material,
	const std::string &cache_path
)
: charge_(charge)
, mass_(mass)
, material_(material)
, cache_path_(cache_path)
, max_energy_(double(mass) * 500.0)
{
	if (charge_ <= 0 || mass_ <= 0) {
		throw std::runtime_error("Invalid projectile for RangeEnergyCalculator.");
	}
	if (!LoadSplines()) {
		BuildSplines();
	}
}

bool RangeEnergyCalculator::LoadSplines() {
	std::unique_ptr<TFile> file(TFile::Open(cache_path_.c_str(), "read"));
	if (!file || file->IsZombie()) return false;

	auto *range_spline = dynamic_cast<TSpline3*>(file->Get("re"));
	auto *energy_spline = dynamic_cast<TSpline3*>(file->Get("er"));
	if (!range_spline || !energy_spline) return false;

	range_energy_func_.reset(dynamic_cast<TSpline3*>(range_spline->Clone("re")));
	energy_range_func_.reset(dynamic_cast<TSpline3*>(energy_spline->Clone("er")));
	if (!range_energy_func_ || !energy_range_func_) return false;

	max_range_ = range_energy_func_->Eval(max_energy_);
	return true;
}

void RangeEnergyCalculator::BuildSplines() {
	std::vector<double> energy(kSplinePoints, 0.0);
	std::vector<double> range(kSplinePoints, 0.0);

	for (int i = 1; i < kSplinePoints; ++i) {
		double fraction = static_cast<double>(i) / double(kSplinePoints - 1);
		double current_energy = max_energy_ * fraction * fraction;
		energy[i] = current_energy;

		catima::Projectile projectile(double(mass_), charge_, 0, current_energy / double(mass_));
		double current_range = RangeGcm2ToUm(catima::range(projectile, material_), material_);
		if (current_range <= range[i - 1]) {
			current_range = range[i - 1] + kRangeEpsilon;
		}
		range[i] = current_range;
	}

	max_range_ = range.back();
	range_energy_func_ = std::make_unique<TSpline3>(
		"re",
		energy.data(),
		range.data(),
		kSplinePoints
	);
	energy_range_func_ = std::make_unique<TSpline3>(
		"er",
		range.data(),
		energy.data(),
		kSplinePoints
	);

	std::filesystem::path path(cache_path_);
	if (!path.parent_path().empty()) {
		std::filesystem::create_directories(path.parent_path());
	}
	TFile file(cache_path_.c_str(), "recreate");
	range_energy_func_->Write("re");
	energy_range_func_->Write("er");
}

double RangeEnergyCalculator::Range(double energy) const {
	if (energy <= 0.0) return 0.0;
	double clamped_energy = std::min(energy, max_energy_);
	return range_energy_func_->Eval(clamped_energy);
}

double RangeEnergyCalculator::Energy(double range) const {
	if (range <= 0.0) return 0.0;
	double clamped_range = std::min(range, max_range_);
	return energy_range_func_->Eval(clamped_range);
}

catima::Material SiliconMaterial() {
	catima::Material material;
	material.add_element(28.0, 14.0, 1.0);
	material.density(kSiliconDensity);
	return material;
}

} // namespace brill
