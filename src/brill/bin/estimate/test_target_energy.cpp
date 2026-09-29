#include <iostream>
#include <string>

#include "external/cxxopts.hpp"
#include "include/Lise++/target_energy_loss.h"
#include "include/config.h"

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

int main(int argc, char **argv) {
	cxxopts::Options options(
		"test_target_energy",
		"Test target energy loss correction for half-target thickness."
	);
	options.add_options()
		("h,help", "Print help information.")
		("z,charge", "Particle atomic number Z.", cxxopts::value<int>(), "Z")
		("a,mass", "Particle mass number A.", cxxopts::value<int>(), "A")
		("E,energy", "Measured kinetic energy after half-target (MeV).", cxxopts::value<double>(), "MeV")
		("t,theta", "Lab angle in degrees (0-70).", cxxopts::value<double>(), "deg")
		("target", "Target type: cd2 or ch2.", cxxopts::value<std::string>(), "target")
		(
			"c,config",
			"Config file path.",
			cxxopts::value<std::string>()->default_value("config.toml"),
			"file"
		);

	auto result = options.parse(argc, argv);
	if (result.count("help")) {
		PrintUsage(options);
		return 0;
	}
	if (
		!result.count("charge") || !result.count("mass") ||
		!result.count("energy") || !result.count("theta") ||
		!result.count("target")
	) {
		std::cerr << "Error: Missing required options.\n";
		PrintUsage(options);
		return 1;
	}

	brill::AppConfig config;
	if (brill::LoadConfig(result["config"].as<std::string>(), config)) {
		return 1;
	}

	int Z = result["charge"].as<int>();
	int A = result["mass"].as<int>();
	double E = result["energy"].as<double>();
	double theta = result["theta"].as<double>();
	std::string target_str = result["target"].as<std::string>();

	if (theta < 0.0 || theta > 70.0) {
		std::cerr << "Error: theta must be in [0, 70] degrees.\n";
		return 1;
	}

	brill::TargetType target;
	if (target_str == "cd2") {
		target = brill::TargetType::kCD2;
	} else if (target_str == "ch2") {
		target = brill::TargetType::kCH2;
	} else {
		std::cerr << "Error: Unknown target '" << target_str
			<< "'. Use 'cd2' or 'ch2'.\n";
		return 1;
	}

	brill::TargetEnergyLoss calculator(config, Z, A);

	double half_thickness_um;
	double E0 = calculator.IncidentEnergy(E, theta, target, half_thickness_um);

	std::cout << "=== Target Energy Loss Correction ===\n";
	std::cout << "Particle: Z=" << Z << ", A=" << A << "\n";
	std::cout << "Target: " << target_str << "\n";
	std::cout << "Theta: " << theta << " deg\n";
	std::cout << "Measured energy E: " << E << " MeV\n";
	std::cout << "Effective half-thickness: " << half_thickness_um << " um\n";
	std::cout << "Incident energy E0: " << E0 << " MeV\n";

	return 0;
}