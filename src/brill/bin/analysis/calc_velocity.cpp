#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cctype>
#include <stdexcept>

#include "include/rebuild/constants.h"
#include "include/rebuild/nuclear_data.h"

constexpr double c = 299792458.0;

static void ParseNucleus(const std::string& input, int& a, std::string& element) {
	size_t i = 0;
	while (i < input.size() && std::isdigit(input[i])) {
		++i;
	}
	if (i == 0) {
		throw std::runtime_error("Invalid nucleus name: " + input);
	}
	a = std::stoi(input.substr(0, i));
	element = input.substr(i);
}

static brill::NuclearData FindNuclearByName(int a, const std::string& element) {
	std::ifstream fin("assets/amdc_ion_2020.txt");
	if (!fin.good()) {
		throw std::runtime_error("Cannot open assets/amdc_ion_2020.txt");
	}
	std::string line;
	std::getline(fin, line);
	brill::NuclearData data;
	while (fin >> data) {
		if (data.a == a && data.name == element) {
			return data;
		}
	}
	throw std::runtime_error("Nucleus not found: " + std::to_string(a) + element);
}

int main() {
	std::cout << "请输入核素 (e.g. 6Li): ";
	std::string nucleus;
	std::getline(std::cin, nucleus);

	std::cout << "请输入能量 (MeV/u): ";
	double energy_input;
	std::cin >> energy_input;

	int a;
	std::string element;
	ParseNucleus(nucleus, a, element);

	brill::NuclearData nd = FindNuclearByName(a, element);

	double mass_u = nd.mass;
	double mass_MeV = mass_u * brill::atomic_mass_unit;

	double E_kin = energy_input * a;
	double E_total = E_kin + mass_MeV;

	double gamma = E_total / mass_MeV;
	double beta = std::sqrt(1.0 - 1.0 / (gamma * gamma));
	double v = beta * c;
	double t_10m = 10.0 / v * 1e9;

	std::cout << "\n==============================" << std::endl;
	std::cout << "核素: " << nd.z << nd.name << " (Z=" << nd.z << ", A=" << nd.a << ")" << std::endl;
	std::cout << "质量: " << mass_u << " u = " << mass_MeV << " MeV/c^2" << std::endl;
	std::cout << "动能: " << E_kin << " MeV (" << energy_input << " MeV/u × " << a << ")" << std::endl;
	std::cout << "总能量: " << E_total << " MeV" << std::endl;
	std::cout << "洛伦兹因子 γ: " << gamma << std::endl;
	std::cout << "β: " << beta << std::endl;
	std::cout << "速度: " << v << " m/s" << std::endl;
	std::cout << "光速百分比: " << beta * 100.0 << "% c" << std::endl;
	std::cout << "飞行 10m 时间: " << t_10m << " ns" << std::endl;
	std::cout << "==============================" << std::endl;

	return 0;
}