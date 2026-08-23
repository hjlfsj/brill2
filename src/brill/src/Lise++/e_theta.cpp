#include "include/Lise++/e_theta.h"

#include <TAxis.h>
#include <TGraph.h>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <cctype>

namespace brill {

static std::vector<std::string> SplitByTab(const std::string& line) {
	std::vector<std::string> tokens;
	std::string token;
	for (char c : line) {
		if (c == '\t') {
			tokens.push_back(token);
			token.clear();
		} else {
			token += c;
		}
	}
	if (!token.empty()) {
		tokens.push_back(token);
	}
	return tokens;
}

static int ParseA(const std::string& nucleus_name) {
	std::string a_str;
	for (char c : nucleus_name) {
		if (std::isdigit(c)) {
			a_str += c;
		} else {
			break;
		}
	}
	if (a_str.empty()) {
		throw std::runtime_error("Cannot parse mass number from nucleus name: " + nucleus_name);
	}
	return std::stoi(a_str);
}

TGraph* LoadEThetaCurve(const std::string& filename, const std::string& nucleus_name) {
	int a = ParseA(nucleus_name);

	std::ifstream fin(filename);
	if (!fin.good()) {
		throw std::runtime_error("Cannot open file: " + filename);
	}

	std::string line;
	while (std::getline(fin, line)) {
		if (line.empty()) continue;
		if (line[0] == '!') continue;
		break;
	}

	std::vector<std::string> headers = SplitByTab(line);

	int col_angle = -1;
	int col_energy = -1;
	for (size_t i = 0; i < headers.size(); ++i) {
		if (headers[i].find(nucleus_name) != std::string::npos) {
			col_energy = (int)i;
			col_angle = (int)i - 1;
			break;
		}
	}
	if (col_energy < 0) {
		throw std::runtime_error(
			"Nucleus '" + nucleus_name + "' not found in header of " + filename);
	}

	std::vector<double> angles;
	std::vector<double> energies;
	while (std::getline(fin, line)) {
		if (line.empty()) continue;
		std::vector<std::string> tokens = SplitByTab(line);
		if ((int)tokens.size() <= col_energy) continue;

		double angle = std::stod(tokens[col_angle]);
		double e_per_u = std::stod(tokens[col_energy]);
		double e_total = e_per_u * a;

		angles.push_back(angle);
		energies.push_back(e_total);
	}

	if (angles.empty()) {
		throw std::runtime_error("No data points found in " + filename);
	}

	TGraph* graph = new TGraph((int)angles.size(), angles.data(), energies.data());
	graph->SetLineColor(kRed);
	graph->SetLineWidth(2);
	graph->SetTitle((nucleus_name + " E-#theta").c_str());
	graph->GetXaxis()->SetTitle("Angle (Lab-deg)");
	graph->GetYaxis()->SetTitle("Total Energy (MeV)");

	return graph;
}

TGraph* LoadThetaThetaCurve(const std::string& filename) {
	std::ifstream fin(filename);
	if (!fin.good()) {
		throw std::runtime_error("Cannot open file: " + filename);
	}

	std::string line;
	while (std::getline(fin, line)) {
		if (line.empty()) continue;
		if (line[0] == '!') continue;
		break;
	}

	std::vector<double> angles_x;
	std::vector<double> angles_y;
	while (std::getline(fin, line)) {
		if (line.empty()) continue;
		std::vector<std::string> tokens = SplitByTab(line);
		if (tokens.size() < 2) continue;
		angles_x.push_back(std::stod(tokens[0]));
		angles_y.push_back(std::stod(tokens[1]));
	}

	if (angles_x.empty()) {
		throw std::runtime_error("No data points found in " + filename);
	}

	TGraph* graph = new TGraph((int)angles_x.size(), angles_x.data(), angles_y.data());
	graph->SetLineColor(kBlue);
	graph->SetLineWidth(2);
	graph->SetTitle("#theta-#theta correlation");
	graph->GetXaxis()->SetTitle("Angle (Lab-deg)");
	graph->GetYaxis()->SetTitle("Angle (Lab-deg)");

	return graph;
}

} // namespace brill