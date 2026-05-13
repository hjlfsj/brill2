#include "include/analysis/tree.h"

#include <iostream>

#include <TFile.h>
#include <TTree.h>

namespace glimmer {

TTree *OpenInputTree(TFile &file, const std::string &path) {
	TTree *tree = dynamic_cast<TTree*>(file.Get(path.c_str()));
	if (!tree) {
		std::cerr << "Error: Get tree " << path << " failed.\n";
	}
	return tree;
}

} // namespace glimmer
