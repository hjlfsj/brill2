#pragma once

#include <string>

class TFile;
class TTree;

namespace glimmer {

TTree *OpenInputTree(TFile &file, const std::string &path);

} // namespace glimmer
