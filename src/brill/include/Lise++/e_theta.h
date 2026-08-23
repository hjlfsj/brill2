#ifndef __BRILL_LISE_ETHETA_H__
#define __BRILL_LISE_ETHETA_H__

#include <string>

class TGraph;

namespace brill {

TGraph* LoadEThetaCurve(const std::string& filename, const std::string& nucleus_name);
TGraph* LoadThetaThetaCurve(const std::string& filename);

} // namespace brill

#endif