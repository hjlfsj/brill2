#include "include/10C+4He/tafcsi.h"

namespace brill {

TafcsiTarget FindTafcsiTarget(const CsiEvent<kMaxTafcsi> &event, int tafd_index) {
	TafcsiTarget result;

	if (tafd_index < 0 || tafd_index >= 6) return result;

	int idx0 = 2 * tafd_index;
	int idx1 = 2 * tafd_index + 1;

	bool v0 = event.valid[idx0];
	bool v1 = event.valid[idx1];

	if (v0 && v1) {
		result.found = true;
		result.channel = (event.energy[idx0] >= event.energy[idx1]) ? idx0 : idx1;
		result.energy = event.energy[result.channel];
	} else if (v0) {
		result.found = true;
		result.channel = idx0;
		result.energy = event.energy[idx0];
	} else if (v1) {
		result.found = true;
		result.channel = idx1;
		result.energy = event.energy[idx1];
	}

	return result;
}

}  // namespace brill