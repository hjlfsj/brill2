#pragma once

#include "include/event/ingot/csi_event.h"

namespace brill {

constexpr int kMaxTafcsi = 12;

struct TafcsiTarget {
	bool found = false;
	int energy = 0;
	int channel = -1;
};

TafcsiTarget FindTafcsiTarget(const CsiEvent<kMaxTafcsi> &event, int tafd_index);

}  // namespace brill