#pragma once

#include "include/config.h"
#include "include/event/ingot/ppac_event.h"
#include "include/event/ppac/ppac_track_event.h"

namespace brill {

struct PpacOffsetParams {
	double delay_x_peak = 0.0;
	double delay_x_sigma = 0.0;
	double delay_y_peak = 0.0;
	double delay_y_sigma = 0.0;
	double position_x_p0 = 0.0;
	double position_x_p1 = 0.0;
	double position_y_p0 = 0.0;
	double position_y_p1 = 0.0;
};

int ReadPpacOffsetParams(const std::string &path, PpacOffsetParams params[3]);

void TrackPpac(
	const PpacEvent &event,
	const PpacConfig &ppac_config,
	const PpacOffsetParams offset[3],
	PpacTrackEvent &track
);

} // namespace brill