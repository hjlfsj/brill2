./build/bin/match_dssd -r $1 t0d1 t0d2 t0d3 t0d4
./build/bin/estimate/estimate_normalize -r $1 t0d1 t0d2 t0d3 t0d4
./build/bin/track_t0 -r $1
./build/bin/rebuild_t0 -r $1
