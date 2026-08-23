#!/bin/bash
# brill2 t0

./build/bin/ppac_normalize -r $1 -t main
./build/bin/sort_beam -r $1 -t main
./build/bin/sort_beam -r $1 -t t1
./build/bin/match_dssd -r $1 t0d1 t0d2 t0d3 t0d4 -t main
./build/bin/match_dssd -r $1 t0d1 t0d2 t0d3 t0d4 -t t1
./build/bin/track_ppac -r $1 -t main --draw
./build/bin/track_ppac -r $1 -t t1 --draw



