#!/bin/bash

mkdir -p build
rm -r build/*
cp data/* build &&
cd build && cmake $BROMDIR/src -DFABM_BASE=$FABMDIR -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=~/local/brom
