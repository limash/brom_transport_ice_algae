#!/bin/bash

rm -r build/* && cp data/* build &&
cd build && cmake $BROMDIR/code -DFABM_BASE=$FABMDIR -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=~/local/brom
