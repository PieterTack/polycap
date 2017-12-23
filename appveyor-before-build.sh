#!/usr/bin/env bash

set -e
set -x

export PKG_CONFIG_PATH=$HOME/install/lib/pkgconfig
export PATH=$HOME/install/bin:$PATH

# install xraylib
wget -q https://xraylib.tomschoonjans.eu/xraylib-3.3.0.tar.gz
tar xfz xraylib-3.3.0.tar.gz
cd xraylib-3.3.0
./configure --prefix=$HOME/install --disable-static
make -j2
make install
cd ..

