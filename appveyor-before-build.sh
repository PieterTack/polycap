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

if test $RNG = "easyRNG" ; then
	wget -q https://github.com/tschoonj/easyRNG/releases/download/easyRNG-1.1/easyRNG-1.1.tar.gz
	tar xfz easyRNG-1.1.tar.gz
	cd easyRNG-1.1
	./configure --prefix=$HOME/install --disable-static
	make -j2
	make install
	cd ..
fi

