#!/usr/bin/env bash

set -e
set -x

export PKG_CONFIG_PATH=$HOME/install/lib/pkgconfig
export PATH=$HOME/install/bin:$PATH

XRL_VERSION=4.0.0
EASYRNG_VERSION=1.2

# install xraylib
wget -q https://xraylib.tomschoonjans.eu/xraylib-${XRL_VERSION}.tar.gz
tar xfz xraylib-${XRL_VERSION}.tar.gz
cd xraylib-${XRL_VERSION}
./configure --prefix=$HOME/install --disable-static
make -j2
make install
cd ..

if test $RNG = "easyRNG" ; then
	wget -q https://github.com/tschoonj/easyRNG/releases/download/easyRNG-${EASYRNG_VERSION}/easyRNG-${EASYRNG_VERSION}.tar.gz
	tar xfz easyRNG-${EASYRNG_VERSION}.tar.gz
	cd easyRNG-${EASYRNG_VERSION}
	./configure --prefix=$HOME/install --disable-static
	make -j2
	make install
	cd ..
fi

