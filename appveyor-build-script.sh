#!/usr/bin/env bash

set -e
set -x

export PKG_CONFIG_PATH=$HOME/install/lib/pkgconfig
export PATH=$HOME/install/bin:$PATH
export CYTHON=cython
export HDF5_CFLAGS="-I/mingw64/include"
export HDF5_LIBS="-L/mingw64/lib -lhdf5"

cd $APPVEYOR_BUILD_FOLDER

autoreconf -fi
export PYTHON=python3
./configure --prefix=$HOME/install --enable-python
make
make check

mkdir build
cd build
meson -Dbuild-documentation=false -Dpython=python3 ..
ninja
ninja test
