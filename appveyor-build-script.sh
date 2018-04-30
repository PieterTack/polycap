#!/usr/bin/env bash

set -e
set -x

export PKG_CONFIG_PATH=$HOME/install/lib/pkgconfig
export PATH=$HOME/install/bin:$PATH
export CYTHON=cython2

cd $APPVEYOR_BUILD_FOLDER

autoreconf -fi
export CFLAGS="-Wall -Werror"
./configure --prefix=$HOME/install --disable-python
make
make check
make distcheck
make distclean

export PYTHON=python2
./configure --prefix=$HOME/install --enable-python
make
make check
make distcheck
make distclean

export PYTHON=python3
./configure --prefix=$HOME/install --enable-python
make
make check
make distcheck
make distclean
