#!/usr/bin/env bash

set -e
set -x

export PKG_CONFIG_PATH=$HOME/install/lib/pkgconfig
export PATH=$HOME/install/bin:$PATH

cd $APPVEYOR_BUILD_FOLDER

autoreconf -fi
export CFLAGS="-Wall -Werror"
./configure --prefix=$HOME/install
make
make check
make distcheck
