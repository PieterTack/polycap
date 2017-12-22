#!/usr/bin/env bash

set -e
set -x

export PKG_CONFIG_PATH=$HOME/install/lib/pkgconfig
export PATH=$HOME/install/bin:$PATH

cd $APPVEYOR_BUILD_FOLDER

autoreconf -fi
export CPPFLAGS="-I/usr/local/include -I$HOME/install/include"
export CFLAGS="-Wno-deprecated -Wno-deprecated-declarations"
./configure --prefix=$HOME/install
make
make check
make distcheck
