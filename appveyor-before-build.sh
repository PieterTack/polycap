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

# install hdf5
curl -L -s -O https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.12/src/hdf5-1.8.12.tar.gz
tar xfz hdf5-1.8.12.tar.gz
cd hdf5-1.8.12
# add support for UTF-8 filenames
curl -L -s -O https://www.dropbox.com/s/gowzeo6vdhjpxnw/hdf5-1.8.12.diff
patch -p1 < hdf5-1.8.12.diff
autoreconf -i
./configure --disable-fortran --disable-cxx --disable-hl --prefix=$HOME/install --disable-static CPPFLAGS=-D_GNU_SOURCE=1
# patch hdf5 -> https://tschoonj.github.io/blog/2014/01/29/building-a-64-bit-version-of-hdf5-with-mingw-w64/
echo "#ifndef H5_HAVE_WIN32_API" >> src/H5pubconf.h
echo "#ifdef WIN32 /* defined for all windows systems */" >> src/H5pubconf.h
echo "#define H5_HAVE_WIN32_API 1" >> src/H5pubconf.h
echo "#endif" >> src/H5pubconf.h
echo "#endif" >> src/H5pubconf.h
echo "#ifndef H5_HAVE_MINGW" >> src/H5pubconf.h
echo "#ifdef __MINGW32__ /*defined for all MinGW compilers */" >> src/H5pubconf.h
echo "#define H5_HAVE_MINGW 1" >> src/H5pubconf.h
echo "#define H5_HAVE_WINDOWS 1" >> src/H5pubconf.h
echo "#endif" >> src/H5pubconf.h
echo "#endif" >> src/H5pubconf.h
echo "#define H5_BUILT_AS_DYNAMIC_LIB 1" >> src/H5pubconf.h
make -j2
make install
cd ..
