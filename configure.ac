AC_INIT([polycap],[0.1],[Pieter.Tack@UGent.be])
AC_PREREQ([2.60])
AC_CONFIG_SRCDIR([src/polycap.c])
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CANONICAL_HOST


LIB_CURRENT=0
LIB_REVISION=0
LIB_AGE=0
AC_SUBST(LIB_CURRENT)
AC_SUBST(LIB_REVISION)
AC_SUBST(LIB_AGE)
LIB_CURRENT_MINUS_AGE=`expr $LIB_CURRENT - $LIB_AGE`
AC_SUBST(LIB_CURRENT_MINUS_AGE)
AC_CONFIG_MACRO_DIR([m4])

AC_USE_SYSTEM_EXTENSIONS

m4_ifdef([AM_PROG_AR],[AM_PROG_AR])
#at least version 2.0 of libtool is required for creating the fortran bindings
LT_PREREQ([2.0.0])
#dlopen is necessary for IDL bindings
LT_INIT([dlopen disable-fast-install win32-dll])
#AC_LIBTOOL_DLOPEN
AC_PROG_LIBTOOL

AC_PROG_CC
AM_PROG_CC_C_O
#this next line may never be reached...
if test -z $CC ; then
	AC_MSG_ERROR([no C compiler was found on the system.])
fi

AC_CANONICAL_HOST

CFLAGS_BKP=$CFLAGS
CFLAGS=
AC_CHECK_FUNC([vasprintf], [], [AC_MSG_ERROR([vasprintf must be present on the system])])
AC_CHECK_FUNC([strdup], [], [AC_MSG_ERROR([strdup must be present on the system])])
CFLAGS=$CFLAGS_BKP

#initialize pkg-config
PKG_PROG_PKG_CONFIG

#search for xraylib and other modules
PKG_CHECK_MODULES([xraylib],libxrl >= 3.3.0)
#search for gsl
PKG_CHECK_MODULES([gsl],[gsl])
#try easyRNG first, gsl as fallback
PKG_CHECK_MODULES([easyRNG],[easyRNG >= 1.1],
	[
	AC_DEFINE([HAVE_EASYRNG], [], [building with easyRNG support])
	],[
	USELESS_ASSIGNMENT=1
	])

AC_LANG_PUSH(C)
AC_OPENMP
AC_LANG_POP(C)

if test $ac_cv_prog_c_openmp = unsupported ; then
	AC_MSG_ERROR([C compiler does not support OpenMP])
fi

#search for HDF5...
AC_PATH_PROGS([H5CC],[h5cc-64 h5cc],["fail"])
if test x$H5CC = xfail ; then
	AC_MSG_ERROR([An installation of HDF5 is required, including the headers])
fi

H5CC_CC=
H5CC_CFLAGS=
H5CC_LIBS=
H5CC_STATIC=

H5CC_PARSED=$($H5CC -show -shlib)

COUNTER="0"

AC_CHECK_PROGS([BASENAME], [gbasename basename],["fail"])
if test x"$BASENAME" = xfail ; then
	AC_MSG_ERROR([basename shell command not found])
fi

AC_CHECK_PROGS([DIRNAME], [gdirname dirname],["fail"])
if test x"$DIRNAME" = xfail ; then
	AC_MSG_ERROR([dirname shell command not found])
fi

for i in $H5CC_PARSED
do
	if test x$COUNTER = "x0" ; then
		H5CC_CC=$($BASENAME $i)
		COUNTER="1"
	elif test x${i:0:2} = "x-L" || test x${i:0:2} = "x-l" ; then
		H5CC_LIBS+="$i "
	elif test x${i:0:2} = "x-I" ; then
		H5CC_CFLAGS+="$i "
	elif test x${i:0:1} = "x/" && test x${i:${#i}-2:2} = "x.a" ; then
		H5CC_STATIC+="$i "
#	else
#		AC_MSG_ERROR([Unrecognized output from $H5FC: $i])
	fi
done

# h5cc sucks something awful and will not find the headers if they are not in a location that's in the preprocessor's search path
# use h5cc' location to add one more folder to the path
H5CC_CFLAGS+=" -I$($DIRNAME $($DIRNAME $H5CC))/include"

AC_SUBST(H5CC_LIBS)
AC_SUBST(H5CC_CFLAGS)
AC_SUBST(H5CC_STATIC)
AC_CONFIG_FILES([ 	Makefile
		 	include/Makefile
			src/Makefile
			tests/Makefile
			example/Makefile
		])
AC_CONFIG_HEADERS([config.h])

AC_OUTPUT
