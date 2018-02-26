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

AC_CONFIG_FILES([ 	Makefile
		 	include/Makefile
			src/Makefile
			tests/Makefile
		])
AC_CONFIG_HEADERS([config.h])

AC_OUTPUT
