AC_INIT([polycap],[0.1],[Pieter.Tack@UGent.be])
AC_PREREQ([2.60])
AC_CONFIG_SRCDIR([src/polycap-description.c])
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
LT_INIT([dlopen disable-fast-install win32-dll disable-static])
#AC_LIBTOOL_DLOPEN
AC_PROG_LIBTOOL

AC_PROG_CC
AM_PROG_CC_C_O
#this next line may never be reached...
if test -z $CC ; then
	AC_MSG_ERROR([no C compiler was found on the system.])
fi

AC_CANONICAL_HOST

OS_WINDOWS=0
case "$host" in
	*mingw*)
		#build dll
		OS_WINDOWS=1
		LDFLAGS_LIBPOLYCAP="-no-undefined"
		#LDFLAGS_FORTRAN="-no-undefined"
		LDFLAGS_PYTHON="-no-undefined -Wl,-subsystem,windows"
		;;
esac
AM_CONDITIONAL([OS_WINDOWS],[test x$OS_WINDOWS = x1])
AC_SUBST(OS_WINDOWS)

#AC_SUBST(LDFLAGS_FORTRAN)
AC_SUBST(LDFLAGS_LIBPOLYCAP)
AC_SUBST(LDFLAGS_PYTHON)

CFLAGS_BKP=$CFLAGS
CFLAGS=
AC_CHECK_FUNC([vasprintf], [], [AC_MSG_ERROR([vasprintf must be present on the system])])
AC_CHECK_FUNC([strdup], [], [AC_MSG_ERROR([strdup must be present on the system])])
AC_CHECK_FUNC([gettimeofday], [], [AC_MSG_ERROR([strdup must be present on the system])])
CFLAGS=$CFLAGS_BKP

# Symbol visibility handling.
#
# Taken from gtksourceview and modified where necessary
HIDDEN_VISIBILITY_CFLAGS=""
case "$host" in
    *-*-mingw*)
        dnl on mingw32 we do -fvisibility=hidden and __declspec(dllexport)
        AC_DEFINE([POLYCAP_EXTERN], [__attribute__((visibility("default"))) __declspec(dllexport) extern],
                  [defines how to decorate public symbols while building])
        HIDDEN_VISIBILITY_CFLAGS="-fvisibility=hidden"
        ;;
    *)
        dnl on other compilers, check if we can do -fvisibility=hidden
        SAVED_CFLAGS="${CFLAGS}"
        CFLAGS="-fvisibility=hidden"
        AC_MSG_CHECKING([for -fvisibility=hidden compiler flag])
        AC_TRY_COMPILE([], [return 0],
                       AC_MSG_RESULT(yes)
                       enable_fvisibility_hidden=yes,
                       AC_MSG_RESULT(no)
                       enable_fvisibility_hidden=no)
        CFLAGS="${SAVED_CFLAGS}"

        AS_IF([test "${enable_fvisibility_hidden}" = "yes"], [
            AC_DEFINE([POLYCAP_EXTERN], [__attribute__((visibility("default"))) extern],
                      [defines how to decorate public symbols while building])
            HIDDEN_VISIBILITY_CFLAGS="-fvisibility=hidden"
        ])
        ;;
esac
AC_SUBST(HIDDEN_VISIBILITY_CFLAGS)



#initialize pkg-config
PKG_PROG_PKG_CONFIG

#search for xraylib and other modules
PKG_CHECK_MODULES([xraylib],libxrl >= 4.0.0)
#search for gsl
PKG_CHECK_MODULES([gsl],[gsl])
#try easyRNG first, gsl as fallback
PKG_CHECK_MODULES([easyRNG],[easyRNG >= 1.1],
	[
	AC_DEFINE([HAVE_EASYRNG], [], [building with easyRNG support])
	AC_SUBST([EASYRNG_PC], "easyRNG")
	],[
	:
	])

AC_LANG_PUSH(C)
AC_OPENMP
AC_LANG_POP(C)

if test $ac_cv_prog_c_openmp = unsupported ; then
	AC_MSG_ERROR([C compiler does not support OpenMP])
fi


AC_ARG_VAR(HDF5_CFLAGS, [C compiler flags for HDF5])
AC_ARG_VAR(HDF5_LIBS, [linker flags for HDF5])
 
if test "x$HDF5_CFLAGS" = "x" && test "x$HDF5_LIBS" = "x" ; then

#search for HDF5...
AC_PATH_PROGS([H5CC],[h5cc-64 h5cc],["fail"])
if test x$H5CC = xfail ; then
	AC_MSG_ERROR([An installation of HDF5 is required, including the headers])
fi

H5CC_CC=
HDF5_CFLAGS=
HDF5_LIBS=

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
		HDF5_LIBS+="$i "
	elif test x${i:0:2} = "x-I" ; then
		HDF5_CFLAGS+="$i "
#	else
#		AC_MSG_ERROR([Unrecognized output from $H5FC: $i])
	fi
done

# h5cc sucks something awful and will not find the headers if they are not in a location that's in the preprocessor's search path
# use h5cc' location to add one more folder to the path
HDF5_CFLAGS+=" -I$($DIRNAME $($DIRNAME $H5CC))/include"

AC_SUBST(HDF5_LIBS)
AC_SUBST(HDF5_CFLAGS)

fi

#
#Python bindings using cython
#
AC_PROG_SED
AC_ARG_ENABLE([python],[AS_HELP_STRING([--disable-python],[build without the pythonbindings])],[enable_python=$enableval],[enable_python=check])
#for this we need to check for cython and for numpy, and python itself obviously
VALID_PYTHON=
AC_ARG_VAR(CYTHON, [the Cython command])
if test "x$enable_python" != xno  ; then
	##check for python
	AM_PATH_PYTHON([3.5], [PYTHON_FOUND=true], [PYTHON_FOUND=false])
	if test "x$PYTHON_FOUND" = xtrue ; then
		PYTHON_CPPFLAGS=
		PYTHON_LDFLAGS=
		AX_PYTHON_DEVEL
		if test "x$PYTHON" = x ; then
			if test "x$enable_python" = xyes ; then
				AC_MSG_ERROR([Incomplete python development package])
			else
				AC_MSG_WARN([Incomplete python development package])
			fi
			VALID_PYTHON_EXE=no
		else
			# check for numpy!
			AC_MSG_CHECKING([for numpy])
			res=`$PYTHON -c "import numpy" &>/dev/null`
			NUMPY=$?
			if test x$NUMPY = x0 ; then
				AC_MSG_RESULT([yes])
				NUMPY_HEADERS_PATH=`$PYTHON -c "from __future__ import print_function; import numpy; \
					print(numpy.get_include())"`
				if test $OS_WINDOWS = 1 ; then
					NUMPY_HEADERS_PATH=`cygpath -u $NUMPY_HEADERS_PATH`
				fi
				NUMPY_HEADERS=-I$NUMPY_HEADERS_PATH
				AC_SUBST(NUMPY_HEADERS)
				VALID_PYTHON_EXE=yes
			else
				VALID_PYTHON_EXE=no
				AC_MSG_RESULT([no])
			fi
		fi
	fi
	if test "x$VALID_PYTHON_EXE" = xyes ; then
		#now check cython
		AC_CHECK_PROGS(CYTHON, [cython-[$PYTHON_VERSION] cython[$PYTHON_VERSION] cython3 cython],["nocython"])
		if test x$CYTHON != "xnocython" ; then
			VALID_PYTHON=yes
			AC_MSG_NOTICE([Building with Python bindings])
		elif test "x$enable_python" = "xyes" ; then
			AC_MSG_ERROR([Cannot build Python bindings])
		else
			AC_MSG_WARN([Cannot build Python-NumPy bindings])
			VALID_PYTHON=no
		fi
		# check for OpenMP
		AC_OPENMP
	fi
fi
AM_CONDITIONAL([ENABLE_PYTHON],[test x$VALID_PYTHON = xyes])

if test "x$VALID_PYTHON" = xyes ; then
AC_SUBST(PYTHON_POLYCAP_VERSION,$VERSION)
AC_SUBST(PYTHONDIR,$pythondir)
AC_SUBST(PKGPYTHONDIR,$pkgpythondir)
AC_SUBST(PYEXECDIR,$pyexecdir)
AC_SUBST(PKGPYEXECDIR,$pkgpyexecdir)
fi

AC_ARG_ENABLE([python-integration],[AS_HELP_STRING([--enable-python-integration],[install the python bindings in the interpreters site-packages folder])],[enable_python_integration=$enableval],[enable_python_integration=check])
if test "x$enable_python_integration" = xyes ; then
	pythondir=$PYTHON_SITE_PKG
	pyexecdir=$PYTHON_SITE_PKG_EXEC
fi

# documentation
AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test -z "$DOXYGEN";
   then AC_MSG_WARN([Doxygen not found - continuing without Doxygen support])
fi
AM_CONDITIONAL([HAVE_DOXYGEN], [test -n "$DOXYGEN"])
AM_COND_IF([HAVE_DOXYGEN], [AC_CONFIG_FILES([docs/Doxyfile])])

AC_CONFIG_FILES([ 	Makefile
		 	include/Makefile
			src/Makefile
			python/Makefile
			tests/Makefile
			example/Makefile
			docs/Makefile
			libpolycap.pc
			polycap.spec
		])
AC_CONFIG_HEADERS([config.h])

AC_OUTPUT
