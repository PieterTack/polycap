AM_CPPFLAGS = -I${top_srcdir}/src -I$(top_srcdir)/include

check_PROGRAMS = test1
TESTS = $(check_PROGRAMS)

test1_SOURCES = test1.c
test1_LDADD = ../src/libpolycap.la
test1_CFLAGS = @OPENMP_CFLAGS@
test1_LDFLAGS = @OPENMP_CFLAGS@
