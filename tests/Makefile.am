AUTOMAKE_OPTIONS = foreign
AM_CFLAGS  = ${DEBUG_CFLAGS} -I$(top_srcdir) -I$(top_srcdir)/dgs -fopenmp

MY_PROGRAM_LDADD= $(top_builddir)/libgghlite.la \
		  $(top_builddir)/dgs/libdgs.la \
		  $(top_builddir)/oz/liboz.la \
                  $(top_builddir)/dgsl/libdgsl.la \
                  $(top_builddir)/flint/libflint.a -lgmp -lmpfr -lm 

test_dgsl_SOURCES = test_dgsl.c
test_dgsl_LDADD = $(MY_PROGRAM_LDADD)
test_dgsl_LDFLAGS = -no-install

test_norm_SOURCES = test_norm.c
test_norm_LDADD = $(MY_PROGRAM_LDADD)
test_norm_LDFLAGS = -no-install

test_invert_SOURCES = test_invert.c
test_invert_LDADD = $(MY_PROGRAM_LDADD)
test_invert_LDFLAGS = -no-install

TESTS = test_dgsl test_norm test_invert
check_PROGRAMS = test_dgsl test_norm test_invert

