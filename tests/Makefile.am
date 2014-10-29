AUTOMAKE_OPTIONS = foreign
AM_CFLAGS  = ${DEBUG_CFLAGS} -I$(top_srcdir) -I$(top_srcdir)/dgs -fopenmp

LDADD = $(top_builddir)/gghlite/libgghlite.la \
        $(top_builddir)/dgs/libdgs.la \
        $(top_builddir)/oz/liboz.la \
        $(top_builddir)/dgsl/libdgsl.la \
        $(top_builddir)/flint/libflint.a -lgmp -lmpfr -lm 

LDFLAGS = -no-install

TESTS = test_dgsl test_norm test_invert test_mul
check_PROGRAMS = $(TESTS)
