AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4 --install
SUBDIRS = flint oz dgs dgsl

AM_CFLAGS  = ${DEBUG_CFLAGS} -I$(srcdir)/dgs/ -I$(srcdir)/flint -fopenmp

lib_LTLIBRARIES = libgghlite.la

MY_LIB_LDFLAGS=-release 0.0.$(RELEASE) -no-undefined

libgghlite_la_SOURCES = gghlite/gghlite.c \
                        gghlite/gghlite_pk.c \
                        gghlite/misc.c \
                        gghlite/lattice_reduction.c \
                        gghlite/gghlite.h \
                        gghlite/misc.h \
                        gghlite/ggh-defs.h \
                        gghlite/ggh-internals.h \
                        gghlite/api.c

libgghlite_la_LDFLAGS = $(MY_LIB_LDFLAGS)

include_HEADERS = gghlite/gghlite.h gghlite/gghlite-defs.h gghlite/gghlite-internals.h gghlite/misc.h

MY_PROGRAM_LDADD= $(srcdir)/libgghlite.la \
		  $(top_builddir)/dgs/libdgs.la \
		  $(top_builddir)/oz/liboz.la \
                  $(top_builddir)/dgsl/libdgsl.la \
                  $(top_builddir)/flint/libflint.a -lgmp -lmpfr -lm 

bin_PROGRAMS = gghlite_params gghlite_instance gghlite_noise_dist nike test_misc gghlite_false_friends bench_dgsl

gghlite_params_SOURCES = applications/gghlite_params.c applications/common.h
gghlite_params_LDADD = $(MY_PROGRAM_LDADD)
gghlite_params_LDFLAGS = -no-install

gghlite_instance_SOURCES = applications/gghlite_instance.c applications/common.h
gghlite_instance_LDADD = $(MY_PROGRAM_LDADD)
gghlite_instance_LDFLAGS = -no-install

gghlite_noise_dist_SOURCES = applications/gghlite_noise_dist.c applications/common.h
gghlite_noise_dist_LDADD = $(MY_PROGRAM_LDADD)
gghlite_noise_dist_LDFLAGS = -no-install

gghlite_false_friends_SOURCES = applications/gghlite_false_friends.c applications/common.h
gghlite_false_friends_LDADD = $(MY_PROGRAM_LDADD)
gghlite_false_friends_LDFLAGS = -no-install

nike_SOURCES = applications/nike.c applications/nike.h
nike_LDADD = $(MY_PROGRAM_LDADD)
nike_LDFLAGS = -no-install

bench_dgsl_SOURCES = applications/bench_dgsl.c applications/bench_dgsl.h
bench_dgsl_LDADD = $(MY_PROGRAM_LDADD)
bench_dgsl_LDFLAGS = -no-install


test_misc_SOURCES = tests/test_misc.c
test_misc_LDADD = $(MY_PROGRAM_LDADD)
test_misc_LDFLAGS = -no-install

check_PROGRAMS = test_dgsl test_resultant

test_dgsl_SOURCES = tests/test_dgsl.c
test_dgsl_LDADD = $(MY_PROGRAM_LDADD)
test_dgsl_LDFLAGS = -no-install

test_resultant_SOURCES = tests/test_resultant.c
test_resultant_LDADD = $(MY_PROGRAM_LDADD)
test_resultant_LDFLAGS = -no-install
