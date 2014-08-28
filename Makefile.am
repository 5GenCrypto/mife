AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4
SUBDIRS = flint dgs

AM_CFLAGS  = ${DEBUG_CFLAGS} -I$(srcdir)/dgs/ -I$(srcdir)/flint #-I$(srcdir)/flint/flint

lib_LTLIBRARIES = libgpv.la libgghlite.la

MY_LIB_LDFLAGS=-release 0.0.$(RELEASE) -no-undefined

libgpv_la_SOURCES = gpv/gpv.c gpv/gso.c
libgpv_la_LDFLAGS = $(MY_LIB_LDFLAGS)

libgghlite_la_SOURCES = gghlite/gghlite.c \
                        gghlite/gghlite_pk.c \
                        gghlite/misc.c \
                        gghlite/gghlite.h \
                        gghlite/misc.h \
                        gghlite/ggh-defs.h \
                        gghlite/ggh-internals.h \
                        gghlite/api.c

libgghlite_la_LDFLAGS = $(MY_LIB_LDFLAGS)

include_HEADERS = gpv/gpv.h gpv/gso.h gghlite/gghlite.h gghlite/misc.h

MY_PROGRAM_LDADD=$(srcdir)/libgghlite.la $(srcdir)/dgs/libdgs.la $(srcdir)/libgpv.la $(srcdir)/flint/libflint.a -lgmp -lmpfr -lm

bin_PROGRAMS = gghlite_params gghlite_instance gghlite_noise_dist nike
gghlite_params_SOURCES = applications/gghlite_params.c applications/common.h
gghlite_params_LDADD = $(MY_PROGRAM_LDADD)
gghlite_params_LDFLAGS = -no-install

gghlite_instance_SOURCES = applications/gghlite_instance.c applications/common.h
gghlite_instance_LDADD = $(MY_PROGRAM_LDADD)
gghlite_instance_LDFLAGS = -no-install

gghlite_noise_dist_SOURCES = applications/gghlite_noise_dist.c applications/common.h
gghlite_noise_dist_LDADD = $(MY_PROGRAM_LDADD)
gghlite_noise_dist_LDFLAGS = -no-install

nike_SOURCES = applications/nike.c applications/nike.h
nike_LDADD = $(MY_PROGRAM_LDADD)
nike_LDFLAGS = -no-install

check_PROGRAMS = test_gpv
test_gpv_SOURCES = tests/test_gpv.c
test_gpv_LDADD = $(srcdir)/dgs/libdgs.la $(srcdir)/libgpv.la $(srcdir)/flint/libflint.a -lgmp -lmpfr -lm
test_gpv_LDFLAGS = -no-install
