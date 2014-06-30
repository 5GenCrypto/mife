AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4
SUBDIRS = flint dgs .

AM_CFLAGS  = ${DEBUG_CFLAGS} -I$(srcdir)/dgs/ -I$(srcdir)/flint -I$(srcdir)/flint/flint

bin_PROGRAMS = gghlite
gghlite_SOURCES = src/gghlite.c
gghlite_LDADD = $(srcdir)/dgs/libdgs.la $(srcdir)/flint/libflint.a -lgmp -lmpfr -lm
gghlite_LDFLAGS = -no-install
