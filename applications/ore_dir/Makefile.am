AUTOMAKE_OPTIONS = foreign
SUBDIRS = jsmn
LDADD = \
	-L$(top_builddir)/mife \
	-L$(top_builddir)/applications/ore_dir/jsmn \
	-L$(top_builddir)/flint \
	-ljsmn -lmife
AM_CFLAGS = ${DEBUG_CFLAGS} -I$(top_srcdir) -I$(top_srcdir)/dgs -I$(top_srcdir)/applications/ore_dir/jsmn -I$(top_srcdir)/mife
bin_PROGRAMS = keygen
keygen_SOURCES = keygen.c parse.h parse.c
