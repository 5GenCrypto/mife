AUTOMAKE_OPTIONS = foreign
SUBDIRS = jsmn
LDADD = -L$(top_builddir)/applications/ore_dir/jsmn -ljsmn
AM_CFLAGS = ${DEBUG_CFLAGS} -I$(top_srcdir) -I$(top_srcdir)/dgs -I$(top_srcdir)/applications/ore_dir/jsmn
bin_PROGRAMS = keygen
keygen_SOURCES = keygen.c matrix.h matrix.c parse.h parse.c
