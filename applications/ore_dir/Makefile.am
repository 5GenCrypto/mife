AUTOMAKE_OPTIONS = foreign
SUBDIRS = jsmn
LDADD = \
	$(top_builddir)/mife/libmife.la \
	-L$(top_builddir)/applications/ore_dir/jsmn \
	-L$(top_builddir)/flint -lflint -ljsmn -lgmp -lmpfr -lm -lssl -lcrypto
AM_CFLAGS = ${DEBUG_CFLAGS} -I$(top_srcdir) -I$(top_srcdir)/dgs -I$(top_srcdir)/applications/ore_dir/jsmn -I$(top_srcdir)/mife -I$(top_srcdir)/flint
bin_PROGRAMS = keygen
keygen_SOURCES = keygen.c parse.h parse.c mife_glue.h mife_glue.c
