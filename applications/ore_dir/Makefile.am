AUTOMAKE_OPTIONS = foreign
LDADD = \
	$(top_builddir)/aesrand/libaesrand.la \
	$(top_builddir)/mife/libmife.la \
	-L$(top_builddir)/mife/jsmn \
	-L$(top_builddir)/flint -lflint -ljsmn -lgmp -lmpfr -lm -lssl -lcrypto
AM_LDFLAGS = -rpath $(abs_top_builddir)/flint
AM_CFLAGS = ${DEBUG_CFLAGS} -I$(top_srcdir) -I$(top_srcdir)/dgs -I$(top_srcdir)/applications/ore_dir/jsmn -I$(top_srcdir)/mife -I$(top_srcdir)/flint
COMMON_SOURCES = mife_glue.h mife_glue.c
bin_PROGRAMS = keygen encrypt eval
keygen_SOURCES  =  keygen.c $(COMMON_SOURCES)
encrypt_SOURCES = encrypt.c $(COMMON_SOURCES)
eval_SOURCES    =    eval.c $(COMMON_SOURCES)
