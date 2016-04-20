AUTOMAKE_OPTIONS = foreign
LDADD = \
	$(top_builddir)/mife/libmife.la \
	-L$(top_builddir)/mife/jsmn \
	$(top_builddir)/gghlite/libgghlite.la \
	-laesrand -lflint -ljsmn -lgmp -lmpfr -lm -lssl -lcrypto -lclt13 -lmmap
AM_LDFLAGS = -rpath $(abs_top_builddir)/flint
AM_CFLAGS = ${DEBUG_CFLAGS} -I$(top_srcdir) -I$(top_srcdir)/dgs -I$(top_srcdir)/flint
bin_PROGRAMS = keygen encrypt eval
keygen_SOURCES  =  keygen.c
encrypt_SOURCES = encrypt.c
eval_SOURCES    =    eval.c
