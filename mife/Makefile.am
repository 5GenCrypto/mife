AUTOMAKE_OPTIONS = foreign -Wall
SUBDIRS = jsmn

MY_SOURCES = mife.c mife_io.c mife_internals.c flint_raw_io.c parse.c \
             mbp_glue.c util.c mbp_types.c jsmn/jsmn.c cmdline.c

AM_CFLAGS = ${DEBUG_CFLAGS} -I$(top_srcdir) -Ijsmn -fopenmp

LDADD = -lmmap -lgghlite -lflint -laesrand -loz -lgmp

# pkgincludesubdir = $(includedir)/src
# pkgincludesub_HEADERS = mife.h mife_io.h mife_internals.h flint_raw_io.h \
# 	mife_defs.h parse.h mbp_glue.h util.h cmdline.h mbp_types.h

# AM_LDFLAGS = -rpath $(abs_top_builddir)/flint

bin_PROGRAMS = keygen encrypt eval
keygen_SOURCES  =  keygen.c $(MY_SOURCES)
encrypt_SOURCES = encrypt.c $(MY_SOURCES)
eval_SOURCES    =    eval.c $(MY_SOURCES)
