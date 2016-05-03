AUTOMAKE_OPTIONS = foreign -Wall
SUBDIRS = jsmn

MY_SOURCES = mife.c mife_io.c mife_internals.c flint_raw_io.c parse.c \
             mbp_glue.c util.c mbp_types.c jsmn/jsmn.c cmdline.c

AM_CFLAGS = ${DEBUG_CFLAGS} -I$(top_srcdir) -Ijsmn -std=c99 -D_DEFAULT_SOURCE -fopenmp

bin_PROGRAMS = keygen encrypt eval
keygen_SOURCES  =  keygen.c $(MY_SOURCES)
encrypt_SOURCES = encrypt.c $(MY_SOURCES)
eval_SOURCES    =    eval.c $(MY_SOURCES)
