AUTOMAKE_OPTIONS = foreign
SUBDIRS = ore_dir
AM_CFLAGS = ${DEBUG_CFLAGS} -I$(top_srcdir) -I$(top_srcdir)/dgs -fopenmp -Wall -g

AM_LDFLAGS = -Wl,-rpath -Wl,$(abs_top_builddir)/flint

LDADD = $(top_builddir)/aesrand/libaesrand.la \
	$(top_builddir)/mife/libmife.la \
	$(top_builddir)/dgs/libdgs.la \
	$(top_builddir)/oz/liboz.la \
	$(top_builddir)/gghlite/libgghlite.la \
        $(top_builddir)/dgsl/libdgsl.la \
	-L$(top_builddir)/flint -lflint \
	-L$(top_builddir)/mife/jsmn -ljsmn \
        -lgmp -lmpfr -lm -lssl -lcrypto -lclt13

LDFLAGS = -no-install

bin_PROGRAMS = gghlite_params \
               gghlite_instance \
               gghlite_noise_dist \
               nike \
               jigsaw \
							 ore \
               gghlite_false_friends \
               bench_dgsl \
               bench_prime_g \
               bench_ideal_span \
               bench_invert \
               bench_co_prime \
               bench_rem

ore_SOURCES = ore.h ore.c ore_dir/mife_glue.h ore_dir/mife_glue.c
