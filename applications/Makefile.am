AUTOMAKE_OPTIONS = foreign
#SUBDIRS = ore_dir
AM_CFLAGS = ${DEBUG_CFLAGS} -I$(top_srcdir) -I$(top_srcdir)/dgs -fopenmp

AM_LDFLAGS = -Wl,-rpath -Wl,$(abs_top_builddir)/flint

LDADD = $(top_builddir)/gghlite/libgghlite.la \
	$(top_builddir)/dgs/libdgs.la \
	$(top_builddir)/oz/liboz.la \
        $(top_builddir)/dgsl/libdgsl.la \
	-L$(top_builddir)/flint -lflint \
        -lgmp -lmpfr -lm

LDFLAGS = -no-install

include_HEADERS = ore.h

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
