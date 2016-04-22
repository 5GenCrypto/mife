AUTOMAKE_OPTIONS = foreign
AM_CFLAGS = ${DEBUG_CFLAGS} -I$(top_srcdir) -I$(top_srcdir)/dgs -fopenmp -Wall -g

# AM_LDFLAGS = -Wl,-rpath -Wl,$(abs_top_builddir)/flint

LDADD = \
    $(top_builddir)/dgs/libdgs.la \
    $(top_builddir)/oz/liboz.la \
    $(top_builddir)/gghlite/libgghlite.la \
    $(top_builddir)/dgsl/libdgsl.la \
    -laesrand -lflint \
    -L$(top_builddir)/mife/jsmn -ljsmn \
    -L$(top_builddir)/mmap \
    -lgmp -lmpfr -lm -lssl -lcrypto -lmmap

LDFLAGS = -no-install

bin_PROGRAMS = gghlite_params \
               gghlite_instance \
               gghlite_noise_dist \
               nike \
               jigsaw \
               gghlite_false_friends \
               bench_dgsl \
               bench_prime_g \
               bench_ideal_span \
               bench_invert \
               bench_co_prime \
               bench_rem
