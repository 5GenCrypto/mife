AUTOMAKE_OPTIONS = foreign
AM_CFLAGS = ${DEBUG_CFLAGS} -I$(top_srcdir) -I$(top_srcdir)/dgs -fopenmp

LDADD = $(top_builddir)/gghlite/libgghlite.la \
	$(top_builddir)/dgs/libdgs.la \
	$(top_builddir)/oz/liboz.la \
        $(top_builddir)/dgsl/libdgsl.la \
        $(top_builddir)/flint/libflint.a -lgmp -lmpfr -lm 

LDFLAGS = -no-install

bin_PROGRAMS = gghlite_params \
               gghlite_instance \
               gghlite_noise_dist \
               nike \
               gghlite_false_friends \
               bench_dgsl \
               bench_prime_g \
               bench_ideal_span
