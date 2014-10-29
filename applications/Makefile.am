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
               bench_dgsl

gghlite_params_SOURCES = gghlite_params.c common.h
gghlite_instance_SOURCES = gghlite_instance.c common.h
gghlite_noise_dist_SOURCES = gghlite_noise_dist.c common.h
gghlite_false_friends_SOURCES = gghlite_false_friends.c common.h
nike_SOURCES = nike.c nike.h
bench_dgsl_SOURCES = bench_dgsl.c bench_dgsl.h
