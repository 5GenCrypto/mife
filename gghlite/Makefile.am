AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = libgghlite.la

include_HEADERS = gghlite.h \
                  gghlite-defs.h \
                  gghlite-internals.h \
                  misc.h

libgghlite_la_SOURCES = gghlite.c \
                        gghlite_pk.c \
                        misc.c \
                        lattice_reduction.c \
                        gghlite.h \
                        misc.h \
                        ggh-defs.h \
                        ggh-internals.h \
                        api.c
libgghlite_la_LDFLAGS = -release 0.0.$(RELEASE) -no-undefined
libgghlite_la_INCLUDEDIR = $(includedir)/gghlite
