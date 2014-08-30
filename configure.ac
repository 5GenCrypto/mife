AC_INIT([gghflint], [20140613])

AM_SILENT_RULES([yes])

AC_CONFIG_AUX_DIR([.])

AC_CONFIG_SRCDIR([gghlite/gghlite.c])

AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([subdir-objects])

AC_CONFIG_HEADERS([gghlite/config.h])

AC_PROG_LIBTOOL()

AC_PROG_CC_C99()

AS_IF([test "$ac_cv_prog_cc_c99" = "no"], AC_MSG_ERROR([C99 support is required but not found.]))

AC_ARG_ENABLE(debug,        [  --enable-debug          Enable assert() statements for debugging.], [enable_debug=yes])
AC_ARG_ENABLE(sloppy-checks,[  --enable-sloppy-checks  Disable expensive tests.], [enable_sloppy_checks=yes])
AC_ARG_ENABLE(heuristics,   [  --enable-heuristics     Enable various heuristic improvements.], [enable_heuristics=yes])

# Debugging support
if test "x$enable_debug" = x"yes"; then
   DEBUG_FLAGS="-g"
   AC_SUBST(DEBUG_FLAGS)
else
   AC_DEFINE(NDEBUG,1,[Define whether debugging is enabled])
fi

# Primality tests

if test "x$enable_sloppy_checks" = x"yes"; then
   AC_DEFINE(GGHLITE_SLOPPY,1,[Define whether expensive tests are disabled])
fi

# Heuristics
if test "x$enable_heuristics" = x"yes"; then
   AC_DEFINE(GGHLITE_HEURISTICS,1,[Define whether to enable performance improvements])
fi

AC_CONFIG_FILES([Makefile])

AC_CONFIG_SUBDIRS([dgs])

# configure FLINT by hand because it doesn't support all options autotools might
# throw at it

cd flint && ./configure --prefix=$(PREFIX) && cd ..

AC_OUTPUT
