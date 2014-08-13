AC_INIT([gghflint], [20140613])

AM_SILENT_RULES([yes])

AC_CONFIG_AUX_DIR([.])

AC_CONFIG_SRCDIR([src/gghlite.c])

AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([subdir-objects])

AC_CONFIG_HEADERS([src/config.h])

AC_PROG_LIBTOOL()

AC_PROG_CC_C99()

AS_IF([test "$ac_cv_prog_cc_c99" = "no"], AC_MSG_ERROR([C99 support is required but not found.]))

# Debugging support
AC_ARG_ENABLE(debug, [  --enable-debug          Enable assert() statements for debugging.], [enable_debug=yes])

if test "x$enable_debug" = x"yes"; then
   DEBUG_FLAGS="-g"
   AC_SUBST(DEBUG_FLAGS)
else
   AC_DEFINE(NDEBUG,1,[Define whether debugging is enabled])
fi

# Primality tests
AC_ARG_ENABLE(primality-test, [  --enable-primality-test    Check if g is prime (slow).], [enable_primality_test=yes])

if test "x$enable_primality_test" = x"yes"; then
   AC_DEFINE(GGHLITE_CHECK_PRIMALITY,1,[Define whether g is checked for primality])
fi

AC_CONFIG_FILES([Makefile])

AC_CONFIG_SUBDIRS([dgs])

# configure FLINT by hand because it doesn't support all options autotools might
# throw at it

cd flint && ./configure --prefix=$(PREFIX) && cd ..

AC_OUTPUT
