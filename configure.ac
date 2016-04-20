AC_INIT([gghflint], [20140613])

AC_CONFIG_AUX_DIR([.])

AC_CONFIG_SRCDIR([gghlite])

AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([subdir-objects])

# support older automake versions
m4_ifdef([AM_SILENT_RULES],
         [AM_SILENT_RULES([yes])],
         [AC_SUBST([AM_DEFAULT_VERBOSITY], [1])])

AC_CONFIG_HEADERS([gghlite/config.h])

m4_pattern_allow([AM_PROG_AR])
AM_PROG_AR()

AC_PROG_LIBTOOL()

AC_PROG_CC_C99()

AS_IF([test "$ac_cv_prog_cc_c99" = "no"], AC_MSG_ERROR([C99 support is required but not found.]))

AC_ARG_ENABLE(debug,        [  --enable-debug          Enable assert() statements for debugging.], [enable_debug=yes])
AC_ARG_ENABLE(heuristics,   [  --enable-heuristics     Enable various heuristic improvements.], [enable_heuristics=yes])

# Debugging support
if test "x$enable_debug" = x"yes"; then
   DEBUG_FLAGS="-g"
   AC_SUBST(DEBUG_FLAGS)
else
   AC_DEFINE(NDEBUG,1,[Define whether debugging is enabled])
fi

AX_VALGRIND_CHECK()

# Heuristics
if test "x$enable_heuristics" = x"yes"; then
   AC_DEFINE(GGHLITE_HEURISTICS,1,[Define whether to enable performance improvements])
fi

AC_SEARCH_LIBS(aes_randinit,aesrand)
if test "x$ac_cv_search_aes_randinit" = "xno"; then
   AC_MSG_ERROR([libaesrand not found])
fi
AC_SEARCH_LIBS(clt_state_init,clt13)
if test "x$ac_cv_search_clt_state_init" = "xno"; then
   AC_MSG_ERROR([libclt13 not found])
fi
AC_SEARCH_LIBS(fmpz_init,flint)
if test "x$ac_cv_search_fmpz_init" = "xno"; then
   AC_MSG_ERROR([libflint not found])
fi

AC_CONFIG_FILES([Makefile oz/Makefile dgsl/Makefile gghlite/Makefile applications/Makefile tests/Makefile applications/ore_dir/Makefile mife/Makefile mmap/Makefile])

AC_CONFIG_SUBDIRS([dgs])

# configure FLINT by hand because it doesn't support all options autotools might
# throw at it

dnl AC_ARG_VAR(FLINT, [If set to no then calling configure on FLINT is skipped])

dnl if test "x$FLINT" != "xno"; then
dnl    cd flint && ./configure --prefix=${prefix} && cd ..
dnl fi   

AC_OUTPUT
