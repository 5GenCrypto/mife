m4_define([version_major], [0])
m4_define([version_minor], [9])
m4_define([version_micro], [0])

AC_INIT([libmife], [version_major.version_minor.version_micro])

AC_CONFIG_HEADERS([mife/config.h])
AC_CONFIG_SRCDIR([./])

AC_CONFIG_AUX_DIR([build/autoconf])
AC_CONFIG_MACRO_DIR([build/autoconf])
AM_INIT_AUTOMAKE([subdir-objects])

AM_PROG_AR
dnl AC_PROG_LIBTOOL
dnl AC_PROG_CC_C99

dnl AS_IF([test "$ac_cv_prog_cc_c99" = "no"], AC_MSG_ERROR([C99 support is required but not found.]))

dnl AC_ARG_ENABLE(debug,        [  --enable-debug          Enable assert() statements for debugging.], [enable_debug=yes])

dnl # Debugging support
dnl if test "x$enable_debug" = x"yes"; then
dnl    DEBUG_FLAGS="-g"
dnl    AC_SUBST(DEBUG_FLAGS)
dnl else
dnl    AC_DEFINE(NDEBUG,1,[Define whether debugging is enabled])
dnl fi

dnl AX_VALGRIND_CHECK()

AC_CONFIG_FILES([Makefile mife/Makefile])

AC_OUTPUT
