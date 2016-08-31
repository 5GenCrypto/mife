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

AC_ARG_ENABLE(debug,        [  --enable-debug          Enable assert() statements for debugging.], [enable_debug=yes])

CFLAGS=                         dnl get rid of default -g -O2
COMMON_CFLAGS="-Wall -Wformat -Wformat-security -Wextra -Wunused \
-Wshadow -Wmissing-prototypes -Wfloat-equal -Wpointer-arith -Wcast-align \
-Wstrict-prototypes -Wredundant-decls -Wendif-labels -Wcast-qual \
-std=gnu11 -Wpedantic"

if test "x$enable_debug" = x"yes"; then
  EXTRA_CFLAGS="-O0 -g"
else
  EXTRA_CFLAGS="-O3"
  AC_DEFINE(NDEBUG,1,[Define whether debugging is enabled])
fi
AC_SUBST(COMMON_CFLAGS)
AC_SUBST(EXTRA_CFLAGS)

AC_CHECK_HEADERS([omp.h])

AC_SEARCH_LIBS(aes_randinit,aesrand)
if test "x$ac_cv_search_aes_randinit" = "xno"; then
  AC_MSG_ERROR([libaesrand not found])
fi
AC_SEARCH_LIBS(clt_state_new,clt13)
if test "x$ac_cv_search_clt_state_new" = "xno"; then
AC_MSG_ERROR([libclt13 not found])
fi
AC_SEARCH_LIBS(fmpz_init,flint)
if test "x$ac_cv_search_fmpz_init" = "xno"; then
  AC_MSG_ERROR([libflint not found])
fi
AC_SEARCH_LIBS(fmpz_mod_poly_oz_intt,oz)
if test "x$ac_cv_search_fmpz_mod_poly_oz_intt" = "xno"; then
  AC_MSG_ERROR([liboz not found])
fi
AC_SEARCH_LIBS(gghlite_init,gghlite)
if test "x$ac_cv_search_gghlite_init" = "xno"; then
  AC_MSG_ERROR([libgghlite not found])
fi
AC_SEARCH_LIBS(__gmpz_init,gmp)
if test "x$ac_cv_search___gmpz_init" = "xno"; then
  AC_MSG_ERROR([libgmp not found])
fi
AC_SEARCH_LIBS(mmap_enc_mat_init,mmap)
if test "x$ac_cv_search_mmap_enc_mat_init" = "xno"; then
  AC_MSG_ERROR([libmmap not found])
fi

AC_CONFIG_FILES([Makefile mife/Makefile])

AC_OUTPUT
