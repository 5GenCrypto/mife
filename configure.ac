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

AC_SEARCH_LIBS(mpz_urandomm_aes,aesrand)
if test "x$ac_cv_search_mpz_urandomm_aes" = "xno"; then
  AC_MSG_ERROR([libaesrand not found])
fi
AC_SEARCH_LIBS(mmap_enc_mat_init,mmap)
if test "x$ac_cv_search_mmap_enc_mat_init" = "xno"; then
  AC_MSG_ERROR([libmmap not found])
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

AC_CONFIG_FILES([Makefile mife/Makefile])

AC_OUTPUT
