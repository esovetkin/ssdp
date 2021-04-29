AC_INIT([libssdp], [1.0])
AC_CONFIG_MACRO_DIRS([m4])
AC_PROG_LIBTOOL()
AC_ENABLE_SHARED
AM_ENABLE_SHARED
LT_INIT()
AM_INIT_AUTOMAKE()
AC_PROG_CC
AC_CONFIG_HEADERS([config.h])
AC_CHECK_HEADERS([readline/readline.h, readline/rlconf.h], [HAVE_READLINE_H=yes], [HAVE_READLINE_H=no], [])
  if test "$HAVE_READLINE_H" = "no"; then
    AC_MSG_ERROR([GNU readline development header is not installed.])
  fi
AC_DEFINE([LIBSSDPVERSION], ["1.0"], [libssdp version string])
AC_CONFIG_FILES([
 Makefile
])
AC_OUTPUT
