AC_INIT([libssdp], [1.1])
AC_CONFIG_MACRO_DIRS([m4])
AC_PROG_LIBTOOL()
AC_ENABLE_SHARED
AM_ENABLE_SHARED
LT_INIT()
AM_INIT_AUTOMAKE()
AC_PROG_CC
AC_CONFIG_SRCDIR([src/sky_dome.c])
AC_CONFIG_HEADERS([src/config.h])
AC_CHECK_HEADERS([readline/readline.h, readline/rlconf.h], [HAVE_READLINE_H=yes], [HAVE_READLINE_H=no], [])
  if test "$HAVE_READLINE_H" = "no"; then
    AC_MSG_ERROR([GNU readline development header is not installed.])
  fi
AC_DEFINE([LIBSSDPVERSION], ["1.1"], [libssdp version string])
AC_CONFIG_FILES([
 Makefile src/Makefile
])
AC_OUTPUT
