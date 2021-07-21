m4_define([ssdp_major_version],    [2])
m4_define([ssdp_minor_version],    [2])
m4_define([ssdp_revis_version],    [4])
m4_define([ssdp_version],[ssdp_major_version.ssdp_minor_version.ssdp_revis_version])
AC_INIT([libssdp], [ssdp_version])
AC_CONFIG_MACRO_DIRS([m4])
AC_PROG_LIBTOOL()
AC_PROG_SED()
AC_PROG_AWK()
AC_PROG_EGREP()
AC_ENABLE_SHARED
AM_ENABLE_SHARED
LT_INIT([win32-dll])
AM_INIT_AUTOMAKE()
AC_PROG_CC
AC_CONFIG_SRCDIR([src/sky_dome.c])
AC_CONFIG_HEADERS([src/config.h])
AC_CHECK_HEADERS([readline/readline.h, readline/rlconf.h], [HAVE_READLINE_H=yes], [HAVE_READLINE_H=no], [])
  if test "$HAVE_READLINE_H" = "no"; then
    AC_MSG_ERROR([GNU readline development header is not installed.])
  fi
AC_ARG_ENABLE([fastatan2],
    AS_HELP_STRING([--disable-fastatan2], [Disable the fast atan2 approximation]))

AS_IF([test "x$enable_fastatan2" != "xno"], [
  AC_DEFINE([FAST_ATAN2], [], [fast atan2 approximation])
])
AC_ARG_ENABLE([openmp],
[  --enable-openmp    enable openmp parallelization],
[case "${enableval}" in
  yes) openmp=true ;;
  no)  openmp=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-openmp]) ;;
esac],[openmp=false])
AM_CONDITIONAL([OPENMP], [test x$openmp = xtrue])



AC_DEFINE([LIBSSDPVERSION], ["[ssdp_version]"], [libssdp version string])
AC_CONFIG_FILES([
 Makefile src/Makefile
])
AC_OUTPUT
