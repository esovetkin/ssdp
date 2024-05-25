m4_define([ssdp_major_version],    [3])
m4_define([ssdp_minor_version],    [0])
m4_define([ssdp_revis_version],    [0])
m4_define([ssdp_libversion],[ssdp_major_version:ssdp_minor_version:ssdp_revis_version])
m4_define([ssdp_version],[ssdp_major_version.ssdp_minor_version.ssdp_revis_version])

AC_INIT([libssdp],[ssdp_version])
AC_CONFIG_MACRO_DIRS([m4])
LT_INIT([])
AC_PROG_SED()
AC_PROG_AWK()
AC_PROG_EGREP()
AC_ENABLE_SHARED
AC_ENABLE_SHARED([])
LT_INIT([win32-dll])
AM_INIT_AUTOMAKE([subdir-objects])
AC_PROG_CC
AC_CONFIG_SRCDIR([src/sky_dome.c])
AC_CONFIG_HEADERS([src/config.h])
AC_CHECK_LIB(m,main)
AC_CHECK_HEADERS([readline/readline.h, readline/rlconf.h], [HAVE_READLINE_H=yes], [HAVE_READLINE_H=no], [])
  if test "$HAVE_READLINE_H" = "no"; then
    AC_MSG_ERROR([GNU readline development header is not installed.])
  fi
AX_BLAS
AX_LAPACK
AX_LIB_GDAL
AC_CHECK_HEADERS([proj.h], [HAVE_PROJ_H=yes], [HAVE_PROJ_H=no], [])
  if test "$HAVE_PROJ_H" = "no"; then
    AC_MSG_ERROR([PROJ development header is not installed.])
  fi
AC_ARG_ENABLE([fastatan2],
    AS_HELP_STRING([--disable-fastatan2], [Disable the fast atan2 approximation]))

AS_IF([test "x$enable_fastatan2" != "xno"], [
  AC_DEFINE([FAST_ATAN2], [], [fast atan2 approximation])
])

AX_LIB_HDF5([serial])
if test "$with_hdf5" = "yes"; then
        CC="$HDF5_CC"
        #AC_MSG_NOTICE([cc "$CC"])
        #AC_MSG_NOTICE([h5cc "$H5CC"])
        #AC_MSG_NOTICE([libs "$HDF5_LIBS"])
        #AC_MSG_NOTICE([cflags "$HDF5_CPPFLAGS"])
        #AC_MSG_NOTICE([ldflags "$HDF5_LDFLAGS"])
else
        AC_MSG_ERROR([Unable to find HDF5, serial HDF5 required!])
fi
AC_ARG_ENABLE([openmp],
[  --enable-openmp    enable openmp parallelization],
[case "${enableval}" in
  yes) openmp=true ;;
  no)  openmp=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-openmp]) ;;
esac],[openmp=false])
AM_CONDITIONAL([OPENMP], [test x$openmp = xtrue])

# Check for FFTW
AC_CHECK_LIB([fftw3], [fftw_plan_dft_1d], [have_fftw=yes], [have_fftw=no])

if test "$have_fftw" = "yes"; then
  AC_DEFINE([HAVE_FFTW], [1], [Have FFTW library])
  AC_SUBST([FFTW_LIBS], [-lfftw3])
else
  AC_MSG_ERROR([FFTW library not found])
fi

AC_ARG_ENABLE([gendaylit],
[  --enable-gendaylit    enable gendaylit-like Perez model implementation],
[case "${enableval}" in
  yes) gendaylit=true ;;
  no)  gendaylit=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-gendaylit]) ;;
esac],[gendaylit=false])
AS_IF([test "x$gendaylit" == "xtrue"], [
  AC_DEFINE([GENDAYLIT], [], [Use the gendaylit Perez model flavor])
])
AC_SUBST(LIBVERSION, [ssdp_libversion])
AC_DEFINE([LIBSSDPVERSION], ["[ssdp_libversion]"], [libssdp version string])
AC_CONFIG_FILES([
 Makefile src/Makefile tests/Makefile
])
AC_OUTPUT
