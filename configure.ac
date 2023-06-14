m4_define([ssdp_major_version],    [3])
m4_define([ssdp_minor_version],    [0])
m4_define([ssdp_revis_version],    [0])
m4_define([ssdp_libversion],[ssdp_major_version:ssdp_minor_version:ssdp_revis_version])
m4_define([ssdp_version],[ssdp_major_version.ssdp_minor_version.ssdp_revis_version])

AC_DEFUN([AX_LIB_GDAL],
[
    dnl If gdal-config path is not given in ---with-gdal option,
    dnl check if it is present in the system anyway
    AC_ARG_WITH([gdal],
        AS_HELP_STRING([--with-gdal=@<:@ARG@:>@],
            [Specify full path to gdal-config script]),
        [ac_gdal_config_path=$withval],
        [gdal_config_system=check])

    dnl if gdal-config is present in the system, fill the ac_gdal_config_path variable with it full path
    AS_IF([test "x$gdal_config_system" = xcheck],
          [AC_PATH_PROG([GDAL_CONFIG], [gdal-config])],
          [AC_PATH_PROG([GDAL_CONFIG], [gdal-config],
              [no], [`dirname $ac_gdal_config_path 2> /dev/null`])]
    )

    if test ! -x "$GDAL_CONFIG"; then
        AC_MSG_ERROR([gdal-config does not exist or it is not an executable file])
            GDAL_CONFIG="no"
            found_gdal="no"
    fi

    GDAL_VERSION=""
    GDAL_CFLAGS=""
    GDAL_LDFLAGS=""
    GDAL_DEP_LDFLAGS=""
    GDAL_OGR_ENABLED=""


    dnl
    dnl Check GDAL library (libgdal)
    dnl

    if test "$GDAL_CONFIG" != "no"; then
        AC_MSG_CHECKING([for GDAL library])

        GDAL_VERSION="`$GDAL_CONFIG --version`"
        GDAL_CFLAGS="`$GDAL_CONFIG --cflags`"
        GDAL_LDFLAGS="`$GDAL_CONFIG --libs`"
        GDAL_DEP_LDFLAGS="`$GDAL_CONFIG --dep-libs`"

        AC_DEFINE([HAVE_GDAL], [1], [Define to 1 if GDAL library are available])

        found_gdal="yes"
    else
        found_gdal="no"
    fi

    AC_MSG_RESULT([$found_gdal])

    if test "$found_gdal" = "yes"; then
        AC_MSG_CHECKING([for OGR support in GDAL library])

        GDAL_OGR_ENABLED="`$GDAL_CONFIG --ogr-enabled`"
        AC_DEFINE([HAVE_GDAL_OGR], [1], [Define to 1 if GDAL library includes OGR support])

        AC_MSG_RESULT([$GDAL_OGR_ENABLED])
    fi

    dnl
    dnl Check if required version of GDAL is available
    dnl

    gdal_version_req=ifelse([$1], [], [], [$1])
    if test "$found_gdal" = "yes" -a -n "$gdal_version_req"; then

        AC_MSG_CHECKING([if GDAL version is >= $gdal_version_req])

        dnl Decompose required version string of GDAL
        dnl and calculate its number representation
        gdal_version_req_major=`expr $gdal_version_req : '\([[0-9]]*\)'`
        gdal_version_req_minor=`expr $gdal_version_req : '[[0-9]]*\.\([[0-9]]*\)'`
        gdal_version_req_micro=`expr $gdal_version_req : '[[0-9]]*\.[[0-9]]*\.\([[0-9]]*\)'`
        if test "x$gdal_version_req_micro" = "x"; then
            gdal_version_req_micro="0"
        fi

        gdal_version_req_number=`expr $gdal_version_req_major \* 1000000 \
                                   \+ $gdal_version_req_minor \* 1000 \
                                   \+ $gdal_version_req_micro`

        dnl Decompose version string of installed GDAL
        dnl and calculate its number representation
        gdal_version_major=`expr $GDAL_VERSION : '\([[0-9]]*\)'`
        gdal_version_minor=`expr $GDAL_VERSION : '[[0-9]]*\.\([[0-9]]*\)'`
        gdal_version_micro=`expr $GDAL_VERSION : '[[0-9]]*\.[[0-9]]*\.\([[0-9]]*\)'`
        if test "x$gdal_version_micro" = "x"; then
            gdal_version_micro="0"
        fi

        gdal_version_number=`expr $gdal_version_major \* 1000000 \
                                   \+ $gdal_version_minor \* 1000 \
                                   \+ $gdal_version_micro`

        gdal_version_check=`expr $gdal_version_number \>\= $gdal_version_req_number`
        if test "$gdal_version_check" = "1"; then
            AC_MSG_RESULT([yes])
        else
            AC_MSG_RESULT([no])
            AC_MSG_ERROR([GDAL $GDAL_VERSION found, but required version is $gdal_version_req])
        fi
    fi

    AC_SUBST(GDAL_VERSION)
    AC_SUBST(GDAL_CFLAGS)
    AC_SUBST(GDAL_LDFLAGS)
    AC_SUBST(GDAL_DEP_LDFLAGS)
    AC_SUBST(GDAL_OGR_ENABLED)
])

AC_INIT([libssdp],[ssdp_version])
AC_CONFIG_MACRO_DIRS([m4])
LT_INIT([])
AC_PROG_SED()
AC_PROG_AWK()
AC_PROG_EGREP()
AC_ENABLE_SHARED
AC_ENABLE_SHARED([])
LT_INIT([win32-dll])
AM_INIT_AUTOMAKE()
AC_PROG_CC
AC_CONFIG_SRCDIR([src/sky_dome.c])
AC_CONFIG_HEADERS([src/config.h])
AC_CHECK_LIB(m,main)
AC_CHECK_HEADERS([readline/readline.h, readline/rlconf.h], [HAVE_READLINE_H=yes], [HAVE_READLINE_H=no], [])
  if test "$HAVE_READLINE_H" = "no"; then
    AC_MSG_ERROR([GNU readline development header is not installed.])
  fi
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
