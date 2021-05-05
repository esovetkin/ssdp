ACLOCAL_AMFLAGS = -I m4
AM_CFLAGS = -Wall

# -- Top-level Makefile.am --
lib_LTLIBRARIES = libssdp.la
libssdp_la_SOURCES = sky_dome.c sky_model.c pv-aoi.c project.c trace.c vector.c ground.c print.c libssdp.c shull.c delaunay.c ll.c error.c sunpos.c sky_dome.h sky_model.h pv-aoi.h project.h trace.h vector.h ground.h print.h libssdp.h shull.h delaunay.h ll.h error.h sunpos.h libssdp.map
CLEANFILES = errorflags.h errormessages.h parsedef.h libssdp_structs.h
libssdp_la_LDFLAGS = -lm -Wl,--version-script=libssdp.map
include_HEADERS = libssdp.h libssdp_structs.h
libssdp_la_CFLAGS = $(AM_CFLAGS) -DLL_CIRCULAR -D_GNU_SOURCE
sky_dome.$(OBJEXT): errorflags.h
sky_model.$(OBJEXT): errorflags.h
project.$(OBJEXT): errorflags.h
trace.$(OBJEXT): errorflags.h
ground.$(OBJEXT): errorflags.h
shull.$(OBJEXT): errorflags.h
delaunay.$(OBJEXT): errorflags.h
ll.$(OBJEXT): errorflags.h
error.$(OBJEXT): errorflags.h errormessages.h
sunpos.$(OBJEXT): errorflags.h
nodist_libssdp_la_SOURCES = errorflags.h errormessages.h libssdp_structs.h


bin_PROGRAMS = ssdp
ssdp_SOURCES = io.c variables.c parser.c configparser.c arrayparser.c topoparser.c simparser.c readlineshell.c util.c ssdp.c io.h variables.h parser.h parserutil.h readlineshell.h util.h gen_parserflags.sh gen_errorflags.sh collect_datastructs.sh gen_ssdpdoc.sh 
ssdp_LDADD = libssdp.la -lreadline
man1_MANS = ssdp.man

BUILT_SOURCES = errorflags.h errormessages.h gen_errorflags.sh libssdp_structs.h collect_datastructs.sh parsedef.h gen_parserflags.sh ssdp.man
MOSTLYCLEANFILES = parsedef.h errorflags.h errormessages.h libssdp_structs.h ssdp.man
parsedef.h: gen_parserflags.sh parser.c configparser.c arrayparser.c topoparser.c simparser.c
	${SHELL} ${top_srcdir}/gen_parserflags.sh ${top_srcdir}/parser.c ${top_srcdir}/configparser.c ${top_srcdir}/arrayparser.c ${top_srcdir}/topoparser.c ${top_srcdir}/simparser.c
errorflags.h errormessages.h: gen_errorflags.sh $(libssdp_la_SOURCES)
	${SHELL} ${top_srcdir}/gen_errorflags.sh ${top_srcdir}/errorflags.h ${top_srcdir}/errormessages.h $(libssdp_la_SOURCES)
libssdp_structs.h: collect_datastructs.sh vector.h sky_dome.h delaunay.h trace.h ground.h project.h libssdp.c
	${SHELL} ${top_srcdir}/collect_datastructs.sh ${top_srcdir}/vector.h ${top_srcdir}/sky_dome.h ${top_srcdir}/delaunay.h ${top_srcdir}/trace.h ${top_srcdir}/ground.h ${top_srcdir}/project.h ${top_srcdir}/libssdp.c
ssdp.man: gen_ssdpdoc.sh parser.c configparser.c arrayparser.c topoparser.c simparser.c
	${SHELL} ${top_srcdir}/gen_ssdpdoc.sh parser.c ${top_srcdir}/configparser.c ${top_srcdir}/arrayparser.c ${top_srcdir}/topoparser.c ${top_srcdir}/simparser.c

