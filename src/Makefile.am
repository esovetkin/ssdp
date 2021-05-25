ACLOCAL_AMFLAGS = -I m4
AM_CFLAGS = -Wall

# -- Top-level Makefile.am --
lib_LTLIBRARIES = libssdp.la
libssdp_la_SOURCES = sky_dome.c sky_model.c pv-aoi.c project.c trace.c vector.c ground.c print.c libssdp.c shull.c delaunay.c ll.c error.c sunpos.c fatan2.c sky_dome.h sky_model.h pv-aoi.h project.h trace.h vector.h ground.h print.h libssdp.h shull.h delaunay.h ll.h error.h sunpos.h fatan2.h winunistd.h libssdp.map
CLEANFILES = errorflags.h errormessages.h parsedef.h libssdp_structs.h
libssdp_la_LDFLAGS = -lm -Wl,--version-script=${srcdir}/libssdp.map -no-undefined
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
ssdp_SOURCES = lio.c variables.c parser.c configparser.c arrayparser.c topoparser.c simparser.c readlineshell.c util.c ssdp.c lio.h variables.h parser.h parserutil.h readlineshell.h util.h winunistd.h gen_parserflags.sh gen_errorflags.sh collect_datastructs.sh gen_ssdpdoc.sh 
ssdp_LDADD = libssdp.la -lreadline
man1_MANS = ssdp.man

BUILT_SOURCES = errorflags.h errormessages.h gen_errorflags.sh libssdp_structs.h collect_datastructs.sh parsedef.h gen_parserflags.sh ssdp.man
MOSTLYCLEANFILES = parsedef.h errorflags.h errormessages.h libssdp_structs.h ssdp.man
parsedef.h: gen_parserflags.sh parser.c configparser.c arrayparser.c topoparser.c simparser.c
	${SHELL} ${srcdir}/gen_parserflags.sh ${srcdir} parser.c configparser.c arrayparser.c topoparser.c simparser.c
errorflags.h errormessages.h: gen_errorflags.sh $(libssdp_la_SOURCES)
	${SHELL} ${srcdir}/gen_errorflags.sh ${srcdir} $(libssdp_la_SOURCES)
libssdp_structs.h: collect_datastructs.sh vector.h sky_dome.h delaunay.h trace.h ground.h project.h libssdp.c
	${SHELL} ${srcdir}/collect_datastructs.sh ${srcdir} vector.h  sky_dome.h delaunay.h trace.h ground.h project.h libssdp.c
ssdp.man: gen_ssdpdoc.sh parser.c configparser.c arrayparser.c topoparser.c simparser.c
	${SHELL} ${srcdir}/gen_ssdpdoc.sh ${srcdir} parser.c configparser.c arrayparser.c topoparser.c simparser.c

