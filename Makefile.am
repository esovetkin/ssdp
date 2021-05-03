ACLOCAL_AMFLAGS = -I m4
AM_CFLAGS = -Wall

# -- Top-level Makefile.am --
lib_LTLIBRARIES = libssdp.la
libssdp_la_SOURCES = sky_dome.c sky_model.c pv-aoi.c project.c trace.c vector.c ground.c print.c libssdp.c shull.c delaunay.c ll.c error.c sunpos.c sky_dome.h sky_model.h pv-aoi.h project.h trace.h vector.h ground.h print.h libssdp.h shull.h delaunay.h ll.h error.h sunpos.h
CLEANFILES = errorflags.h errormessages.h parsedef.h libssdp_structs.h
libssdp_la_LDFLAGS = -lm -Wl,--version-script=libssdp.map
include_HEADERS = libssdp.h
libssdp_la_CFLAGS = $(AM_CFLAGS) -DLL_CIRCULAR -D_GNU_SOURCE


bin_PROGRAMS = ssdp
ssdp_SOURCES = io.c variables.c parser.c configparser.c arrayparser.c topoparser.c simparser.c readlineshell.c util.c ssdp.c io.h variables.h parser.h parserutil.h readlineshell.h util.h 
ssdp_LDADD = libssdp.la -lreadline
man1_MANS = ssdp.man

BUILT_SOURCES = errorflags.h errormessages.h gen_errorflags.sh libssdp_structs.h collect_datastructs.sh parsedef.h gen_parserflags.sh ssdp.man

parsedef.h: gen_parserflags.sh parser.c configparser.c arrayparser.c topoparser.c simparser.c
	cd ${top_srcdir} && ${SHELL} gen_parserflags.sh parser.c configparser.c arrayparser.c topoparser.c simparser.c
errorflags.h errormessages.h: gen_errorflags.sh $(libssdp_la_SOURCES)
	cd ${top_srcdir} && ${SHELL} gen_errorflags.sh  $(libssdp_la_SOURCES)
libssdp_structs.h: collect_datastructs.sh vector.h sky_dome.h delaunay.h trace.h ground.h project.h libssdp.c
	cd ${top_srcdir} && ${SHELL} collect_datastructs.sh vector.h sky_dome.h delaunay.h trace.h ground.h project.h libssdp.c
ssdp.man: gen_ssdpdoc.sh parser.c configparser.c arrayparser.c topoparser.c simparser.c
	cd ${top_srcdir} && ${SHELL} gen_ssdpdoc.sh parser.c configparser.c arrayparser.c topoparser.c simparser.c

