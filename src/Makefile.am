AUTOMAKE_OPTIONS = subdir-objects

ACLOCAL_AMFLAGS = -I m4
AM_CFLAGS = -Wall


lib_LTLIBRARIES = libssdp.la
libssdp_la_SOURCES = sky_dome.c sky_model.c pv-aoi.c project.c trace.c vector.c ground.c print.c libssdp.c shull.c delaunay.c ll.c error.c sunpos.c freespa/freespa.c filterimage/filter.c filterimage/floatimage.c fatan2.c coordinates.c epsg.c edt.c minfill.c topogdal.c pngout.c h5io.c iset.c rtree.c sobolseq.c random_fail_malloc.c sky_dome.h sky_model.h pv-aoi.h project.h trace.h vector.h ground.h print.h libssdp.h shull.h delaunay.h ll.h error.h sunpos.h fatan2.h winunistd.h coordinates.h epsg.h edt.h minfill.h topogdal.h pngout.h h5io.h iset.h rtree.h sobolseq.h soboldata.h random_fail_malloc.h libssdp.map
CLEANFILES = errorflags.h errormessages.h parsedef.h libssdp_structs.h
libssdp_la_LDFLAGS = -lm -lproj -lpng $(LAPACK_LIBS) $(BLAS_LIBS) $(FFTW_LIBS) $(GDAL_LDFLAGS) -Wl,--version-script=${srcdir}/libssdp.map -no-undefined -version-number $(LIBVERSION)
include_HEADERS = libssdp.h libssdp_structs.h
libssdp_la_CFLAGS = $(AM_CFLAGS) $(GDAL_CFLAGS) $(HDF5_CPPFLAGS) -DLL_CIRCULAR -D_GNU_SOURCE -D'ERRORSTATE=ssdp_error_state'
sky_dome.$(OBJEXT): errorflags.h
sky_model.$(OBJEXT): errorflags.h
project.$(OBJEXT): errorflags.h
trace.$(OBJEXT): errorflags.h
ground.$(OBJEXT): errorflags.h
shull.$(OBJEXT): errorflags.h
delaunay.$(OBJEXT): errorflags.h
ll.$(OBJEXT): errorflags.h
error.$(OBJEXT): errorflags.h errormessages.h
freespa/freespa.$(OBJEXT):
sunpos.$(OBJEXT): errorflags.h freespa/freespa.h
nodist_libssdp_la_SOURCES = errorflags.h errormessages.h libssdp_structs.h


bin_PROGRAMS = ssdp
ssdp_SOURCES = lio.c variables.c parser.c configparser.c arrayparser.c topoparser.c simparser.c readlineshell.c util.c epsg.c coordinates.c pngout.c h5io.c iset.c rtree.c sobolseq.c random_fail_malloc.c ssdp.c lio.h variables.h parser.h parserutil.h readlineshell.h util.h winunistd.h epsg.h coordinates.h gen_parserflags.sh gen_errorflags.sh collect_datastructs.sh gen_ssdpdoc.sh
if OPENMP
ssdp_LDFLAGS= -lreadline -fopenmp -lm -lproj -lpng $(LAPACK_LIBS) $(BLAS_LIBS) $(FFTW_LIBS) $(HDF5_LIBS) $(GDAL_LDFLAGS) $(HDF5_LDFLAGS)
else
ssdp_LDFLAGS= -lreadline -lm -lproj -lpng $(LAPACK_LIBS) $(BLAS_LIBS) $(FFTW_LIBS) $(HDF5_LIBS) $(GDAL_LDFLAGS) $(HDF5_LDFLAGS)
endif
ssdp_CFLAGS = $(AM_CFLAGS) $(GDAL_CFLAGS) $(HDF5_CPPFLAGS) -DOPENMP -fopenmp -D'ERRORSTATE=ssdp_error_state'
ssdp_LDADD = libssdp.la
man1_MANS = ssdp.man
# scripts in extra dist
EXTRA_DIST = gen_errorflags.sh collect_datastructs.sh gen_parserflags.sh freespa/freespa.c freespa/freespa.h freespa/freespa_tables.h freespa/freespa_dt_table.h filterimage/filter.c filterimage/floatimage.c filterimage/filter.h filterimage/floatimage.h

BUILT_SOURCES = errorflags.h errormessages.h libssdp_structs.h parsedef.h ssdp.man

MOSTLYCLEANFILES = parsedef.h errorflags.h errormessages.h libssdp_structs.h ssdp.man
parsedef.h: gen_parserflags.sh parser.c configparser.c arrayparser.c topoparser.c simparser.c
	${SHELL} ${srcdir}/gen_parserflags.sh ${srcdir} parser.c configparser.c arrayparser.c topoparser.c simparser.c
# errorflags and errormessages are created with the same script
# when using parallel make we must prevent make trying to build both headers at the same time
# as this leads to race conditions, this ensures make will not try to build both at the same time
errorflags.h: gen_errorflags.sh $(libssdp_la_SOURCES)
	${SHELL} ${srcdir}/gen_errorflags.sh ${srcdir} $(libssdp_la_SOURCES) filterimage/floatimage_io.c
errormessages.h: errorflags.h
	@test -f $@ || rm -f errorflags.h
	@test -f $@ || $(MAKE) $(AM_MAKEFLAGS) errorflags.h
libssdp_structs.h: collect_datastructs.sh vector.h sky_dome.h delaunay.h trace.h ground.h project.h libssdp.c
	${SHELL} ${srcdir}/collect_datastructs.sh ${srcdir} vector.h  sky_dome.h delaunay.h trace.h ground.h project.h libssdp.c
ssdp.man: gen_ssdpdoc.sh parser.c configparser.c arrayparser.c topoparser.c simparser.c
	${SHELL} ${srcdir}/gen_ssdpdoc.sh ${srcdir} parser.c configparser.c arrayparser.c topoparser.c simparser.c

