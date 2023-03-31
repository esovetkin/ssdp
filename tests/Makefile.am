bin_PROGRAMS = edt epsg coords topogdal

edt_SOURCES = $(top_builddir)/src/edt.c
edt_CFLAGS = -DRUNTEST

epsg_SOURCES =  $(top_builddir)/src/epsg.c
epsg_CFLAGS = -DRUNTEST
epsg_LDFLAGS = -lproj
epsg_LDADD = libcoordinates.la

coords_SOURCES = $(top_builddir)/src/coordinates.c
coords_CFLAGS = -DRUNTEST
coords_LDFLAGS = -lproj
coords_LDADD = libepsg.la

topogdal_SOURCES = $(top_builddir)/src/topogdal.c $(top_builddir)/src/error.c
topogdal_CFLAGS = -DRUNTEST
topogdal_LDFLAGS = -lproj $(GDAL_LDFLAGS)
topogdal_LDADD = libepsg.la libcoordinates.la

noinst_LTLIBRARIES = libcoordinates.la libepsg.la
libcoordinates_la_SOURCES = $(top_builddir)/src/coordinates.c
libepsg_la_SOURCES = $(top_builddir)/src/epsg.c

TESTS=edt epsg coords topogdal
