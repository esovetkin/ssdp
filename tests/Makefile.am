bin_PROGRAMS = edt epsg coords topogdal pngout h5io

edt_SOURCES = $(top_builddir)/src/edt.c
edt_CFLAGS = -DRUNTEST

epsg_SOURCES =  $(top_builddir)/src/epsg.c
epsg_CFLAGS = -DRUNTEST
epsg_LDFLAGS = -lproj
epsg_LDADD = libcoordinates.la

pngout_SOURCES = $(top_builddir)/src/pngout.c
pngout_CFLAGS = -DRUNTEST
pngout_LDFLAGS = -lpng

h5io_SOURCES = $(top_builddir)/src/h5io.c
h5io_CFLAGS = -DRUNTEST $(HDF5_CPPFLAGS)
h5io_LDFLAGS = $(HDF5_LIBS) $(HDF5_LDFLAGS)

coords_SOURCES = $(top_builddir)/src/coordinates.c
coords_CFLAGS = -DRUNTEST
coords_LDFLAGS = -lproj
coords_LDADD = libepsg.la

topogdal_SOURCES = $(top_builddir)/src/topogdal.c $(top_builddir)/src/error.c
topogdal_CFLAGS = -DRUNTEST $(GDAL_CFLAGS)
topogdal_LDFLAGS = -lproj $(GDAL_LDFLAGS)
topogdal_LDADD = libedt.la libepsg.la libcoordinates.la

noinst_LTLIBRARIES = libcoordinates.la libepsg.la libedt.la
libcoordinates_la_SOURCES = $(top_builddir)/src/coordinates.c
libepsg_la_SOURCES = $(top_builddir)/src/epsg.c
libedt_la_SOURCES = $(top_builddir)/src/edt.c

TESTS=edt epsg coords topogdal pngout h5io
