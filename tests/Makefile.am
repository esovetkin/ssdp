AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = ssdptest_edt ssdptest_minfill ssdptest_epsg ssdptest_coords ssdptest_topogdal ssdptest_pngout ssdptest_h5io ssdptest_hashmap ssdptest_rtree ssdptest_sobolseq ssdptest_lcg ssdptest_skydome

ssdptest_edt_SOURCES = $(top_builddir)/src/edt.c
ssdptest_edt_CFLAGS = -DRUNTEST

ssdptest_minfill_SOURCES = $(top_builddir)/src/minfill.c
ssdptest_minfill_CFLAGS = -DRUNTEST
ssdptest_minfill_LDADD = libhashmap.la

ssdptest_epsg_SOURCES =  $(top_builddir)/src/epsg.c
ssdptest_epsg_CFLAGS = -DRUNTEST
ssdptest_epsg_LDADD = libcoordinates.la -lproj

ssdptest_pngout_SOURCES = $(top_builddir)/src/pngout.c
ssdptest_pngout_CFLAGS = -DRUNTEST
ssdptest_pngout_LDFLAGS = -lpng

ssdptest_h5io_SOURCES = $(top_builddir)/src/h5io.c
ssdptest_h5io_CFLAGS = -DRUNTEST $(HDF5_CPPFLAGS)
ssdptest_h5io_LDFLAGS = $(HDF5_LIBS) $(HDF5_LDFLAGS)

ssdptest_hashmap_SOURCES = $(top_builddir)/src/hashmap.c $(top_builddir)/src/murmurhash3.c
ssdptest_hashmap_CFLAGS = -DRUNTEST

ssdptest_rtree_SOURCES = $(top_builddir)/src/rtree.c
ssdptest_rtree_CFLAGS = -DRUNTEST

ssdptest_sobolseq_SOURCES = $(top_builddir)/src/sobolseq.c
ssdptest_sobolseq_CFLAGS = -DRUNTEST

ssdptest_lcg_SOURCES = $(top_builddir)/src/lcg.c
ssdptest_lcg_CFLAGS = -fopenmp -DRUNTEST

ssdptest_coords_SOURCES = $(top_builddir)/src/coordinates.c
ssdptest_coords_CFLAGS = -DRUNTEST
ssdptest_coords_LDADD = libepsg.la -lproj

ssdptest_topogdal_SOURCES = $(top_builddir)/src/topogdal.c $(top_builddir)/src/error.c
ssdptest_topogdal_CFLAGS = -fopenmp -DRUNTEST $(GDAL_CFLAGS)
ssdptest_topogdal_LDADD = libedt.la libepsg.la libcoordinates.la -lproj $(GDAL_LDFLAGS)

ssdptest_skydome_SOURCES = $(top_builddir)/src/sky_dome.c $(top_builddir)/src/error.c $(top_builddir)/src/vector.c $(top_builddir)/src/print.c $(top_builddir)/src/fatan2.c
ssdptest_skydome_CFLAGS = -DRUNTEST $(HDF5_CPPFLAGS)
ssdptest_skydome_LDFLAGS = $(HDF5_LIBS) $(HDF5_LDFLAGS)
ssdptest_skydome_LDADD = libh5io.la

noinst_LTLIBRARIES = libcoordinates.la libepsg.la libedt.la libhashmap.la libh5io.la
libcoordinates_la_SOURCES = $(top_builddir)/src/coordinates.c
libepsg_la_SOURCES = $(top_builddir)/src/epsg.c
libedt_la_SOURCES = $(top_builddir)/src/edt.c
libh5io_la_SOURCES = $(top_builddir)/src/h5io.c
libhashmap_la_SOURCES = $(top_builddir)/src/hashmap.c $(top_builddir)/src/murmurhash3.c

TESTS=ssdptest_edt ssdptest_minfill ssdptest_epsg ssdptest_coords ssdptest_topogdal ssdptest_pngout ssdptest_h5io ssdptest_hashmap ssdptest_rtree ssdptest_sobolseq ssdptest_lcg ssdptest_skydome
