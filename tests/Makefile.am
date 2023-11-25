bin_PROGRAMS = ssdptest_edt ssdptest_minfill ssdptest_epsg ssdptest_coords ssdptest_topogdal ssdptest_pngout ssdptest_h5io ssdptest_iset ssdptest_rtree ssdptest_sobolseq ssdptest_lcg

ssdptest_edt_SOURCES = $(top_builddir)/src/edt.c
ssdptest_edt_CFLAGS = -DRUNTEST
ssdptest_edt_LDADD = libfmem.la

ssdptest_minfill_SOURCES = $(top_builddir)/src/minfill.c
ssdptest_minfill_CFLAGS = -DRUNTEST
ssdptest_minfill_LDADD = libiset.la libfmem.la

ssdptest_epsg_SOURCES =  $(top_builddir)/src/epsg.c
ssdptest_epsg_CFLAGS = -DRUNTEST
ssdptest_epsg_LDFLAGS = -lproj
ssdptest_epsg_LDADD = libcoordinates.la libfmem.la

ssdptest_pngout_SOURCES = $(top_builddir)/src/pngout.c
ssdptest_pngout_CFLAGS = -DRUNTEST
ssdptest_pngout_LDFLAGS = -lpng
ssdptest_pngout_LDADD = libfmem.la

ssdptest_h5io_SOURCES = $(top_builddir)/src/h5io.c
ssdptest_h5io_CFLAGS = -DRUNTEST $(HDF5_CPPFLAGS)
ssdptest_h5io_LDFLAGS = $(HDF5_LIBS) $(HDF5_LDFLAGS)
ssdptest_h5io_LDADD = libfmem.la

ssdptest_iset_SOURCES = $(top_builddir)/src/iset.c
ssdptest_iset_CFLAGS = -DRUNTEST
ssdptest_iset_LDADD = libfmem.la

ssdptest_rtree_SOURCES = $(top_builddir)/src/rtree.c
ssdptest_rtree_CFLAGS = -DRUNTEST

ssdptest_sobolseq_SOURCES = $(top_builddir)/src/sobolseq.c
ssdptest_sobolseq_CFLAGS = -DRUNTEST

ssdptest_lcg_SOURCES = $(top_builddir)/src/lcg.c
ssdptest_lcg_CFLAGS = -fopenmp -DRUNTEST

ssdptest_coords_SOURCES = $(top_builddir)/src/coordinates.c
ssdptest_coords_CFLAGS = -DRUNTEST
ssdptest_coords_LDFLAGS = -lproj
ssdptest_coords_LDADD = libepsg.la libfmem.la

ssdptest_topogdal_SOURCES = $(top_builddir)/src/topogdal.c $(top_builddir)/src/error.c
ssdptest_topogdal_CFLAGS = -DRUNTEST $(GDAL_CFLAGS)
ssdptest_topogdal_LDFLAGS = -lproj $(GDAL_LDFLAGS)
ssdptest_topogdal_LDADD = libedt.la libepsg.la libcoordinates.la libfmem.la

noinst_LTLIBRARIES = libcoordinates.la libepsg.la libedt.la libiset.la libfmem.la
libcoordinates_la_SOURCES = $(top_builddir)/src/coordinates.c
libepsg_la_SOURCES = $(top_builddir)/src/epsg.c
libedt_la_SOURCES = $(top_builddir)/src/edt.c
libiset_la_SOURCES = $(top_builddir)/src/iset.c
libfmem_la_SOURCES = $(top_builddir)/src/random_fail_malloc.c

TESTS=ssdptest_edt ssdptest_minfill ssdptest_epsg ssdptest_coords ssdptest_topogdal ssdptest_pngout ssdptest_h5io ssdptest_iset ssdptest_rtree ssdptest_sobolseq ssdptest_lcg
