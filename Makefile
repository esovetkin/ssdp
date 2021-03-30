libSRC=sky_dome.c sky_model.c pv-aoi.c project.c vector.c ground.c util.c libssdp.c shull.c delaunay.c ll.c error.c sunpos.c
libOBJ=sky_dome.o sky_model.o pv-aoi.o project.o vector.o ground.o util.o libssdp.o shull.o delaunay.o ll.o error.o sunpos.o
libHDR=sky_dome.h sky_model.h pv-aoi.h project.h vector.h ground.h util.h libssdp.h shull.h delaunay.h ll.h error.h sunpos.h

SRC=io.c ssdp.c
OBJ=io.o ssdp.o
HDR=io.h
CC=gcc
# CC=clang
# CFLAGS=-Og -flto -g -Wall -pedantic -fPIC
CFLAGS=-O3 -flto -fPIC -march=native
# Ofast can speed up code by some 20% over O3, however, do not use 
# -funsafe-math-optimizations, it breaks the code.
# if you feel lucky you can use:
# CFLAGS=-Ofast -flto -fno-unsafe-math-optimizations -march=native -fPIC
LFLAGS=-lm -L.
libssdp_la_LDFLAGS = -lm -Wl,--version-script=libssdp.map

all: ssdp
ssdp: libssdp.so $(OBJ)
	$(CC) ssdp.c util.o io.o -o ssdp $(LFLAGS) -lssdp
io.o:project.h sky_dome.h vector.h util.h
libssdp.so: Makefile libssdp.map $(libOBJ)
	$(CC) -shared -o libssdp.so $(libOBJ) $(libssdp_la_LDFLAGS)
libssdp.o: Makefile libssdp.c $(libHDR) error.h errorflags.h
pv-aoi.o:pv-aoi.c
project.o:pv-aoi.h project.h sky_dome.h vector.h util.h error.h errorflags.h project.c
sky_dome.o:sky_dome.h vector.h util.h error.h errorflags.h sky_dome.c
sky_model.o:sky_dome.h util.h error.h errorflags.h sky_model.c
vector.o:sky_dome.h util.h error.h errorflags.h vector.c
ground.o:sky_dome.h vector.h util.h errorflags.h ground.c
delaunay.o:ll.h shull.h error.h errorflags.h delaunay.c
shull.o:ll.h shull.h error.h errorflags.h shull.c 
	$(CC) -c -DLL_CIRCULAR -D_GNU_SOURCE $(CFLAGS) -o shull.o shull.c
ll.o:ll.h error.h ll.c
	$(CC) -c -DLL_CIRCULAR -D_GNU_SOURCE $(CFLAGS) -o ll.o ll.c
error.o: error.h errorflags.h errormessages.h
errorflags.h errormessages.h: gen_errorflags.sh $(libSRC)
	${SHELL} gen_errorflags.sh  $(libSRC)
clean:
	rm *.o ssdp libssdp.so errorflags.h errormessages.h
