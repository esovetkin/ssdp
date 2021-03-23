libSRC=ssdp.c sky_dome.c sky_model.c project.c vector.c ground.c util.c libssdp.c shull.c delaunay.c ll.c
libOBJ=ssdp.o sky_dome.o sky_model.o project.o vector.o ground.o util.o libssdp.o shull.o delaunay.o ll.o
libHDR=       sky_dome.h sky_model.h project.h vector.h ground.h util.h libssdp.h shull.h delaunay.h ll.h

SRC=io.c ssdp.c
OBJ=io.o ssdp.o
HDR=io.h
CC=gcc
# CFLAGS=-Og -flto -g -Wall -pedantic -fPIC
CFLAGS=-O3 -flto -fPIC
LFLAGS=-lm -L.
libssdp_la_LDFLAGS = -lm -Wl,--version-script=libssdp.map

all: ssdp
ssdp: libssdp.so $(OBJ)
	$(CC) ssdp.c util.o io.o -o ssdp $(LFLAGS) -lssdp
libssdp.so: Makefile libssdp.map $(libOBJ)
	$(CC) -shared -o libssdp.so $(libOBJ) $(libssdp_la_LDFLAGS)
libssdp.o: Makefile libssdp.c $(libHDR)
project.o:project.h sky_dome.h vector.h util.h
io.o:project.h sky_dome.h vector.h util.h
sky_dome.o:sky_dome.h vector.h util.h
sky_model.o:sky_dome.h util.h
vector.o:sky_dome.h util.h
ground.o:sky_dome.h vector.h util.h
delaunay.o:ll.h shull.h
shull.o:ll.h shull.h
	$(CC) -c -DLL_CIRCULAR -D_GNU_SOURCE $(CFLAGS) -o shull.o shull.c
ll.o:ll.h
	$(CC) -c -DLL_CIRCULAR -D_GNU_SOURCE $(CFLAGS) -o ll.o ll.c
clean:
	rm *.o ssdp
