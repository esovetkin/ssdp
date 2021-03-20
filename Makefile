SRC=ssdp.c sky_dome.c sky_model.c project.c vector.c ground.c io.c util.c libssdp.c
OBJ=ssdp.o sky_dome.o sky_model.o project.o vector.o ground.o io.o util.o libssdp.o
HDR=       sky_dome.h sky_model.h project.h vector.h ground.h io.h util.h libssdp.h
CC=gcc
# CFLAGS=-Og -flto -g -Wall -pedantic -fPIC
CFLAGS=-O3 -flto -fPIC
LFLAGS=-lm -L.
libssdp_la_LDFLAGS = -lm -Wl,--version-script=libssdp.map

all: ssdp
ssdp: libssdp.so util.o io.o 
	$(CC) ssdp.c util.o io.o -o ssdp $(LFLAGS) -lssdp
libssdp.so: Makefile libssdp.map $(OBJ)
	$(CC) -shared -o libssdp.so $(OBJ) $(libssdp_la_LDFLAGS)
libssdp.o: Makefile libssdp.c $(HDR)
project.o:project.h sky_dome.h vector.h util.h
io.o:project.h sky_dome.h vector.h util.h
sky_dome.o:sky_dome.h vector.h util.h
sky_model.o:sky_dome.h util.h
vector.o:sky_dome.h util.h
ground.o:sky_dome.h vector.h util.h
clean:
	rm *.o ssdp
