SRC=ssdp.c sky_dome.c sky_model.c project.c vector.c ground.c io.c util.c
OBJ=ssdp.o sky_dome.o sky_model.o project.o vector.o ground.o io.o util.o
HDR=       sky_dome.h sky_model.h project.h vector.h ground.h io.h util.h
CC=gcc
# CFLAGS=-Og -flto -g -Wall -pedantic
CFLAGS=-O3 -flto
LFLAGS=-lm

all: ssdp
ssdp: Makefile $(SRC) $(HDR) $(OBJ)
	$(CC) -o ssdp $(OBJ) $(LFLAGS)
project.o:project.h sky_dome.h vector.h util.h
io.o:project.h sky_dome.h vector.h util.h
sky_dome.o:sky_dome.h vector.h util.h
sky_model.o:sky_dome.h util.h
vector.o:sky_dome.h util.h
ground.o:sky_dome.h vector.h util.h
clean:
	rm *.o ssdp
