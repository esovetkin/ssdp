libSRC=sky_dome.c sky_model.c pv-aoi.c project.c trace.c vector.c ground.c print.c libssdp.c shull.c delaunay.c ll.c error.c sunpos.c
libOBJ=sky_dome.o sky_model.o pv-aoi.o project.o trace.o vector.o ground.o print.o libssdp.o shull.o delaunay.o ll.o error.o sunpos.o
libHDR=sky_dome.h sky_model.h pv-aoi.h project.h trace.h vector.h ground.h print.h libssdp.h shull.h delaunay.h ll.h error.h sunpos.h

SRC=io.c variables.c parser.c readlineshell.c util.c ssdp.c
OBJ=io.o variables.o parser.o readlineshell.o util.o
HDR=io.h variables.h parser.h readlineshell.h util.h 
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

all: libssdp.so ssdp
io.o:project.h sky_dome.h vector.h print.h
libssdp.so: Makefile libssdp.map $(libOBJ)
	$(CC) -shared -o libssdp.so $(libOBJ) $(libssdp_la_LDFLAGS)
libssdp.o: Makefile libssdp.c $(libHDR) error.h errorflags.h libssdp_structs.h
libssdp_structs.h: collect_datastructs.sh vector.h sky_dome.h delaunay.h trace.h ground.h project.h libssdp.c
	${SHELL} collect_datastructs.sh vector.h sky_dome.h delaunay.h trace.h ground.h project.h libssdp.c
pv-aoi.o:pv-aoi.c
project.o:pv-aoi.h project.h sky_dome.h vector.h print.h error.h errorflags.h project.c
sky_dome.o:sky_dome.h vector.h print.h error.h errorflags.h sky_dome.c
sky_model.o:sky_dome.h print.h error.h errorflags.h sky_model.c
vector.o:sky_dome.h print.h error.h errorflags.h vector.c
ground.o:sky_dome.h vector.h print.h errorflags.h ground.c
delaunay.o:ll.h shull.h error.h errorflags.h delaunay.c
shull.o:ll.h shull.h error.h errorflags.h shull.c 
	$(CC) -c -DLL_CIRCULAR -D_GNU_SOURCE $(CFLAGS) -o shull.o shull.c
ll.o:ll.h error.h ll.c
	$(CC) -c -DLL_CIRCULAR -D_GNU_SOURCE $(CFLAGS) -o ll.o ll.c
error.o: error.h errorflags.h errormessages.h
errorflags.h errormessages.h: gen_errorflags.sh $(libSRC)
	${SHELL} gen_errorflags.sh  $(libSRC)
ssdp: libssdp.so $(OBJ)
	$(CC) ssdp.c $(OBJ) -o ssdp $(LFLAGS) -lssdp -lreadline
io.o:project.h sky_dome.h vector.h print.h util.h 
parser.o: parsedef.h parser.h variables.h util.h 
	$(CC) -c -D_GNU_SOURCE $(CFLAGS) -o parser.o parser.c
readlineshell.o: parser.h
parsedef.h: gen_parserflags.sh parser.c
	${SHELL} gen_parserflags.sh parser.c
clean:
	rm *.o ssdp libssdp.so errorflags.h errormessages.h parsedef.h
