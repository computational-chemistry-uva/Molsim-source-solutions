CC          = gcc
CLINKER     = gcc

CFLAGS      = -Wall -O2 -g
LIBS        = -lm
DEPEND= makedepend

SRC     = system.c ewald.c fourierspace.c realspace.c lattice.c ran_uniform.c
OBJS    = system.o ewald.o fourierspace.o realspace.o lattice.o ran_uniform.o
EXECS   = ewald.exe

default: ewald

all: ewald

ewald:$(OBJS)
	$(CLINKER) $(OPTFLAGS) -o $(EXECS) $(OBJS) $(LIBS)

clean:
	/bin/rm -f *.o *~ $(EXECS)

.c.o:
	$(CC) $(CFLAGS) -c $*.c

ewald.o: system.h
fourierspace.o: system.h  
lattice.o: system.h  
realspace.o: system.h  
system.o: system.h
