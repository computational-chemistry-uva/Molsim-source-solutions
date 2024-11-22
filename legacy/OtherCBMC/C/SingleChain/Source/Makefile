CC          = gcc
CLINKER     = gcc

CFLAGS      = -Wall -O2 -g 
LIBS        = -lm
DEPEND= makedepend

SRC     = grow.c mcloop.c polymer.c ran_uniform.c sample.c system.c writepdb.c
OBJS    = grow.o mcloop.o polymer.o ran_uniform.o sample.o system.o writepdb.o
EXECS   = polymer.exe

default: polymer

all: polymer

polymer:$(OBJS)
	$(CLINKER) $(OPTFLAGS) -o $(EXECS) $(OBJS) $(LIBS)

clean:
	/bin/rm -f *.o *~ $(EXECS)

.c.o:
	$(CC) $(CFLAGS) -c $*.c

angle.o: system.h
grow.o: system.h 
mcloop.o: system.h 
polymer.o: system.h 
ran_uniform.o: system.h 
repulsion.o: system.h 
sample.o: system.h 
system.o: system.h 
writepdb.o: system.h
