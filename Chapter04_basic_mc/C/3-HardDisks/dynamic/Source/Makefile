CC          = gcc
CLINKER     = gcc

CFLAGS      = -Wall -O2 -g
LIBS        = -lm
DEPEND= makedepend

SRC     = hard2.c sample.c system.c writepdb.c ran_uniform.c
OBJS    = hard2.o sample.o system.o writepdb.o ran_uniform.o
EXECS   = hard2.exe

default: hard2

all: hard2

hard2:$(OBJS)
	$(CLINKER) $(OPTFLAGS) -o $(EXECS) $(OBJS) $(LIBS)

clean:
	/bin/rm -f *.o *~ $(EXECS)

.c.o:
	$(CC) $(CFLAGS) -c $*.c

hard2.o: system.h 
sample.o: system.h
system.o: system.h 
writepdb.o: system.h 
