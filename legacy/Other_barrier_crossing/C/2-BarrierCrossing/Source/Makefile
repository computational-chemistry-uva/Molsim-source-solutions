CC          = gcc
CLINKER     = gcc

CFLAGS      = -Wall -O2 -g 
LIBS        = -lm
DEPEND= makedepend

SRC     = barrier2.c force.c integrate.c mdloop.c ran_uniform.c readdat.c sample.c system.c
OBJS    = barrier2.o force.o integrate.o mdloop.o ran_uniform.o readdat.o sample.o system.o
EXECS   = barrier2.exe

default: barrier2

all: barriers

barrier2:$(OBJS)
	$(CLINKER) $(OPTFLAGS) -o $(EXECS) $(OBJS) $(LIBS)

clean:
	/bin/rm -f *.o *~ $(EXECS)

.c.o:
	$(CC) $(CFLAGS) -c $*.c

barrier2.o: system.h
force.o: system.h
integrate.o: system.h 
mdloop.o: system.h 
ran_uniform.o: system.h 
readdat.o: system.h 
sample.o: system.h 
system.o: system.h
