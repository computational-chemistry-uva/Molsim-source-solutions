CC          = gcc
CLINKER     = gcc

CFLAGS      = -Wall -O2 -g
LIBS        = -lm
DEPEND= makedepend

SRC        = boltzmann.c
OBJS       = boltzmann.o
EXECS      = boltzmann.exe

default: boltzmann

all: boltzmann

boltzmann:$(OBJS)
	$(CLINKER) $(OPTFLAGS) -o $(EXECS) $(OBJS) $(LIBS)

clean:
	/bin/rm -f *.o *~ $(EXECS)

.c.o:
	$(CC) $(CFLAGS) -c $*.c
