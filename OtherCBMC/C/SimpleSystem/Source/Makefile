CC          = gcc
CLINKER     = gcc

CFLAGS      = -Wall -O2 -g 
LIBS        = -lm
DEPEND= makedepend

SRC     = cbmc.c ran_uniform.c
OBJS    = cbmc.o ran_uniform.o
EXECS   = cbmc.exe

default: cbmc

all: cbmc

cbmc:$(OBJS)
	$(CLINKER) $(OPTFLAGS) -o $(EXECS) $(OBJS) $(LIBS)

clean:
	/bin/rm -f *.o *~ $(EXECS)

.c.o:
	$(CC) $(CFLAGS) -c $*.c
