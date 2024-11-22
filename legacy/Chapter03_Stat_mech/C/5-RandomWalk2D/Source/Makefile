CC          = gcc
CLINKER     = gcc

CFLAGS      = -Wall -O2 -g
LIBS        = -lm
DEPEND= makedepend

SRC      = random2d.c ran_uniform.c sample.c system.c
OBJS     = random2d.o ran_uniform.o sample.o system.o
EXECS    = random2d.exe

default: random2d

all: random2d

random2d:$(OBJS)
	$(CLINKER) $(OPTFLAGS) -o $(EXECS) $(OBJS) $(LIBS)

clean:
	/bin/rm -f *.o *~ $(EXECS)

.c.o:
	$(CC) $(CFLAGS) -c $*.c

ran_uniform.o: system.h
random2d.o: system.h
sample.o: system.h  
system.o: system.h
