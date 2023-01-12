FCOMPFLAGS    =	-O2 -Wall -g -ffixed-line-length-132
CCOMPFLAGS    =	-O2 -Wall -g

FFLAGS        =	$(FCOMPFLAGS)
CFLAGS        = $(CCOMPFLAGS)
LDFLAGS       =	$(FCOMPFLAGS)

LD            =	gfortran
FC            =	gfortran
CC            = cc

MAKEFILE      =	Makefile
PROG1         =	mc-npt
EXECS         =	mc-npt.exe

INCS	      =	system.inc

OBJ1	      =	mc-npt.o writepdb.o ran_uniform.o sstmm.o

all:		$(PROG1)

$(PROG1)::	$(INCS)
		@/bin/rm -f $(OBJS) core

$(PROG1)::	$(OBJ1) $(OBJS) $(MAKEFILE)
		@$(LD) $(LDFLAGS) $(OBJ1) -o $(EXECS)

clean:;		@rm -f $(OBJ1) core $(EXECS)
