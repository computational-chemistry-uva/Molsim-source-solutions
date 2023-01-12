FCOMPFLAGS    =	-O2 -Wall -g
CCOMPFLAGS    =	-O2 -Wall -g

FFLAGS        =	$(FCOMPFLAGS)
CFLAGS        = $(CCOMPFLAGS)
LDFLAGS       =	$(FCOMPFLAGS)

LD            =	gfortran
FC            =	gfortran
CC            =	cc

MAKEFILE      =	Makefile
PROGRAM       =	barrier2
EXECS         =	barrier2.exe

INCS	      =	system.inc

OBJS	      =	barrier2.o          \
		force.o             \
		integrate.o         \
		mdloop.o            \
		ran_uniform.o       \
		ran_vel.o           \
		ran_gauss.o         \
		readdat.o           \
		sample.o            \
		sstmm.o

all:		$(PROGRAM)

$(PROGRAM)::	$(INCS)
		@/bin/rm -f $(OBJS) core

$(PROGRAM)::	$(OBJS) $(MAKEFILE)
		@$(LD) $(LDFLAGS) $(OBJS) -o $(EXECS)

clean:;		@rm -f $(OBJS) core $(EXECS)
