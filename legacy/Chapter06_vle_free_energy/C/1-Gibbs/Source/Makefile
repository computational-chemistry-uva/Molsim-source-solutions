CCOMPFLAGS    =	-Wall -g -O2

CFLAGS        = $(CCOMPFLAGS)
LIBS          = -lm

LD            =	gcc
CC            =	gcc

MAKEFILE      =	Makefile
PROGRAM       =	Gibbs
EXECS         =	Gibbs.exe

INCS	      =	chem.h        \
		conf.h        \
		parameter.h   \
		potential.h   \
		system.h system_2.h 

OBJS	      =	adjust.o      \
		ener.o        \
		eneri.o       \
		gibbs.o       \
		init_chem.o   \
		lattice.o     \
		mcmove.o      \
		mcswap.o      \
		mcvol.o       \
		ran_uniform.o \
		readdat.o     \
		sample.o      \
		store.o       \
		toterg.o
#		sstmm.o       \

all:		Gibbs

$(PROGRAM)::	$(INCS)
		@/bin/rm -f $(OBJS) core

$(PROGRAM)::	$(OBJS) $(MAKEFILE)
		@$(LD) $(CFLAGS) $(OBJS) -o $(EXECS) $(LIBS)

clean:;		@rm -f $(OBJS) core $(EXECS)
