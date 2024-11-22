FCOMPFLAGS    =	-O2 -Wall -g -ffixed-line-length-132
CCOMPFLAGS    =	-O2 -Wall -g

FFLAGS        =	$(FCOMPFLAGS)
CFLAGS        = $(CCOMPFLAGS)
LDFLAGS       =	$(FCOMPFLAGS)

LD            =	gfortran
FC            =	gfortran
CC            =	cc

MAKEFILE      =	Makefile
PROGRAM       =	Gibbs
EXECS         =	Gibbs.exe

INCS	      =	chem.inc      \
		conf.inc      \
		parameter.inc \
		potential.inc \
		system.inc

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
		sstmm.o       \
		toterg.o

all:		$(PROGRAM)

$(PROGRAM)::	$(INCS)
		@/bin/rm -f $(OBJS) core

$(PROGRAM)::	$(OBJS) $(MAKEFILE)
		@$(LD) $(LDFLAGS) $(OBJS) -o $(EXECS)

clean:;		@rm -f $(OBJS) core $(EXECS)
