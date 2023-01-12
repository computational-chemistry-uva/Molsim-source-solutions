FCOMPFLAGS    =	-O2 -Wall -g -ffixed-line-length-132
CCOMPFLAGS    =	-O2 -Wall -g

FFLAGS        =	$(FCOMPFLAGS)
CFLAGS        = $(CCOMPFLAGS)
LDFLAGS       =	$(FCOMPFLAGS)

LD            =	gfortran
FC            =	gfortran
CC            =	cc

MAKEFILE      =	Makefile
PROGRAM       =	mc_nvt
EXECS         =	mc_nvt.exe

INCS	      =	system.inc conf.inc parameter.inc potential.inc

OBJS	      =	adjust.o          \
		energyparticle.o  \
		lattice.o         \
		mc_nvt.o          \
		mcmove.o          \
		readdat.o         \
		sample.o          \
		sstmm.o           \
		store.o           \
		energysystem.o    \
		writepdb.o        \
		ran_uniform.o

all:		$(PROGRAM)

$(PROGRAM)::	$(INCS)
		@/bin/rm -f $(OBJS) core

$(PROGRAM)::	$(OBJS) $(MAKEFILE)
		@$(LD) $(LDFLAGS) $(OBJS) -o $(EXECS)

clean:;		@rm -f $(OBJS) core $(EXECS)
