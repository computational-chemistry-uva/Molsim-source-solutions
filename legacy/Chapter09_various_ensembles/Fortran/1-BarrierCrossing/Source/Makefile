FCOMPFLAGS    =	-O2 -Wall -g -ffixed-line-length-132
CCOMPFLAGS    =	-O2 -Wall -g

FFLAGS        =	$(FCOMPFLAGS)
CFLAGS        = $(CCOMPFLAGS)
LDFLAGS       =	$(FCOMPFLAGS)

LD            =	gfortran
FC            =	gfortran
CC            =	cc

MAKEFILE      =	Makefile
PROGRAM       =	barrier1
EXECS         =	barrier1.exe

INCS	      =	system.inc

OBJS	      =	force.o             \
		integrate_and.o     \
		integrate_mc.o      \
		integrate_nve.o     \
		integrate_res.o     \
		integrate_ber.o     \
		barrier1.o          \
		mdloop.o            \
		ran_uniform.o       \
		ran_vel.o           \
		ran_gauss.o         \
		readdat.o           \
		sample_diff.o       \
		sample_prof.o       \
		sstmm.o

all:		$(PROGRAM)

$(PROGRAM)::	$(INCS)
		@/bin/rm -f $(OBJS) core

$(PROGRAM)::	$(OBJS) $(MAKEFILE)
		@$(LD) $(LDFLAGS) $(OBJS) -o $(EXECS)

clean:;		@rm -f $(OBJS) core $(EXECS)
