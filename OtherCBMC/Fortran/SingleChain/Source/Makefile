FCOMPFLAGS    =	-O
CCOMPFLAGS    =	-O

FFLAGS        =	$(FCOMPFLAGS)
CFLAGS        = $(CCOMPFLAGS)
LDFLAGS       =	$(FCOMPFLAGS)

LD            =	gfortran
FC            =	gfortran
CC            =	cc

MAKEFILE      =	Makefile
PROGRAM       =	polymer
EXECS         =	polymer.exe

INCS	      =	system.inc

OBJS	      =	angle.o          \
		grow.o           \
		mcloop.o         \
		polymer.o        \
		ran_uniform.o    \
		ran_sphere.o     \
		repulsion.o      \
		sample.o         \
		sstmm.o          \
		writepdb.o

all:		$(PROGRAM)

$(PROGRAM)::	$(INCS)
		@/bin/rm -f $(OBJS) core

$(PROGRAM)::	$(OBJS) $(MAKEFILE)
		@$(LD) $(LDFLAGS) $(OBJS) -o $(EXECS)

clean:;		@rm -f $(OBJS) core $(EXECS)
