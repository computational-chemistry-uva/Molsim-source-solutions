FCOMPFLAGS    =	-O -ffixed-line-length-132
CCOMPFLAGS    =	-O

FFLAGS        =	$(FCOMPFLAGS)
CFLAGS        = $(CCOMPFLAGS)
LDFLAGS       =	$(FCOMPFLAGS)

LD            =	gfortran
FC            =	gfortran
CC            =	cc

MAKEFILE      =	Makefile
PROGRAM       =	cbmc
EXECS         =	cbmc.exe

INCS	      =	

OBJS	      =	cbmc.o              \
		ran_uniform.o       \
		sstmm.o

all:		$(PROGRAM)

$(PROGRAM)::	$(INCS)
		@/bin/rm -f $(OBJS) core

$(PROGRAM)::	$(OBJS) $(MAKEFILE)
		@$(LD) $(LDFLAGS) $(OBJS) -o $(EXECS)

clean:;		@rm -f $(OBJS) core $(EXECS)
