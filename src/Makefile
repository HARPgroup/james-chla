#

#FC = pgf90
FC = gfortran
LN = $(F95)


#

# options for different environments
# (compiler, linker, compilation flags, libraries, etc)
#

## set environment 
##  pgi compiler
## ifeq ($(XCOMP),ifort)
##  LFLAGS = -Bstatic
  FPP_FLAGS = -Mpreprocess

##  ifeq ($(XARCH),amd64b)
##      FFLAGS = -fastsse -O3 -tp piledriver  -m64  -mcmodel=medium
##  else
      FFLAGS = -fastsse -O2 -Mipa=fast -tp amd64  -mcmodel=medium
##  endif
 
##  LIBS =

##endif

#
# set preprocessing options (and combine with FFLAGS)
#


FFLAGS := $(FPP_FLAGS) $(FFLAGS)

#
# comment out the following line if you don't want the build dependent
# upon the makefile, otherwise use the name of the makefile
#
MAKEFILE = Makefile

# the files
EXEC=jameswq
#
# the actual build commands
#
#.SUFFIXES: .o .for

default: $(EXEC)

$(EXEC): *.for $(MAKEFILE)
	$(FC)  $(FFLAGS) $(LFLAGS) *.for $(LIBS) -o $(EXEC)

#
# how to clean _all_ up
#
clean:
	rm *.o *.mod $(EXEC)

