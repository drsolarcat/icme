# BEGIN OF MAKEFILE
# ===============================================
# The compiler
# CF = ifort
# CF = i686-pc-linux-gnu-g95
# CF = f95
CF = gfortran
CF90 = $(CF)

%.o: %.mod
# ===============================================
# Flags for debuggin or for maximum performance
# debugging for "ifort"
#FFLAGS = -check all -free -fpic -ftz -g -warn all -fpe0 -traceback
# debugging for "g95" or "gfortran"
#FFLAGS = -Wall -ffpe-trap=list -fbacktrace -fdump-core
# performance for all
#FFLAGS = -O2
PRECISION = -fdefault-real-8 -fdefault-integer-8
CF90FLAGS = $(PRECISION) $(FFLAGS) -g

# Limit of memory usage
RLIMIT_STACK = 262145

# ===============================================
# list of object files
OBJECTS = TOMS642.o LM.o \
	types.o constants.o extblas.o resample.o datetime.o \
	config.o data.o mva.o  gsr.o cm.o

# ===============================================
# list of linked libraries
GSL_LIB_PATH = -L/home/isavnin/opt/lib64
#GSL_LIB_PATH = -L/usr/lib64
GSL_INCLUDE_PATH = -I/home/isavnin/opt/include/$(CF90)
#GSL_INCLUDE_PATH = -I/usr/include
GSL = $(GSL_LIB_PATH) $(GSL_INCLUDE_PATH) -lgsl -lgslcblas -lfgsl_$(CF90)
#GSL = $(GSL_LIB_PATH) $(GSL_INCLUDE_PATH) -lgsl -lgslcblas -lfgsl

# ===============================================
all: icme

# 3rd party programs
TOMS642.o: TOMS642.F90
	$(CF90) -c $(CF90FLAGS) TOMS642.F90
LM.o: LM.F90
	$(CF90) -c $(CF90FLAGS) LM.F90

# secondary modules
types.o: types.f90
	$(CF90) -c $(CF90FLAGS) types.f90
constants.o: constants.f90
	$(CF90) -c $(CF90FLAGS) constants.f90
extblas.o: extblas.f90
	$(CF90) -c $(CF90FLAGS) extblas.f90
resample.o: resample.f90
	$(CF90) -c $(CF90FLAGS) resample.f90
datetime.o: datetime.f90
	$(CF90) -c $(CF90FLAGS) datetime.f90

# primary modules
config.o: types.o config.f90
	$(CF90) -c $(CF90FLAGS) config.f90
data.o: datetime.o types.o data.f90
	$(CF90) -c $(CF90FLAGS) data.f90 $(GSL)
mva.o: mva.f90 types.o extblas.o
	$(CF90) -c $(CF90FLAGS) mva.f90 $(GSL)
gsr.o: extblas.o constants.o types.o mva.o gsr.f90
	$(CF90) -c $(CF90FLAGS) gsr.f90 $(GSL)
cm.o: types.o constants.o cm.f90
	$(CF90) -c $(CF90FLAGS) cm.f90 $(GSL)

# executables
icme: $(OBJECTS)
	$(CF90) -o icme $(CF90FLAGS) icme.f90 $(OBJECTS) $(GSL)
data_read: data.o types.o data_read.f90
	$(CF90) -o data_read $(CF90FLAGS) types.o data.o datetime.o data_read.f90

# cleaning everything
clean:
	rm $(OBJECTS)
	rm icme
	rm data_read
# END OF MAKEFILE

