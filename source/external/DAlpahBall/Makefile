
NAME1=DAlphaBall
DIR=.

# # FOR INTEL COMPILER
# #
# CPU=icc
# NAMEFUL1=$(DIR)/$(NAME1).$(CPU)
# FC = /opt/intel/fce/10.1.008/bin/ifort
# FFLAGS = -c -u -O -cpp
# CC = /opt/intel/cce/10.1.013/bin/icc
# LD = /opt/intel/cce/10.1.013/bin/icc
# CFLAGS = -c -O -Dintel
# LDFLAGS = -O -static
# LIBRARIES= -lgmp -lm -lifcore
# LIBS=-L/opt/intel/fce/10.1.008/lib/



# # FOR GCC COMPILER LINUX 
# # might need to repalce gcc4 wtih gcc
# #
# CPU=gcc
# NAMEFUL1=$(DIR)/$(NAME1).$(CPU)
# FC = gfortran
# FFLAGS = -c -u -O3 -x f77-cpp-input -fsecond-underscore
# CC = gcc4
# LD = gcc4
# CFLAGS = -c -O3
# LDFLAGS = -O3 #-static
# LIBRARIES= -lgmp  -lm -lgfortran



# # FOR GCC COMPILER MAC
# #
# CPU=macgccdebug
# NAMEFUL1=$(DIR)/$(NAME1).$(CPU)
# FC = gfortran
# FFLAGS = -c -O0 -ggdb -x f77-cpp-input -fsecond-underscore
# CC = gcc
# LD = gcc
# CFLAGS = -c -O0 -ggdb
# LDFLAGS = -O0 -ggdb
# LIBRARIES= -lgmp  -lm -lgfortran
# INCS=
# LIBS=

CPU=gcc
NAMEFUL1=$(DIR)/$(NAME1).$(CPU)
FC = gfortran
FFLAGS = -c -O3 -x f77-cpp-input -fsecond-underscore
CC = gcc
LD = gcc
CFLAGS = -c -O3
LDFLAGS = -O3
LIBRARIES= -lgmp  -lm -lgfortran
INCS=-I/usr/local/include
LIBS=-L/usr/local/lib


.c.o :
	$(CC) $(CFLAGS) $(INCS) $<

.f.o :
	$(FC) $(FFLAGS) $<

OBJECTS1 = \
$(NAME1).o \
sos_minor_gmp.o alf_tools_gmp.o binding.o \
adjust.o alfcx.o alfcx_tools.o cavballs.o delcx.o truncate_real.o \
measure_tools.o surface_tools.o surface_only.o surface_deriv_coord.o \
volume_tools.o volume_only.o volume_deriv_coord.o \
vector.o write_simplices.o

$(NAMEFUL1) : $(OBJECTS1)
	$(LD) -o $(NAMEFUL1) $(LDFLAGS) $(OBJECTS1) $(LIBS) $(LIBRARIES)

all: $(OBJECTS1) 
	$(LD) -o $(NAMEFUL1) $(LDFLAGS) $(OBJECTS1) $(LIBS) $(LIBRARIES)

clean:
	rm -f *.o $(NAMEFUL1)

$(OBJECTS) : defines.h gmpvar.h

# gcc4 -static -o ../bin/DAlphaBall.gccstatic -O3  DAlphaBall.o sos_minor_gmp.o alf_tools_gmp.o binding.o adjust.o alfcx.o alfcx_tools.o cavballs.o delcx.o truncate_real.o measure_tools.o surface_tools.o surface_only.o surface_deriv_coord.o volume_tools.o volume_only.o volume_deriv_coord.o vector.o write_simplices.o /usr/lib/gcc/x86_64-redhat-linux/4.1.0/libgfortran.a -lgmp  -lm 
