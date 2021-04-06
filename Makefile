# compiler
FC = gfortran

# compile flags
FCFLAGS = -g -c -fbacktrace -fbackslash -fno-align-commons -fbounds-check -std=legacy -I./SOFTWARE
# link flags
FLFLAGS =-llapack -lblas -L ./SOFTWARE -lSparseBLAS_GNU #/home/javi/.local/lapack/libblas.a /home/javi/.local/lapack/liblapack.a
# FLFLAGS =  -Wl,--start-group /opt/intel/oneapi/mkl/2021.2.0/lib/intel64/libmkl_gf_lp64.a /opt/intel/oneapi/mkl/2021.2.0/lib/intel64/libmkl_sequential.a /opt/intel/oneapi/mkl/2021.2.0/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl#-I${F95ROOT}/include/intel64/ilp64 -fdefault-integer-8  -m64  -I"${MKLROOT}/include"#-I/usr/lib/x86_64-linux-gnu/lapack -lblas -L/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -llapack

# source files and objects
SRCS = \
ogpf.o \
progress_mod.o \
number_types_mod.o \
param_mod.o \
point_mod.o \
quadrilateral_mod.o \
matrix_mod.o \
out_mod.o \
fem.o

# program name
PROGRAM = blah

all: $(PROGRAM)

$(PROGRAM): $(SRCS)
	$(FC) -o $@ $^ $(FLFLAGS)

%.o: %.f90
	$(FC) $(FCFLAGS)  -o $@ $<

clean:
	rm -f *.o *.mod