FC = ifort
#FC = mpif90
ODIR := obj

FFLAGS_INC := -I./
#FCFLAGS := $(FFLAGS_INC) -O3 -Wall -pedantic -J$(ODIR) # gfortran
FCFLAGS := $(FFLAGS_INC) -O3 -module $(ODIR) # ifort
#FCFLAGS := $(FFLAGS_INC) -O3 -module $(ODIR) -assume byterecl -DFORCE_VERCTORIZATION  # mpif90 
LDFLAGS :=

PROG := interp_xyz
OBJ_ := main.o \
        shared_par.o \
        hex_nodes.o \
        anchor_point_index.o \
        read_mesh.o \
        hpsort_index.o \
        cube2xyz.o \
        xyz2cube_bounded.o \
				lagrange_poly.o

OBJ = $(patsubst %,$(ODIR)/%,$(OBJ_))

all : $(PROG)

# compile modules first
$(ODIR)/%.o: %.F90
	$(FC) -c $< -o $@ $(FCFLAGS)

# compile subroutines and main
$(ODIR)/%.o: %.f90 $(ODIR)/shared_par.o
	$(FC) -c $< -o $@ $(FCFLAGS)

# link
$(PROG): $(OBJ)
	$(FC) -o $@ $^ $(FCFLAGS) $(LDFLAGS)

.PHONY: clean

clean :
	\rm $(PROG) $(ODIR)/*
