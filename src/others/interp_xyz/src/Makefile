FC = ifort
#FC = gfortran 

ODIR = obj
SDIR = src
IDIR = inc
EDIR = bin

FFLAGS_INC := -I$(IDIR)
#FCFLAGS := $(FFLAGS_INC) -Wall -pedantic -J$(ODIR) # gfortran
FCFLAGS := $(FFLAGS_INC) -mod $(ODIR) # ifort
#FCFLAGS := $(FFLAGS_INC) -O3 -module $(ODIR) -assume byterecl # mpif90 
LDFLAGS :=

PROG := $(EDIR)/interp_xyz
OBJ_ := interp_xyz.o
MOD_ := constants_module.o mesh_module.o
SHARED_ :=  anchor_point_index.o \
            cube2xyz.o           \
            hex_nodes.o          \
            hpsort_index.o       \
            lagrange_poly.o      \
            strtok.o             \
            xyz2cube_bounded.o


#------------------------------------------

OBJ = $(patsubst %,$(ODIR)/%,$(OBJ_))
MOD = $(patsubst %,$(ODIR)/%,$(MOD_))
SHARED = $(patsubst %,$(ODIR)/%,$(SHARED_))

all : $(PROG)

$(PROG): $(SHARED) $(MOD) $(OBJ)
	$(FC) -o $@ $^ $(FCFLAGS) $(LDFLAGS)

$(ODIR)/%.o: $(SDIR)/%.f90
	$(FC) -c $< -o $@ $(FCFLAGS)

# explicit specified dependencies
$(OBJ): $(SHARED) $(MOD)
$(ODIR)/mesh_module.o: $(SHARED) $(ODIR)/constants_module.o
$(SHARED) : $(ODIR)/constants_module.o

.PHONY: clean

clean :
	\rm $(PROG) $(ODIR)/*
