#FC = ifort
FC = gfortran 

ODIR = obj
SDIR = src
IDIR = inc
EDIR = bin

FFLAGS_INC := -I$(IDIR)
FCFLAGS := $(FFLAGS_INC) -Wall -pedantic -J$(ODIR) # gfortran
#FCFLAGS := $(FFLAGS_INC) -mod $(ODIR) # ifort
#FCFLAGS := $(FFLAGS_INC) -O3 -module $(ODIR) -assume byterecl # mpif90 
LDFLAGS :=

PROG := $(EDIR)/interp_xyz
OBJ_ := interp_xyz.o
MOD_ := constants_module.o mesh_module.o

#------------------------------------------

OBJ = $(patsubst %,$(ODIR)/%,$(OBJ_))
MOD = $(patsubst %,$(ODIR)/%,$(MOD_))

all : $(PROG)

$(PROG): $(MOD) $(OBJ) 
	$(FC) -o $@ $^ $(FCFLAGS) $(LDFLAGS)

$(ODIR)/%.o: $(SDIR)/%.f90
	$(FC) -c $< -o $@ $(FCFLAGS)

# explicit specified dependencies
$(OBJ): $(MOD)
$(ODIR)/mesh_module.o: $(ODIR)/constants_module.o

.PHONY: clean

clean :
	\rm $(PROG) $(ODIR)/*
