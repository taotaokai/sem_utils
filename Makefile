#FC = ifort
#FC = gfortran 
FC = mpif90

objdir = obj
bindir = bin
incdir = include
srcdir = src
module_dir = ${srcdir}/module
shared_dir = ${srcdir}/shared

FCFLAGS = -g -fbounds-check -Wall -pedantic # gfortran
#FCFLAGS = -O2
#FCFLAGS += $(FFLAGS_INC) -mod $(ODIR) # ifort
#FCFLAGS += $(FFLAGS_INC) -O3 -module $(ODIR) -assume byterecl # mpif90 
FCFLAGS += -I$(incdir) -J$(objdir)
LDFLAGS =

program =
object =
module_ = sem_constants_mod sem_io_mod sem_mesh_mod \
		  sem_parallel_mod sem_utils_mod
shared_ = gll_library

#------------------------------------------

module_obj = $(patsubst %,$(objdir)/%.o, $(module_))
shared_obj = $(patsubst %,$(objdir)/%.o, $(shared_))

all : $(shared_obj) $(module_obj)

$(shared_obj) : 
	$(FC) -c $(shared_dir)/$(patsubst %.o,%.f90,$(@F)) -o $@ $(FCFLAGS)

$(module_obj) : $(shared_obj)
	$(FC) -c $(module_dir)/$(patsubst %.o,%.f90,$(@F)) -o $@ $(FCFLAGS)

#$(module_obj) : $(shared_obj)
#	$(FC) -c $< -o $@ $(FCFLAGS)

# explicit specified dependencies
#$(module): $(shared)
#$(OBJ): $(SHARED) $(MOD)
#$(SHARED) : $(ODIR)/constants_module.o

.PHONY: clean

clean :
	\rm $(bindir)/* $(objdir)/*
