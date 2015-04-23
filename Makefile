#FC = ifort
#FC = gfortran 
FC = mpif90

ODIR = obj
SDIR = src
IDIR = inc
EDIR = bin

FFLAGS_INC := -I$(IDIR)
FCFLAGS := $(FFLAGS_INC) -Wall -pedantic -J$(ODIR) # gfortran
#FCFLAGS := $(FFLAGS_INC) -mod $(ODIR) # ifort
#FCFLAGS := $(FFLAGS_INC) -O3 -module $(ODIR) -assume byterecl # mpif90 
LDFLAGS :=

PROG_ := xsum_kernels_cijkl_to_iso \
				 xmask_kernel \
				 xadd_dmodel_iso_to_tiso \
				 xsem_binary_op \
 				 xget_dmodel_iso_lbfgs \
				 xinfo_max_dlnv_tiso_plus_iso
OBJ_ := sum_kernels_cijkl_to_iso.o \
				mask_kernel.o \
				add_dmodel_iso_to_tiso.o \
				sem_binary_op.o \
				get_dmodel_iso_lbfgs.o \
				info_max_dlnv_tiso_plus_iso.o
MOD_ := constants_module.o sem_IO_module.o sem_tomography_module.o parallel_module.o
SHARED_ := gll_library.o

#------------------------------------------
PROG = $(patsubst %,$(EDIR)/%,$(PROG_))
OBJ = $(patsubst %,$(ODIR)/%,$(OBJ_))
MOD = $(patsubst %,$(ODIR)/%,$(MOD_))
SHARED = $(patsubst %,$(ODIR)/%,$(SHARED_))

all : $(PROG)

$(EDIR)/xsum_kernels_cijkl_to_iso: $(SHARED) $(MOD) $(ODIR)/sum_kernels_cijkl_to_iso.o
	$(FC) -o $@ $^ $(FCFLAGS) $(LDFLAGS)

$(EDIR)/xmask_kernel: $(SHARED) $(MOD) $(ODIR)/mask_kernel.o
	$(FC) -o $@ $^ $(FCFLAGS) $(LDFLAGS)

$(EDIR)/xadd_dmodel_iso_to_tiso: $(SHARED) $(MOD) $(ODIR)/add_dmodel_iso_to_tiso.o
	$(FC) -o $@ $^ $(FCFLAGS) $(LDFLAGS)

$(EDIR)/xget_dmodel_iso_lbfgs: $(SHARED) $(MOD) $(ODIR)/get_dmodel_iso_lbfgs.o
	$(FC) -o $@ $^ $(FCFLAGS) $(LDFLAGS)

$(EDIR)/xsem_binary_op: $(SHARED) $(MOD) $(ODIR)/sem_binary_op.o
	$(FC) -o $@ $^ $(FCFLAGS) $(LDFLAGS)

$(EDIR)/xinfo_max_dlnv_tiso_plus_iso: $(SHARED) $(MOD) $(ODIR)/info_max_dlnv_tiso_plus_iso.o
	$(FC) -o $@ $^ $(FCFLAGS) $(LDFLAGS)

$(ODIR)/%.o: $(SDIR)/%.f90
	$(FC) -c $< -o $@ $(FCFLAGS)

# explicit specified dependencies
$(MOD): $(SHARED)
$(OBJ): $(SHARED) $(MOD)
#$(SHARED) : $(ODIR)/constants_module.o

.PHONY: clean

clean :
	\rm $(EDIR)/* $(ODIR)/*
