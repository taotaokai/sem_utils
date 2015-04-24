FC = ifort

FFLAGS_INC := -I./include
FFLAGS := $(FFLAGS_INC)
LDFLAGS :=

PROG := xinterp_mesh 
SRCS := xinterp_mesh.f90 \
				get_program_args.f90 \
				read_mesh.f90 \
 				cube2xyz_hex333.f90 \
				xyz2cube_hex333.f90 \
				gll_library.f90 \
				index_hex333_in_spec555.f90 \
				lagrange.f90
				
OBJS := $(patsubst %.f90,%.o, $(SRCS))

all : $(PROG) 

$(PROG) : $(OBJS)
	$(FC) $(OBJS) -o $(PROG)

xinterp_mesh.o : xinterp_mesh.f90
	$(FC) $(FFLAGS) -c xinterp_mesh.f90

get_program_args.o : get_program_args.f90
	$(FC) $(FFLAGS) -c get_program_args.f90

read_mesh.o : read_mesh.f90
	$(FC) $(FFLAGS) -c read_mesh.f90

cube2xyz_hex333.o : cube2xyz_hex333.f90
	$(FC) $(FFLAGS) -c cube2xyz_hex333.f90

xyz2cube_hex333.o : xyz2cube_hex333.f90
	$(FC) $(FFLAGS) -c xyz2cube_hex333.f90

gll_library.o : gll_library.f90
	$(FC) $(FFLAGS) -c gll_library.f90

index_hex333_in_spec555.o : index_hex333_in_spec555.f90
	$(FC) $(FFLAGS) -c index_hex333_in_spec555.f90

lagrange.o : lagrange.f90
	$(FC) $(FFLAGS) -c lagrange.f90

clean :
	/bin/rm $(PROG) $(OBJS) 
