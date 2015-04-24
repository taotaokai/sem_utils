FC = gfortran 

IDIR = ../../include
EDIR = ./

FFLAGS_INC := -I$(IDIR)
FCFLAGS := $(FFLAGS_INC) -Wall -pedantic 
LDFLAGS :=

PROG_ := test

all : $(PROG)

.PHONY: clean

clean :
	\rm $(EDIR)/* $(ODIR)/*
