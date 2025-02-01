CC   = mpicc
GCC  = mpicc
LINKER = $(CC)

ifeq ($(ENABLE_OPENMP),true)
OPENMP   = -fopenmp
endif

VERSION  = --version
CFLAGS   = -Ofast -ffreestanding -std=c99 $(OPENMP)
LFLAGS   = $(OPENMP)
DEFINES  = -D_GNU_SOURCE
INCLUDES =
LIBS     = -lm
