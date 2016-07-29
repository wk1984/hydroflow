FFLAGS = -mcmodel=medium -lnetcdff -I/usr/local/include
FC = gfortran
MAIN = driver

all: exec

debug: FFLAGS += -g -fbounds-check -Wall -Wtabs -fbacktrace -finit-real=nan
debug: exec

SOURCES = $(wildcard *_mod.f03)
OBJECTS = $(SOURCES:.f03=.o)


exec: $(MAIN).f03 $(OBJECTS)
	$(FC) $(FFLAGS) -o $(MAIN) $^

ctrl_mod.o: ncutils_mod.o
files_mod.o: ncutils_mod.o ctrl_mod.o dem_mod.o

%.o: %.f03
	$(FC) $(FFLAGS) -c $< -o $@


clean:
	rm $(MAIN) *.o *.mod