CC=gfortran
CFLAGS=-O3

LIBS=-lnetcdff
DEPS = datetime_module.f90 flip_to_cf.f90
OBJ = datetime_module.o flip_to_cf.o

%.o: %.f90 $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

mpas_to_latlon: $(OBJ)
	$(CC) -o ../bin/$@ $^ $(CFLAGS) $(LIBS) mpas_to_latlon.f90

clean:
	rm *.o *.mod
