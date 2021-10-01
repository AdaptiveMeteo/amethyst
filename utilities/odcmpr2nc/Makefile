odcmpr2nc: odcmpr2nc.f90 
        #gfortran -02 -g -o asolar_source_function.o asolar_source_function.f90
	gfortran -O2 -g -o odcmpr2nc odcmpr2nc.f90 -I/usr/include  -lnetcdf -lnetcdff -ffast-math -ffree-form -fno-second-underscore -frecord-marker=4 -funroll-loops -ggdb
#       with conda
#	gfortran -O2 -g  -o odcmpr2nc odcmpr2nc.f90 -I${CONDA_PREFIX}/include  -L${CONDA_PREFIX}/lib -Wl,-rpath=${CONDA_PREFIX}/lib -lnetcdf -lnetcdff

clean:
	rm -f odcmpr2nc

