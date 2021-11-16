FORT = gfortran
#NICE = -r8 -zero
NICE = -fdefault-real-8 
# LIB = /home/bryngel2/packages/fftw/fftw-3.3.4/lib/libfftw3.a /home/bryngel2/packages/BLAS/blas_LINUX.a /home/bryngel2/packages/lapack-3.5.0/liblapack.a /home/bryngel2/packages/BLAS/blas_LINUX.a 
LIB =  /usr/local/Cellar/fftw/3.3.8_2/lib/libfftw3.a -L/usr/local/opt/lapack/lib -llapack -lblas -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib

# /usr/local/Cellar/lapack/3.9.0_1/lib/liblapack.3.9.0.dylib /usr/local/Cellar/lapack/3.9.0_1/lib/libblas.3.9.0.dylib

#LIB = /usr/local/lib/libfftw3.a
#LIB = /export/home/das7/libfftw3.a /home/bryngel2/packages/BLAS/blas_LINUX.a /home/bryngel2/packages/lapack-3.5.0/liblapack.a /home/bryngel2/packages/BLAS/blas_LINUX.a
#LIB  =  /users/spencerbryngelson/downloads/lapack/liblapack.a /users/spencerbryngelson/downloads/lapack/librefblas.a /usr/local/lib/libfftw3.a
#LIB  = -lfftw3 -L/home/jfreund/Code/FS/V1/lib -lblas -llapack -lblas
#LIB  = -lfftw3 -L/home/jfreund/lib/LAPACK -lblas -llapack -lblas
OPTS = -O2
## Note: Does not compile with -O0 gfortran for some reason... 
#PROF = -pg
DEBUG = -g -C -ffixed-line-length-none -fallow-argument-mismatch 
#-Wall -fcheck=all -g -fbacktrace 


# DEBUG = -g -C -ffixed-line-length-none -fallow-argument-mismatch -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -std=legacy  -pedantic  -fbacktrace

OBJECTS = prms.o data.o rk.o main.o ic.o membrane.o stokes.o out.o surffor.o \
              dPackgmres.o ei.o area.o nlist.o pme.o qflow.o

%.o: %.f
	$(FORT) -c $(NICE) $(OPTS) $(PROF) $(DEBUG) $<

%.o: %.f90
	$(FORT) -c $(NICE) $(OPTS) $(PROF) $(DEBUG) $<

all: stokes

stokes:   $(OBJECTS) makefile
	$(FORT) -o stokes $(PROF) $(OBJECTS) $(LIB)

main.o: main.f90 prms.o data.o rk.o membrane.o out.o ic.o area.o ei.o pme.o stokes.o qflow.o makefile

out.o:  out.f90 prms.o data.o membrane.o stokes.o makefile

stokes.o: stokes.f90 surffor.o prms.o data.o pme.o membrane.o nlist.o dPackgmres.o ei.o qflow.o makefile

ic.o:   ic.f90 prms.o data.o makefile

membrane.o: membrane.f90 prms.o data.o makefile

area.o: area.f90 prms.o data.o membrane.o makefile

rk.o: rk.f90 prms.o data.o stokes.o makefile

surffor.o: surffor.f90 prms.o data.o nlist.o membrane.o makefile

pme.o: pme.f90 prms.o data.o nlist.o makefile

qflow.o: qflow.f90 prms.o makefile

data.o: data.f90 prms.o makefile

nlist.o: nlist.f90 prms.o makefile

prms.o: prms.f90 makefile

ei.o:  ei.f90 makefile

dPackgmres.o: dPackgmres.f makefile

clean:
	rm -f flow *.mod $(OBJECTS)
