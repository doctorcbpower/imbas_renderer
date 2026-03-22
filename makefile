OPT += -DCIC
OPT += -DIMAGE_DIMENSIONX=4096
OPT += -DIMAGE_DIMENSIONY=4096
#OPT += -DLONG_IDS
#OPT += -DVELOCITIES
#OPT += -DIMAGE
#OPT += -DDEBUG
OPT += -DENABLE_OPENMP
#OPT += -DENABLE_MPI
#OPT += -DSCATTER_DECOMPOSITION
OPT += -DNONPERIODIC
#OPT += -DKERNEL_SMOOTHING
#OPT += -DSMOOTH_GRID=512

#COMPILE_ON_SYSTEM="MacBook"
#COMPILE_ON_SYSTEM="MacPro"
#COMPILE_ON_SYSTEM="OzSTAR"
COMPILE_ON_SYSTEM="Setonix"

ifeq ($(COMPILE_ON_SYSTEM),"Setonix")
CC=cc -fopenmp
PNG_HOME=/software/projects/pawsey1164/cpower
PNG_LIBS=$(PNG_HOME)/lib
PNG_INCL=$(PNG_HOME)/include
PNG_OPTS=-lpng
HDF5_INCL=${HDF5_DIR}/include -D H5_USE_16_API
HDF5_LIBS=${HDF5_DIR}/lib
HDF5_OPTS=-lhdf5
endif

ifeq ($(COMPILE_ON_SYSTEM),"OzSTAR")
CC=mpicc
FFTW_LIBS=${EBROOTFFTW}/lib
FFTW_INCL=${EBROOTFFTW}/include
FFTW_OPTS=-ldrfftw_mpi -ldrfftw -ldfftw_mpi -ldfftw -DDOUBLE_FFTW
PNG_LIBS=/opt/local/lib
PNG_INCL=/opt/local/include
PNG_OPTS=-lpng
HDF5_INCL=${EBROOTHDF5}/include -D H5_USE_16_API
HDF5_LIBS=${EBROOTHDF5}/lib
HDF5_OPTS=-lhdf5
endif

ifeq ($(COMPILE_ON_SYSTEM),"MacBook")
CC=/opt/homebrew/bin/gcc-15 -fopenmp
FFTW_LIBS=/opt/local/lib
FFTW_INCL=/opt/local/include
FFTW_OPTS=-ldrfftw_mpi -ldrfftw -ldfftw_mpi -ldfftw -DDOUBLE_FFTW
PNG_LIBS=/opt/homebrew/lib
PNG_INCL=/opt/homebrew/include
PNG_OPTS=-lpng
HDF5_INCL=/opt/homebrew/include -D H5_USE_16_API
HDF5_LIBS=/opt/homebrew/lib
HDF5_OPTS=-lhdf5 -lz
endif

ifeq ($(COMPILE_ON_SYSTEM),"MacPro")
CC=mpicc -g
FFTW_LIBS=/usr/local/lib
FFTW_INCL=/usr/local/include
FFTW_OPTS=-ldrfftw_mpi -ldrfftw -ldfftw_mpi -ldfftw -DDOUBLE_FFTW
PNG_LIBS=/usr/local/lib
PNG_INCL=/usr/local/include
PNG_OPTS=-lpng
HDF5_INCL=/usr/local/include -D H5_USE_16_API
HDF5_LIBS=/usr/local/lib
HDF5_OPTS=-lhdf5 -lz
endif

ifeq ($(COMPILE_ON_SYSTEM),"Magnus")
CC=cc
FFTW_LIBS=${FFTW_DIR}
FFTW_INCL=${FFTW_INC}
FFTW_OPTS=-ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw -DDOUBLE_FFTW
PNG_LIBS=
PNG_INCL=
PNG_OPTS=
HDF5_INCL=${HDF5_DIR}/include
HDF5_LIBS=${HDF5_DIR}/lib
HDF5_OPTS=-lhdf5
endif

OPTS=-lm $(OPT) $(PNG_OPTS) $(HDF5_OPTS)

render_image.exe : render_image.c io.c find_neighbours.c make_tree.c walk_tree.c kernels.c split_across_tasks.c select_particles.c header.c smooth_to_mesh.c deposit_sph_2d.c write_to_ppm.c colormap.c
	$(CC) -o render_image.exe -ffast-math -O3 $(OPTS) -I$(HDF5_INCL) -I$(PNG_INCL) render_image.c io.c find_neighbours.c make_tree.c walk_tree.c kernels.c split_across_tasks.c select_particles.c header.c smooth_to_mesh.c deposit_sph_2d.c write_to_ppm.c colormap.c -L$(HDF5_LIBS) $(HDF5_OPTS) -L$(PNG_LIBS) $(PNG_OPTS) 

clean: 
	rm render_image.exe
	touch makefile
