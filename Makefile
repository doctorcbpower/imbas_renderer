OPT += -DCIC
OPT += -DIMAGE_DIMENSIONX=4096
OPT += -DIMAGE_DIMENSIONY=4096
#OPT += -DMAX_H_VOXELS=8
OPT += -DMAX_H_DEPOSIT_PX=32
OPT += -DMIN_RHO_FRAC=0.05
#OPT += -DLONG_IDS
#OPT += -DVELOCITIES
#OPT += -DIMAGE
#OPT += -DDEBUG
OPT += -DENABLE_OPENMP
#OPT += -DENABLE_MPI
#OPT += -DSCATTER_DECOMPOSITION
#OPT += -DNONPERIODIC
OPT += -DKERNEL_SMOOTHING
#OPT += -DSMOOTH_GRID=512

COMPILE_ON_SYSTEM="MacBook"
#COMPILE_ON_SYSTEM="MacPro"
#COMPILE_ON_SYSTEM="OzSTAR"
#COMPILE_ON_SYSTEM="Setonix"
#COMPILE_ON_SYSTEM="Naranjo"

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
CC=mpicc -fopenmp
PNG_LIBS=/opt/local/lib
PNG_INCL=/opt/local/include
PNG_OPTS=-lpng
HDF5_INCL=${EBROOTHDF5}/include -D H5_USE_16_API
HDF5_LIBS=${EBROOTHDF5}/lib
HDF5_OPTS=-lhdf5
endif

ifeq ($(COMPILE_ON_SYSTEM),"Naranjo")
CC=mpicc -fopenmp
PNG_PREFIX=/home/bsreedhar/software/modules/libpng/1.6.44
PNG_LIBS=$(PNG_PREFIX)/lib
PNG_INCL=$(PNG_PREFIX)/include
PNG_OPTS=-lpng
HDF5_INCL=${EBROOTHDF5}/include -D H5_USE_16_API
HDF5_LIBS=${EBROOTHDF5}/lib
HDF5_OPTS=-lhdf5
endif

ifeq ($(COMPILE_ON_SYSTEM),"MacBook")
CC=/opt/homebrew/bin/gcc-15 -fopenmp
PNG_LIBS=/opt/homebrew/lib
PNG_INCL=/opt/homebrew/include
PNG_OPTS=-lpng
HDF5_INCL=/opt/homebrew/include -D H5_USE_16_API
HDF5_LIBS=/opt/homebrew/lib
HDF5_OPTS=-lhdf5 -lz
endif

ifeq ($(COMPILE_ON_SYSTEM),"MacPro")
CC=mpicc -g -fopenmp
PNG_LIBS=/usr/local/lib
PNG_INCL=/usr/local/include
PNG_OPTS=-lpng
HDF5_INCL=/usr/local/include -D H5_USE_16_API
HDF5_LIBS=/usr/local/lib
HDF5_OPTS=-lhdf5 -lz
endif

ifeq ($(COMPILE_ON_SYSTEM),"Magnus")
CC=cc -fopenmp
PNG_LIBS=
PNG_INCL=
PNG_OPTS=
HDF5_INCL=${HDF5_DIR}/include
HDF5_LIBS=${HDF5_DIR}/lib
HDF5_OPTS=-lhdf5
endif

INCL = -Isrc/core \
       -Isrc/io \
       -Isrc/tree \
       -Isrc/sph \
       -Isrc/render \
       -Isrc/parallel

SRCS = src/core/render_image.c \
       src/io/io.c \
       src/io/select_particles.c \
       src/tree/make_tree.c \
       src/tree/walk_tree.c \
       src/sph/kernels.c \
       src/sph/find_neighbours.c \
       src/sph/deposit_sph_2d.c \
       src/sph/smooth_to_mesh.c \
       src/render/colormap.c \
       src/render/write_to_ppm.c \
       src/parallel/split_across_tasks.c \
       src/parallel/header.c

OPTS = -lm $(OPT) $(PNG_OPTS) $(HDF5_OPTS)

render_image : $(SRCS)
	$(CC) -o render_image -ffast-math -O3 $(OPTS) \
	    -I$(HDF5_INCL) -I$(PNG_INCL) $(INCL) \
	    $(SRCS) \
	    -L$(HDF5_LIBS) $(HDF5_OPTS) -L$(PNG_LIBS) $(PNG_OPTS)

clean:
	rm -f render_image
	touch Makefile
