# =============================================================================
# render_image — Makefile
# =============================================================================

# ---------------------------------------------------------------------------
# Feature flags  (comment/uncomment to toggle)
# ---------------------------------------------------------------------------
OPT += -DCIC
OPT += -DCIC_STRIP_HEIGHT=0
OPT += -DIMAGE_DIMENSIONX=1024
OPT += -DIMAGE_DIMENSIONY=1024
#OPT += -DMAX_H_VOXELS=8
OPT += -DMAX_H_DEPOSIT_PX=32
OPT += -DMIN_RHO_FRAC=1.05
#OPT += -DLONG_IDS
#OPT += -DVELOCITIES
#OPT += -DIMAGE
#OPT += -DDEBUG
OPT += -DENABLE_OPENMP
#OPT += -DENABLE_MPI
#OPT += -DSCATTER_DECOMPOSITION
#OPT += -DNONPERIODIC
#OPT += -DKERNEL_SMOOTHING
#OPT += -DSMOOTH_GRID=512

# ---------------------------------------------------------------------------
# Target system  (uncomment exactly one)
# ---------------------------------------------------------------------------
#SYSTEM ?= MacBook
#SYSTEM ?= MacPro
#SYSTEM ?= OzSTAR
#SYSTEM ?= Setonix
SYSTEM ?= Naranjo
#SYSTEM ?= Magnus

# ---------------------------------------------------------------------------
# Per-system compiler + library configuration
# ---------------------------------------------------------------------------
ifeq ($(SYSTEM),Setonix)
    CC       := cc -fopenmp
    PNG_HOME := /software/projects/pawsey1164/cpower
    PNG_LIBS := $(PNG_HOME)/lib
    PNG_INCL := $(PNG_HOME)/include
    PNG_OPTS := -lpng
    HDF5_INCL := $(HDF5_DIR)/include -DH5_USE_16_API
    HDF5_LIBS := $(HDF5_DIR)/lib
    HDF5_OPTS := -lhdf5
endif

ifeq ($(SYSTEM),OzSTAR)
    CC       := mpicc -fopenmp
    PNG_LIBS := /opt/local/lib
    PNG_INCL := /opt/local/include
    PNG_OPTS := -lpng
    HDF5_INCL := $(EBROOTHDF5)/include -DH5_USE_16_API
    HDF5_LIBS := $(EBROOTHDF5)/lib
    HDF5_OPTS := -lhdf5
endif

ifeq ($(SYSTEM),Naranjo)
    CC       := mpicc -fopenmp
    PNG_HOME := /home/bsreedhar/software/modules/libpng/1.6.44
    PNG_LIBS := $(PNG_HOME)/lib
    PNG_INCL := $(PNG_HOME)/include
    PNG_OPTS := -lpng
    HDF5_INCL := $(EBROOTHDF5)/include -DH5_USE_16_API
    HDF5_LIBS := $(EBROOTHDF5)/lib
    HDF5_OPTS := -lhdf5
endif

ifeq ($(SYSTEM),MacBook)
    CC       := /opt/homebrew/bin/gcc-15 -fopenmp
    export MACOSX_DEPLOYMENT_TARGET := 26.4
    PNG_LIBS := /opt/homebrew/lib
    PNG_INCL := /opt/homebrew/include
    PNG_OPTS := -lpng
    HDF5_INCL := /opt/homebrew/include -DH5_USE_16_API
    HDF5_LIBS := /opt/homebrew/lib
    HDF5_OPTS := -lhdf5 -lz
    YAML_INCL := /opt/homebrew/include
    YAML_LIBS := /opt/homebrew/lib
    YAML_OPTS := -lyaml
endif

ifeq ($(SYSTEM),MacPro)
    CC       := mpicc -g -fopenmp
    PNG_LIBS := /usr/local/lib
    PNG_INCL := /usr/local/include
    PNG_OPTS := -lpng
    HDF5_INCL := /usr/local/include -DH5_USE_16_API
    HDF5_LIBS := /usr/local/lib
    HDF5_OPTS := -lhdf5 -lz
    YAML_INCL := /usr/local/include
    YAML_LIBS := /usr/local/lib
    YAML_OPTS := -lyaml
endif

ifeq ($(SYSTEM),Magnus)
    CC       := cc -fopenmp
    PNG_LIBS :=
    PNG_INCL :=
    PNG_OPTS :=
    HDF5_INCL := $(HDF5_DIR)/include
    HDF5_LIBS := $(HDF5_DIR)/lib
    HDF5_OPTS := -lhdf5
endif

# Fall back to pkg-config for YAML on any system that didn't set it above
YAML_INCL ?= $(shell pkg-config --cflags-only-I yaml-0.1 2>/dev/null | sed 's/-I//g')
YAML_LIBS ?= $(shell pkg-config --libs-only-L  yaml-0.1 2>/dev/null | sed 's/-L//g')
YAML_OPTS ?= $(shell pkg-config --libs-only-l  yaml-0.1 2>/dev/null)

# ---------------------------------------------------------------------------
# Compilation flags
# ---------------------------------------------------------------------------
#CFLAGS  := -Wall -Wextra -ffast-math -O3 -std=c99 $(OPT)
CFLAGS  := -Wall -Wextra -ffast-math -O3 -std=gnu17 $(OPT)

# ---------------------------------------------------------------------------
# Include paths
# ---------------------------------------------------------------------------
INCL := \
    -Isrc/core   \
    -Isrc/io     \
    -Isrc/tree   \
    -Isrc/sph    \
    -Isrc/render \
    -Isrc/parallel

ALL_INCL := $(INCL) \
    -I$(HDF5_INCL) \
    -I$(PNG_INCL)  \
    $(if $(YAML_INCL),-I$(YAML_INCL))

# ---------------------------------------------------------------------------
# Sources, objects, deps  (auto-discovered)
# ---------------------------------------------------------------------------
SRC_DIR   := src
BUILD_DIR := build
BIN       := render_image

SRCS := $(shell find $(SRC_DIR) -name '*.c')
OBJS := $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o, $(SRCS))
DEPS := $(OBJS:.o=.d)

-include $(DEPS)

# ---------------------------------------------------------------------------
# Targets
# ---------------------------------------------------------------------------
.PHONY: all force clean help

all: $(BIN)

$(BIN): $(OBJS)
	$(CC) -o $@ $^ \
	    -L$(HDF5_LIBS) $(HDF5_OPTS) \
	    -L$(PNG_LIBS)  $(PNG_OPTS)  \
	    $(if $(YAML_LIBS),-L$(YAML_LIBS)) $(YAML_OPTS) \
	    -lm

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(ALL_INCL) -MMD -MP -c $< -o $@

force:
	touch $(SRCS)
	$(MAKE)

clean:
	rm -rf $(BUILD_DIR) $(BIN)

help:
	@echo "Usage:  make [SYSTEM=<target>]"
	@echo ""
	@echo "Available systems:  MacBook  MacPro  OzSTAR  Setonix  Naranjo  Magnus"
	@echo "Default:            MacBook"
	@echo ""
	@echo "Examples:"
	@echo "  make                        # build for MacBook"
	@echo "  make SYSTEM=Setonix         # build for Setonix"
	@echo "  make clean                  # remove build artefacts"
