# ==========================================
# CUDA
# ==========================================
CUDA_HOME ?= /usr/local/cuda-12
NVCC      := $(CUDA_HOME)/bin/nvcc
CUDA_INC  := -I$(CUDA_HOME)/include
CUDA_LIB  := -L$(CUDA_HOME)/lib64
# Default fatbin covers common Turing/Ampere targets and keeps PTX for newer GPUs.
NVCC_ARCH ?= -gencode arch=compute_75,code=sm_75 \
        -gencode arch=compute_80,code=sm_80 \
        -gencode arch=compute_86,code=sm_86 \
        -gencode arch=compute_86,code=compute_86

# ==========================================
# Compiler
# ==========================================
CXX := g++
REPO_ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))

# ==========================================
# Debug toggle
# Usage: make DEBUG=1
# ==========================================
DEBUG ?= 0
# Default released builds use single precision.
USE_DOUBLE ?= 1
LINEINFO ?= 0
BUILD_DIR ?= build
# ==========================================
# libpll
# ==========================================
PLL_INC_DIR ?= /usr/local/include
PLL_LIB_DIR ?= /usr/local/lib
PLL_LINK_FLAGS := -L$(PLL_LIB_DIR) -lpll

# ==========================================
# System LAPACK / BLAS
# ==========================================
LAPACK_LIBS := -L/usr/lib/x86_64-linux-gnu \
               -llapack -lblas -lm

# ==========================================
# Include paths
# ==========================================
INCLUDES := -I. \
            -Isrc \
            -Isrc/pmatrix \
            -I$(PLL_INC_DIR) \
            -I$(CUDA_HOME)/include \
            -I$(CONDA_PREFIX)/include

BOOST_LINK_FLAGS := -lboost_program_options
ifdef CONDA_PREFIX
BOOST_LINK_FLAGS := -L$(CONDA_PREFIX)/lib $(BOOST_LINK_FLAGS)
endif

# ==========================================
# Flags
# ==========================================
CXXFLAGS := -O3 -std=c++17 -fno-omit-frame-pointer $(INCLUDES)

NVCCFLAGS := -O3 -std=c++17 -ccbin $(CXX) \
             -Xcompiler -fno-omit-frame-pointer -rdc=true \
             $(INCLUDES) $(CUDA_INC) $(NVCC_ARCH)

ifeq ($(DEBUG),1)
    CXXFLAGS  += -g -O0
    NVCCFLAGS += -G -O0 -lineinfo
endif

ifeq ($(LINEINFO),1)
    NVCCFLAGS += -lineinfo
endif

ifeq ($(USE_DOUBLE),1)
    CXXFLAGS  += -DMLIPPER_USE_DOUBLE
    NVCCFLAGS += -DMLIPPER_USE_DOUBLE
endif

# ==========================================
# Sources
# ==========================================
CPU_SRCS  := src/pmatrix/pmat.cpp \
             src/io/input_validation.cpp \
             src/spr/local_spr.cpp \
             src/model_utils.cpp \
             src/io/jplace.cpp \
             src/io/tree_newick.cpp \
             src/tree/tree_generation.cpp \
             src/io/parse_file.cpp \
             src/main.cpp \
             src/msa_preprocess.cpp

CUDA_SRCS := src/tree/tree_generation_device.cu \
             src/placement/placement.cu \
             src/placement/derivative.cu \
             src/pmatrix/pmat_gpu.cu \
             src/likelihood/partial_likelihood.cu \
             src/likelihood/root_likelihood.cu

CPU_OBJS  := $(CPU_SRCS:.cpp=.o)
CUDA_OBJS := $(CUDA_SRCS:.cu=.o)
CPU_OBJS  := $(patsubst %,$(BUILD_DIR)/%,$(CPU_OBJS))
CUDA_OBJS := $(patsubst %,$(BUILD_DIR)/%,$(CUDA_OBJS))
OBJS      := $(CPU_OBJS) $(CUDA_OBJS)

TARGET := MLIPPER
RUN_ARGS ?=
.DEFAULT_GOAL := $(TARGET)

# ==========================================
# Rules
# ==========================================
# Build-config stamp: forces recompilation when toggling Makefile options (e.g. USE_DOUBLE/DEBUG)
BUILD_CONFIG := $(BUILD_DIR)/.build_config

FORCE:
.PHONY: FORCE

$(BUILD_CONFIG): FORCE
	@mkdir -p $(dir $(BUILD_CONFIG))
	@printf "DEBUG=%s\nUSE_DOUBLE=%s\nLINEINFO=%s\nCUDA_HOME=%s\nPLL_INC_DIR=%s\nPLL_LIB_DIR=%s\nNVCC_ARCH=%s\n" "$(DEBUG)" "$(USE_DOUBLE)" "$(LINEINFO)" "$(CUDA_HOME)" "$(PLL_INC_DIR)" "$(PLL_LIB_DIR)" "$(NVCC_ARCH)" > $(BUILD_CONFIG).tmp
	@cmp -s $(BUILD_CONFIG).tmp $(BUILD_CONFIG) 2>/dev/null || mv $(BUILD_CONFIG).tmp $(BUILD_CONFIG)
	@rm -f $(BUILD_CONFIG).tmp

# ---- C++ -> .o ----
$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ---- CUDA -> .o ----
$(BUILD_DIR)/%.o: %.cu
	@mkdir -p $(dir $@)
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

# ---- Link ----
$(TARGET): $(OBJS)
	@mkdir -p $(dir $@)
	$(NVCC) -ccbin $(CXX) $(NVCC_ARCH) -o $@ $^ \
		$(BOOST_LINK_FLAGS) \
		$(PLL_LINK_FLAGS) \
		$(CUDA_LIB) -lcudart -lcuda \
		-ltbb \
		$(LAPACK_LIBS)

# Crude dependency tracking: force rebuild when core headers change.
$(OBJS): $(BUILD_CONFIG) src/tree/tree.hpp src/io/input_validation.hpp src/spr/local_spr.hpp src/model_utils.hpp src/io/jplace.hpp src/io/parse_file.hpp src/util/precision.hpp src/io/tree_newick.hpp src/util/mlipper_util.h src/placement/derivative.cuh src/likelihood/partial_likelihood.cuh src/likelihood/root_likelihood.cuh src/placement/placement.cuh
debug:
	$(MAKE) clean
	$(MAKE) DEBUG=1 $(TARGET)

run: $(TARGET)
	./$(TARGET) $(RUN_ARGS)

cuda-gdb: $(TARGET)
	cuda-gdb --args ./$(TARGET) $(RUN_ARGS)

clean:
	rm -rf $(BUILD_DIR)

double:
	$(MAKE) clean
	$(MAKE) USE_DOUBLE=1 $(TARGET)

float:
	$(MAKE) clean
	$(MAKE) USE_DOUBLE=0 $(TARGET)

.PHONY: clean debug run cuda-gdb double float
