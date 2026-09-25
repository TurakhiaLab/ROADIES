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
CC  := gcc
REPO_ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
DIPPER_DIR ?= third_party/dipper
DIPPER_SRC_DIR ?= $(DIPPER_DIR)/src

# ==========================================
# Debug toggle
# Usage: make DEBUG=1
# ==========================================
DEBUG ?= 0
# Production builds use double precision by default.
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
# coraxlib (pinned model/optimization primitives)
# ==========================================
CORAX_REV := a2c8c6984438ed1e7aadd354f8396a1ddbdc3b69
CORAX_DIR := $(BUILD_DIR)/vendor/coraxlib-$(CORAX_REV)
CORAX_BUILD_DIR := $(CORAX_DIR)/build
CORAX_LIB := $(CORAX_BUILD_DIR)/src/corax/libcorax.a

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
            -I$(DIPPER_SRC_DIR) \
            -I$(PLL_INC_DIR) \
            -I$(CORAX_DIR)/src \
            -I$(CUDA_HOME)/include \
            -I$(CONDA_PREFIX)/include

BOOST_LINK_FLAGS := -lboost_program_options
ifdef CONDA_PREFIX
BOOST_LINK_FLAGS := -L$(CONDA_PREFIX)/lib $(BOOST_LINK_FLAGS)
endif
ZLIB_LINK_FLAGS := -lz

# ==========================================
# Flags
# ==========================================
CXXFLAGS := -O3 -std=c++17 -fno-omit-frame-pointer $(INCLUDES)
DEPFLAGS := -MMD -MP

NVCCFLAGS := -O3 -std=c++17 -ccbin $(CXX) \
             -Xcompiler -fno-omit-frame-pointer -rdc=true --extended-lambda \
             $(INCLUDES) $(CUDA_INC) $(NVCC_ARCH)

CXXFLAGS += -DTBB_PREVIEW_GLOBAL_CONTROL=1
NVCCFLAGS += -DTBB_PREVIEW_GLOBAL_CONTROL=1

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
             src/gpu/gpu_admission.cpp \
             src/io/input_validation.cpp \
             src/spr/nni.cpp \
             src/spr/local_spr.cpp \
             src/spr/local_spr_scoring.cpp \
             src/spr/local_spr_topology.cpp \
             src/util/model_utils.cpp \
             src/optimize/model_optimization_backend.cpp \
             src/optimize/model_parameter_optimizer.cpp \
             src/io/CLI.cpp \
             src/mlipper_session.cpp \
             src/io/jplace.cpp \
             src/io/tree_newick.cpp \
             src/tree/tree_generation.cpp \
             src/io/parse_file.cpp \
             src/main.cpp \
             src/util/msa_preprocess.cpp \
             src/tree/divide_and_conquer.cpp \
             src/workflow/dipper_starting_tree.cpp \

CUDA_SRCS := src/tree/tree_generation_device.cu \
             src/placement/placement.cu \
             src/placement/derivative.cu \
             src/pmatrix/pmat_gpu.cu \
             src/likelihood/partials.cu \
             src/likelihood/placement_likelihood.cu \
             src/likelihood/root_loglikelihood.cu

DIPPER_CPU_SRCS := $(DIPPER_SRC_DIR)/twoBitCompressor.cpp \
                   $(DIPPER_SRC_DIR)/fourBitCompressor.cpp \
                   $(DIPPER_SRC_DIR)/tree.cpp \
                   $(DIPPER_SRC_DIR)/divide_and_conquer/mash.cpp \
                   $(DIPPER_SRC_DIR)/divide_and_conquer/placement_close_k.cpp

DIPPER_CUDA_SRCS := $(DIPPER_SRC_DIR)/placement.cu \
                    $(DIPPER_SRC_DIR)/mash.cu \
                    $(DIPPER_SRC_DIR)/matrix_reader.cu \
                    $(DIPPER_SRC_DIR)/MSA.cu \
                    $(DIPPER_SRC_DIR)/neighborJoining.cu \
                    $(DIPPER_SRC_DIR)/placement_close_k.cu \
                    $(DIPPER_SRC_DIR)/divide_and_conquer/msa.cu

CPU_OBJS  := $(CPU_SRCS:.cpp=.o)
CUDA_OBJS := $(CUDA_SRCS:.cu=.o)
CPU_OBJS  := $(patsubst %,$(BUILD_DIR)/%,$(CPU_OBJS))
CUDA_OBJS := $(patsubst %,$(BUILD_DIR)/%,$(CUDA_OBJS))
DIPPER_CPU_OBJS := $(patsubst %,$(BUILD_DIR)/%,$(DIPPER_CPU_SRCS:.cpp=.o))
DIPPER_CUDA_OBJS := $(patsubst %,$(BUILD_DIR)/%,$(DIPPER_CUDA_SRCS:.cu=.o))
DIPPER_DC_MASH_CUDA_OBJ := $(BUILD_DIR)/dipper_divide_and_conquer_mash_cuda.o
DIPPER_DC_PLACEMENT_CLOSE_K_CUDA_OBJ := $(BUILD_DIR)/dipper_divide_and_conquer_placement_close_k_cuda.o
OBJS      := $(CPU_OBJS) $(CUDA_OBJS)

ALL_OBJS := $(OBJS) \
            $(DIPPER_CPU_OBJS) \
            $(DIPPER_CUDA_OBJS) \
            $(DIPPER_DC_MASH_CUDA_OBJ) \
            $(DIPPER_DC_PLACEMENT_CLOSE_K_CUDA_OBJ)
DEPFILES := $(ALL_OBJS:.o=.d)

TARGET := MLIPPER
MSA_PREPROCESS_TEST := $(BUILD_DIR)/tests/msa_preprocess_test
MODEL_UTILS_TEST := $(BUILD_DIR)/tests/model_utils_test
INPUT_VALIDATION_TEST := $(BUILD_DIR)/tests/input_validation_test
TREE_IO_TEST := $(BUILD_DIR)/tests/tree_io_test
TOPOLOGY_REFINEMENT_TEST := $(BUILD_DIR)/tests/topology_refinement_test
GPU_ADMISSION_TEST := $(BUILD_DIR)/tests/gpu_admission_test
CHECKED_SIZE_TEST := $(BUILD_DIR)/tests/checked_size_test
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
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -MF $(@:.o=.d) -MT $@ -c $< -o $@

$(CORAX_LIB):
	@if [ ! -d "$(CORAX_DIR)/.git" ]; then \
		mkdir -p $(dir $(CORAX_DIR)); \
		git clone --quiet --branch v2.0.2 \
			https://codeberg.org/Exelixis-Lab/coraxlib.git $(CORAX_DIR); \
	fi
	@test "$$(git -C $(CORAX_DIR) rev-parse HEAD)" = "$(CORAX_REV)"
	cmake -S $(CORAX_DIR) -B $(CORAX_BUILD_DIR) \
		-DCMAKE_BUILD_TYPE=Release \
		-DCORAX_BUILD_TESTS=OFF \
		-DCORAX_BUILD_BENCHMARKS=OFF \
		-DCORAX_BUILD_DOCS=OFF
	cmake --build $(CORAX_BUILD_DIR) --target corax --parallel 4

# These translation units include corax public headers, so the pinned source
# tree must exist before they are compiled in a parallel build.
$(BUILD_DIR)/src/util/model_utils.o $(BUILD_DIR)/src/optimize/model_parameter_optimizer.o: | $(CORAX_LIB)

$(DIPPER_DC_MASH_CUDA_OBJ): $(DIPPER_SRC_DIR)/divide_and_conquer/mash.cu
	@mkdir -p $(dir $@)
	$(NVCC) $(NVCCFLAGS) $(DEPFLAGS) -MF $(@:.o=.d) -MT $@ -c $< -o $@

$(DIPPER_DC_PLACEMENT_CLOSE_K_CUDA_OBJ): $(DIPPER_SRC_DIR)/divide_and_conquer/placement_close_k.cu
	@mkdir -p $(dir $@)
	$(NVCC) $(NVCCFLAGS) $(DEPFLAGS) -MF $(@:.o=.d) -MT $@ -c $< -o $@

# ---- CUDA -> .o ----
$(BUILD_DIR)/%.o: %.cu
	@mkdir -p $(dir $@)
	$(NVCC) $(NVCCFLAGS) $(DEPFLAGS) -MF $(@:.o=.d) -MT $@ -c $< -o $@

# ---- Link ----
$(TARGET): $(ALL_OBJS) $(CORAX_LIB)
	@mkdir -p $(dir $@)
	$(NVCC) -ccbin $(CXX) $(NVCC_ARCH) -o $@ $^ \
			$(BOOST_LINK_FLAGS) \
			$(ZLIB_LINK_FLAGS) \
			$(PLL_LINK_FLAGS) \
			$(CUDA_LIB) -lcudart -lcuda -lnvidia-ml \
			-ltbb \
			$(LAPACK_LIBS)

# Each object tracks the headers it actually includes. Missing dependency files
# are expected before the first build.
$(ALL_OBJS): $(BUILD_CONFIG)
-include $(DEPFILES)

debug:
	$(MAKE) clean
	$(MAKE) DEBUG=1 $(TARGET)

run: $(TARGET)
	./$(TARGET) $(RUN_ARGS)

$(MSA_PREPROCESS_TEST): tests/msa_preprocess_test.cpp src/util/msa_preprocess.cpp src/util/msa_preprocess.hpp src/tree/tree.hpp $(BUILD_CONFIG)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) tests/msa_preprocess_test.cpp src/util/msa_preprocess.cpp -o $@ -ltbb

$(MODEL_UTILS_TEST): tests/model_utils_test.cpp src/util/model_utils.cpp src/util/model_utils.hpp src/io/parse_file.hpp src/tree/tree.hpp $(BUILD_CONFIG) $(CORAX_LIB)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) tests/model_utils_test.cpp src/util/model_utils.cpp $(CORAX_LIB) -o $@

$(INPUT_VALIDATION_TEST): tests/input_validation_test.cpp src/io/input_validation.cpp src/io/input_validation.hpp src/io/parse_file.hpp $(BUILD_CONFIG)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) tests/input_validation_test.cpp src/io/input_validation.cpp -o $@ $(PLL_LINK_FLAGS)

$(TREE_IO_TEST): tests/tree_io_test.cpp src/io/parse_file.cpp src/io/tree_newick.cpp src/io/jplace.cpp src/io/parse_file.hpp src/io/tree_newick.hpp src/io/jplace.hpp $(BUILD_CONFIG)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) tests/tree_io_test.cpp src/io/parse_file.cpp src/io/tree_newick.cpp src/io/jplace.cpp -o $@ $(PLL_LINK_FLAGS)

$(TOPOLOGY_REFINEMENT_TEST): tests/topology_refinement_test.cpp src/spr/local_spr_topology.cpp src/spr/nni.cpp src/tree/divide_and_conquer.cpp src/spr/local_spr_internal.hpp src/spr/nni.hpp src/tree/divide_and_conquer.hpp src/tree/tree_topology_utils.hpp $(BUILD_CONFIG)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -ffunction-sections -fdata-sections \
		tests/topology_refinement_test.cpp \
		src/spr/local_spr_topology.cpp src/spr/nni.cpp \
		src/tree/divide_and_conquer.cpp \
		-Wl,--gc-sections -o $@

$(GPU_ADMISSION_TEST): tests/gpu_admission_test.cpp src/gpu/gpu_admission.cpp src/gpu/gpu_admission.hpp $(BUILD_CONFIG)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) tests/gpu_admission_test.cpp src/gpu/gpu_admission.cpp \
		-o $@ $(CUDA_LIB) -lcudart -lnvidia-ml -lrt

$(CHECKED_SIZE_TEST): tests/checked_size_test.cpp src/util/checked_size.hpp $(BUILD_CONFIG)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) tests/checked_size_test.cpp -o $@

test: $(MSA_PREPROCESS_TEST) $(MODEL_UTILS_TEST) $(INPUT_VALIDATION_TEST) $(TREE_IO_TEST) $(TOPOLOGY_REFINEMENT_TEST) $(GPU_ADMISSION_TEST) $(CHECKED_SIZE_TEST)
	./$(MSA_PREPROCESS_TEST)
	./$(MODEL_UTILS_TEST)
	./$(INPUT_VALIDATION_TEST)
	./$(TREE_IO_TEST)
	./$(TOPOLOGY_REFINEMENT_TEST)
	./$(GPU_ADMISSION_TEST)
	./$(CHECKED_SIZE_TEST)
	python3 tests/compare_workflow_outputs_test.py

test-gpu-integration: $(TARGET)
	bash tests/run_workflow_regression.sh

test-workflows: test-gpu-integration

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

.PHONY: clean debug run test test-gpu-integration test-workflows cuda-gdb double float
