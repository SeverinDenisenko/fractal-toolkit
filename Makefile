# Generic parameters
PROJ_NAME := fractaltoolkit

# Compiler and flags
FC       := gfortran
CFLAGS   := -O2 -g -fbacktrace -Wall -Wextra -fcheck=all -std=f2023
LDFLAGS  :=
AR       := ar
ARFLAGS  := rcs

# Directories
SRC_DIR  := src
APP_DIR  := app
TEST_DIR := test
OBJ_DIR  := obj
MOD_DIR  := mod
LIB_DIR  := lib
BIN_DIR  := bin

# Module search path and output
CFLAGS   += -I$(MOD_DIR) -J$(MOD_DIR)

# Library sources in dependency order (modules first)
LIB_SRCS := $(SRC_DIR)/version.f90 \
			$(SRC_DIR)/precision.f90 \
			$(SRC_DIR)/constants.f90 \
			$(SRC_DIR)/complex.f90 \
			$(SRC_DIR)/fourier.f90 \
			$(SRC_DIR)/solvers.f90 \
			$(SRC_DIR)/integers.f90 \
			$(SRC_DIR)/conv.f90 \
			$(SRC_DIR)/frac.f90 \
			$(SRC_DIR)/generators.f90 \
			$(SRC_DIR)/spectra.f90 \
			$(SRC_DIR)/hurst.f90 \
			$(SRC_DIR)/io.f90

LIB_OBJS := $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(LIB_SRCS))
LIB      := $(LIB_DIR)/lib$(PROJ_NAME).a

# Main programs (automatically discovered from main_*.f90)
APP_SRCS := $(wildcard $(APP_DIR)/main_*.f90)
APP_EXES := $(patsubst $(APP_DIR)/main_%.f90,$(BIN_DIR)/%,$(APP_SRCS))

# Test programs (automatically discovered from test_*.f90)
TEST_SRCS := $(wildcard $(TEST_DIR)/test_*.f90)
TEST_EXES := $(patsubst $(TEST_DIR)/%.f90,$(BIN_DIR)/%,$(TEST_SRCS))

# Stdlib
STDLIB_CFLAGS := `pkg-config --cflags fortran_stdlib`
STDLIB_LIBS   := $(shell pkg-config --libs fortran_stdlib)
# On macOS pkg-config emits a .../Accelerate.framework path, which ld cannot mmap().
# Сonvert it to the -framework Accelerate flag.
STDLIB_LIBS   := $(patsubst %/Accelerate.framework,-framework Accelerate,$(STDLIB_LIBS))

CFLAGS  += $(STDLIB_CFLAGS)
LDFLAGS += $(STDLIB_LIBS)

# Default target
all: $(APP_EXES)

# Build static library
$(LIB): $(LIB_OBJS) | $(LIB_DIR)
	$(AR) $(ARFLAGS) $@ $^

# Compile library sources
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(CFLAGS) -c $< -o $@

# Compile main program (depends on library objects so module files exist)
$(OBJ_DIR)/%.o: $(APP_DIR)/%.f90 $(LIB_OBJS) | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(CFLAGS) -c $< -o $@

# Link main executables
$(BIN_DIR)/%: $(OBJ_DIR)/main_%.o $(LIB) | $(BIN_DIR)
	$(FC) $(LDFLAGS) -o $@ $< -L$(LIB_DIR) -l$(PROJ_NAME)

# Build test executables (depends on library objects so module files exist)
$(BIN_DIR)/%: $(TEST_DIR)/%.f90 $(LIB_OBJS) $(LIB) | $(BIN_DIR)
	$(FC) $(CFLAGS) $(LDFLAGS) -o $@ $< -L$(LIB_DIR) -l$(PROJ_NAME)

# Run all tests
test: $(TEST_EXES)
	@for t in $(TEST_EXES); do \
	    echo "Running $$t"; \
	    $$t || exit 1; \
	done
	@echo "All tests passed."

# Explicit module dependencies (for parallel builds)
$(OBJ_DIR)/constants.o: $(OBJ_DIR)/precision.o
$(OBJ_DIR)/complex.o: $(OBJ_DIR)/constants.o
$(OBJ_DIR)/io.o: $(OBJ_DIR)/precision.o
$(OBJ_DIR)/solvers.o: $(OBJ_DIR)/precision.o
$(OBJ_DIR)/fourier.o: $(OBJ_DIR)/precision.o $(OBJ_DIR)/complex.o $(OBJ_DIR)/integers.o
$(OBJ_DIR)/hurst.o: $(OBJ_DIR)/precision.o $(OBJ_DIR)/solvers.o 
$(OBJ_DIR)/integers.o: | $(OBJ_DIR) $(MOD_DIR)
$(OBJ_DIR)/conv.o: $(OBJ_DIR)/precision.o $(OBJ_DIR)/fourier.o $(OBJ_DIR)/integers.o
$(OBJ_DIR)/frac.o: $(OBJ_DIR)/precision.o $(OBJ_DIR)/conv.o
$(OBJ_DIR)/generators.o: $(OBJ_DIR)/precision.o $(OBJ_DIR)/frac.o
$(OBJ_DIR)/spectra.o: $(OBJ_DIR)/precision.o $(OBJ_DIR)/integers.o

# Create directories
$(OBJ_DIR) $(MOD_DIR) $(LIB_DIR) $(BIN_DIR):
	mkdir -p $@

clean:
	rm -rf $(OBJ_DIR) $(MOD_DIR) $(LIB_DIR) $(BIN_DIR)

.PHONY: all test clean
