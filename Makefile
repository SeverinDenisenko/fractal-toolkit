# Generic parameters
PROJ_NAME := fractime
EXEC_NAME := fractime

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
LIB_SRCS := $(SRC_DIR)/precision.f90 \
			$(SRC_DIR)/fit.f90 \
            $(SRC_DIR)/int_math.f90 \
			$(SRC_DIR)/conv.f90 \
			$(SRC_DIR)/frac.f90 \
            $(SRC_DIR)/gen_series.f90 \
			$(SRC_DIR)/spectra.f90

LIB_OBJS := $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(LIB_SRCS))
LIB      := $(LIB_DIR)/lib$(PROJ_NAME).a

# Main program
APP_SRC  := $(APP_DIR)/main.f90
APP_OBJ  := $(OBJ_DIR)/main.o
EXE      := $(BIN_DIR)/$(EXEC_NAME)

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
all: $(EXE)

# Build static library
$(LIB): $(LIB_OBJS) | $(LIB_DIR)
	$(AR) $(ARFLAGS) $@ $^

# Compile library sources
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(CFLAGS) -c $< -o $@

# Compile main program (depends on library objects so module files exist)
$(APP_OBJ): $(APP_SRC) $(LIB_OBJS) | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(CFLAGS) -c $< -o $@

# Link main executable
$(EXE): $(APP_OBJ) $(LIB) | $(BIN_DIR)
	$(FC) $(LDFLAGS) -o $@ $(APP_OBJ) -L$(LIB_DIR) -l$(PROJ_NAME)

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
# TODO

# Create directories
$(OBJ_DIR) $(MOD_DIR) $(LIB_DIR) $(BIN_DIR):
	mkdir -p $@

clean:
	rm -rf $(OBJ_DIR) $(MOD_DIR) $(LIB_DIR) $(BIN_DIR)

.PHONY: all test clean
