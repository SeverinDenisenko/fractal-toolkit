# Generic parameters
PROJ_NAME := fractaltoolkit

# Compiler and flags
FC       := gfortran
CFLAGS   := -O2 -g -fbacktrace -ffpe-trap=invalid,overflow -Wall -Wextra -fcheck=all -std=f2023
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
			$(SRC_DIR)/stat.f90 \
			$(SRC_DIR)/autoreg.f90 \
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

# Python extension module
PY_DIR  := python
PY_MOD  := fractaltoolkit
# The Accelerate framework flag goes through LDFLAGS below; f2py rejects
# -framework/-Wl options on its command line.
PY_STDLIB_LIBS := $(patsubst -framework Accelerate,,$(STDLIB_LIBS))

# Pip package
PYTHON     ?= python3
PKG_STAGE  := build/pip
DIST_DIR   := dist
# Extract version string from the Fortran `ver` subroutine in src/version.f90
PY_VERSION := $(shell sed -n "s/.*v = '\([^']*\)'.*/\1/p" $(SRC_DIR)/version.f90)

# Default target
all: $(APP_EXES)

.PHONY: python
python: $(LIB)
	mkdir -p $(PY_DIR)
	cd $(PY_DIR) && \
	f90wrap -k kindmap.f2cmap -m $(PY_MOD) $(addprefix ../,$(LIB_SRCS)) && \
	FFLAGS="-I$(CURDIR)/mod" LDFLAGS="-Wl,-framework,Accelerate" \
	f2py-f90wrap -c -m _$(PY_MOD) f90wrap_*.f90 $(CURDIR)/$(LIB) $(PY_STDLIB_LIBS) && \
	test -n "$$_$(PY_MOD)"*.so

# Build a pip-installable wheel from the compiled extension.
# Stages the generated wrapper and the .so into a package layout, then runs
# the PEP 517 build configured in pyproject.toml.
.PHONY: pip
pip: python
	mkdir -p $(PKG_STAGE)/$(PY_MOD)
	sed -e 's/^import _$(PY_MOD)$$/from . import _$(PY_MOD)/' \
		$(PY_DIR)/$(PY_MOD).py > $(PKG_STAGE)/$(PY_MOD)/__init__.py
	cp $(PY_DIR)/_$(PY_MOD)*.so $(PKG_STAGE)/$(PY_MOD)/
	printf '__version__ = "%s"\n' "$(PY_VERSION)" > $(PKG_STAGE)/$(PY_MOD)/_version.py
	$(PYTHON) -m pip wheel --no-deps --wheel-dir $(DIST_DIR) .
	@echo "Wheel built: $$(ls $(DIST_DIR)/$(PROJ_NAME)-*.whl)"

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
$(OBJ_DIR)/autoreg.o: $(OBJ_DIR)/precision.o $(OBJ_DIR)/constants.o
$(OBJ_DIR)/stat.o: $(OBJ_DIR)/precision.o
$(OBJ_DIR)/solvers.o: $(OBJ_DIR)/precision.o
$(OBJ_DIR)/fourier.o: $(OBJ_DIR)/precision.o $(OBJ_DIR)/complex.o $(OBJ_DIR)/integers.o $(OBJ_DIR)/constants.o
$(OBJ_DIR)/spectra.o: $(OBJ_DIR)/precision.o $(OBJ_DIR)/integers.o $(OBJ_DIR)/autoreg.o
$(OBJ_DIR)/hurst.o: $(OBJ_DIR)/precision.o $(OBJ_DIR)/solvers.o $(OBJ_DIR)/spectra.o
$(OBJ_DIR)/conv.o: $(OBJ_DIR)/precision.o $(OBJ_DIR)/fourier.o $(OBJ_DIR)/integers.o
$(OBJ_DIR)/frac.o: $(OBJ_DIR)/precision.o $(OBJ_DIR)/conv.o
$(OBJ_DIR)/generators.o: $(OBJ_DIR)/precision.o $(OBJ_DIR)/frac.o

# Create directories
$(OBJ_DIR) $(MOD_DIR) $(LIB_DIR) $(BIN_DIR):
	mkdir -p $@

clean:
	rm -rf $(OBJ_DIR) $(MOD_DIR) $(LIB_DIR) $(BIN_DIR)
	rm -rf $(DIST_DIR) $(PKG_STAGE)
	rm -f $(PY_DIR)/f90wrap_*.f90 $(PY_DIR)/fractaltoolkit.py $(PY_DIR)/*.so
	rm -rf $(PY_DIR)/__pycache__

.PHONY: all test clean
