# GNU Make / GNU Fortran. External dependencies are never downloaded implicitly.
.DEFAULT_GOAL := all
ifeq ($(origin FC),default)
FC = gfortran
endif
PYTHON ?= python3
HOPPET_CONFIG ?= hoppet-config
LHAPDF_CONFIG ?= lhapdf-config
CHAPLIN_LIBS ?= -lchaplin
FFLAGS ?= -O2 -fPIC -ffree-line-length-none -fallow-argument-mismatch
LDFLAGS ?=
PDF_SET ?= NNPDF40_nnlo_as_01180
BUILD := build
MODDIR := $(BUILD)/mod
OBJDIR := $(BUILD)/obj
SOURCES := $(wildcard src/*.f90) $(wildcard src/*.f)
OBJECTS := $(patsubst src/%.f90,$(OBJDIR)/%.o,$(filter %.f90,$(SOURCES))) $(patsubst src/%.f,$(OBJDIR)/%.o,$(filter %.f,$(SOURCES)))
LIBOBJECTS := $(filter-out $(OBJDIR)/jetvheto.o,$(OBJECTS))
INCLUDES = -J$(MODDIR) -I$(MODDIR) -Isrc -Idata/TMDs_ptj $(shell $(HOPPET_CONFIG) --fflags)
LIBS = $(shell $(HOPPET_CONFIG) --ldflags) $(shell $(LHAPDF_CONFIG) --ldflags) $(CHAPLIN_LIBS)
TESTS := rad radiator coefficients prefactor n3ll profile

.PHONY: all check check-fast check-full check-cli check-algebra regenerate-rad clean dist
all: jetvheto

$(BUILD)/dependencies.mk: scripts/fortran_dependencies.py $(SOURCES) Makefile
	mkdir -p $(BUILD)
	$(PYTHON) scripts/fortran_dependencies.py $(SOURCES) > $@
ifneq ($(filter clean dist,$(MAKECMDGOALS)),)
else
-include $(BUILD)/dependencies.mk
endif

$(OBJDIR) $(MODDIR):
	mkdir -p $@
$(OBJDIR)/%.o: src/%.f90 | $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) $(INCLUDES) -c $< -o $@
$(OBJDIR)/%.o: src/%.f | $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) $(INCLUDES) -c $< -o $@
$(OBJDIR)/rapidity_n3ll.o: src/rad_grid.inc
$(OBJDIR)/coefficient_functions_deltaptj.o: $(wildcard data/TMDs_ptj/*.txt)
jetvheto: $(OBJECTS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJECTS) $(LIBS)

tests/test_%: tests/test_%.f90 $(LIBOBJECTS)
	$(FC) $(FFLAGS) $(INCLUDES) $(LDFLAGS) $< $(LIBOBJECTS) $(LIBS) -o $@
check-fast: $(addprefix tests/test_,$(filter-out n3ll profile,$(TESTS)))
	./tests/test_rad
	$(PYTHON) tests/test_rad_import.py
	./tests/test_radiator
	./tests/test_coefficients data/TMDs_ptj
	./tests/test_prefactor
check-full: check-fast tests/test_n3ll tests/test_profile tests/test_small_r
	./tests/test_n3ll data/TMDs_ptj $(PDF_SET)
	./tests/test_profile data/TMDs_ptj $(PDF_SET)
	./tests/test_small_r $(PDF_SET)
	./tests/test_small_r $(PDF_SET) DY
check-cli: jetvheto
	PDF_SET=$(PDF_SET) $(PYTHON) tests/test_cli.py
	PDF_SET=$(PDF_SET) $(PYTHON) tests/test_anew_orders.py
	PDF_SET=$(PDF_SET) $(PYTHON) tests/test_rapidity_truncation.py
	PDF_SET=$(PDF_SET) $(PYTHON) tests/test_small_r_cli.py
check: check-full check-cli
check-algebra:
	$(PYTHON) tests/check_rapidity_truncation_algebra.py
regenerate-rad:
	$(PYTHON) scripts/import_rad.py data/RapidityAnomalousDimension.wl src/rad_grid.inc
dist:
	$(PYTHON) scripts/package_source.py
clean:
	$(PYTHON) scripts/clean_build.py
