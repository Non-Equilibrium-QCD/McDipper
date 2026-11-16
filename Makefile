CXX           = g++
LD            = g++
OMPFLAGS      = -fopenmp
OMPLIBS       = -fopenmp
HDF5LIBDIR    =

# macOS: Apple clang needs libomp (brew install libomp) and -Xpreprocessor
ifeq ($(shell uname -s),Darwin)
  LIBOMP     ?= $(shell brew --prefix libomp 2>/dev/null)
  HDF5PREFIX := $(shell brew --prefix hdf5 2>/dev/null)
  OMPFLAGS    = -Xpreprocessor -fopenmp -I$(LIBOMP)/include -I$(HDF5PREFIX)/include
  OMPLIBS     = -L$(LIBOMP)/lib -lomp
  HDF5LIBDIR  = -L$(HDF5PREFIX)/lib
  CXX         = clang++
  LD          = clang++
endif

CXXFLAGS      = -std=c++17 -Wall -fPIC -O3 $(OMPFLAGS) -I/usr/local/include `gsl-config --cflags`
DEPFLAGS      = -MMD -MP
LDFLAGS       = -O3 $(OMPLIBS)

LIBS          = -L/usr/local/lib `gsl-config --libs` -lcuba -lm $(HDF5LIBDIR) -lhdf5 `lhapdf-config --cflags --ldflags`

BUILD = $(PWD)/build
MCDIP = $(PWD)

vpath %.cpp src
objdir     = obj
tabdir     = tabs

SRC        = main.cpp routines.cpp config.cpp config_mv.cpp pdfs_lhapdf.cpp nucleus.cpp event.cpp charges.cpp dipole_ipsat.cpp dipole_mv.cpp model_gbw.cpp model_ipsat.cpp model_mv.cpp random.cpp h5_out.cpp
SRCH       = Hankel.cpp
OBJS       = $(patsubst %.cpp,$(objdir)/%.o,$(SRC))

TARGET	   = mcDipper
#------------------------------------------------------------------------------
$(TARGET): $(OBJS) $(TABS)
		$(LD)  $(LDFLAGS) $^ -o $@ $(LIBS)
		@echo "$@ done"


$(OBJS): | $(objdir)

$(objdir):
	@mkdir -p $(objdir)

$(TABS): | $(tabdir)

$(tabdir):
	@mkdir -p $(tabdir)

obj/%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -c $< -o $@

# Rebuild a .o when any header it includes changes. Without this, editing
# params_gen.h / params_ipsat.h / routines.h leaves a stale binary that does
# not reflect the change, with no warning.
-include $(OBJS:.o=.d)

buildfldr:
	mkdir $(BUILD)
	cp -r "$(PWD)/src/" "$(BUILD)/src/"
	cp -r "$(PWD)/configs/" "$(BUILD)/configs/"
	cp -r "$(PWD)/ipsat/" "$(BUILD)/ipsat/"
	CP "$(PWD)/Makefile" "$(BUILD)/Makefile"

builder:
		if [ ! -d $(BUILD) ]; then make buildfldr; fi
		cd $(BUILD)



clean:
		@rm -f $(OBJS) $(OBJS:.o=.d) $(TARGET)

again:
		make clean
		make


Hankel:
	$(CXX) $(CXXFLAGS) -o Hankel.exe src/Hankel.cpp $(LIBS)

PDFs:
	$(CXX) $(CXXFLAGS) -o pdfs.exe src/make_pdfs.cpp $(LIBS) `lhapdf-config --cflags --ldflags`

renew:
	clear
	make
