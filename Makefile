CXX           = g++
CXXFLAGS      = -std=c++17 -Wall -fPIC -O3
LD            = g++
LDFLAGS       = -O3

LIBS          = -I//opt/rh/devtoolset-9/root/usr/local/include/ -L/opt/rh/devtoolset-9/root/usr/local/lib -lgsl -lgslcblas -I/home/gsl/2.7.1/include -L/home/gsl/2.7.1/lib -I/home/zhujie/Cuba-4.2.2 -L/home/zhujie/Cuba-4.2.2 -lcuba -lm `lhapdf-config --cflags --ldflags`#-lcuba -lm

BUILD = $(PWD)/build
MCDIP = $(PWD)

vpath %.cpp src
objdir     = obj
tabdir     = tabs

SRC        = main.cpp config.cpp pdfs_lhapdf.cpp nucleus.cpp event.cpp charges.cpp dipole_ipsat.cpp model_gbw.cpp model_ipsat.cpp random.cpp
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
	$(CXX) $(CXXFLAGS) -c $< -o $@

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
		@rm -f $(OBJS) $(TARGET)

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
