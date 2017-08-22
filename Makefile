UNAME := $(shell uname)
NVCC := $(shell command -v nvcc 2> /dev/null)

DIR_MAIN       = ./
DIR_SRC        = $(DIR_MAIN)rhic
DIR_H          = $(DIR_MAIN)include/
DIR_BUILD      = $(DIR_MAIN)build/
DIR_OBJ        = $(DIR_BUILD)rhic

DEBUG =
FLOWTRACE =
LDFLAGS=
CFLAGS = $(DEBUG) $(OPTIMIZATION) $(FLOWTRACE) $(OPTIONS) -Wno-comment

ifdef NVCC
COMPILER = nvcc
OPTIMIZATION = -O5
OPTIONS := $(OPTIONS) --relocatable-device-code=true -Wno-deprecated-gpu-targets
LINK_OPTIONS := $(LINK_OPTIONS) --cudart static --relocatable-device-code=true -link -Wno-deprecated-gpu-targets
endif
ifndef NVCC
COMPILER = gcc
OPTIMIZATION = -O3
endif

ifeq ($(UNAME), Linux)
LIBS = -lm -lgsl -lgslcblas -lconfig -lgtest
endif
ifeq ($(UNAME), Darwin)
LIBS = -L /usr/local/lib -lm -lgsl -lgslcblas -lconfig -lgtest -largp -lc++
endif

INCLUDES =  -I /usr/local/include -I rhic/rhic-core/src/include -I rhic/rhic-harness/src/main/include -I rhic/rhic-trunk/src/include -I rhic/rhic-harness/src/include

CPP := $(shell find $(DIR_SRC) -name '*.cpp')
CPP_OBJ  = $(CPP:$(DIR_SRC)%.cpp=$(DIR_OBJ)%.o)
OBJ = $(CPP_OBJ)

EXE = cpu-vh

$(EXE): $(OBJ)
	echo "Linking:   $@ ($(COMPILER))"
	$(COMPILER) $(LINK_OPTIONS) -o $@ $^ $(LIBS) $(INCLUDES)
	echo "Testing:   $@"
	$(DIR_MAIN)$(EXE) --test

$(DIR_OBJ)%.o: $(DIR_SRC)%.cpp
	@[ -d $(DIR_OBJ) ] || find rhic/rhic-core rhic/rhic-harness rhic/rhic-trunk -type d -exec mkdir -p ./build/{} \;
	@echo "Compiling: $< ($(COMPILER))"
	$(COMPILER) $(CFLAGS) $(INCLUDES) -c -o $@ $<

all: $(EXE)

test: $(EXE)
	echo "Testing:   $(EXE)"
	$(DIR_MAIN)$(EXE) --test

hydro: $(EXE)
	echo "Running Hydro:   $(EXE)"
	rm -rf output
	mkdir output
	$(DIR_MAIN)$(EXE) --config=rhic-conf --output=output --hydro

clean:
	@echo "Object files and executable deleted"
	if [ -d "$(DIR_OBJ)" ]; then rm -rf $(EXE) $(DIR_OBJ)/*; rmdir $(DIR_OBJ); rmdir $(DIR_BUILD); fi

.SILENT:
