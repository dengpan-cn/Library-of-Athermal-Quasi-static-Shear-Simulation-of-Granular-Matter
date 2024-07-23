# Define targets
EXE = a.out

cSRC := $(wildcard *.c)
cOBJ := $(patsubst %.c,build/%.o,$(cSRC))

cppSRC := $(wildcard *.cpp)
cppOBJ := $(patsubst %.cpp,build/%.o,$(cppSRC))

SRC := $(cSRC) $(cppSRC)
OBJ = $(cOBJ) $(cppOBJ)

SourceFileName = $(SRC) $(wildcard *.h) $(wildcard *.hpp)
CompilePath = $(shell pwd)

# Define compiler and its settings
ifneq ($(shell which icc 2>/dev/null),)
    # use the Intel icc compiler for C codes, and -qopenmp for OpenMP
    CC = icc -std=c99
    CFLAGS += -D_GNU_SOURCE -D __SourceFileName__="\"$(SourceFileName)\"" -D __CompilePath__="\"$(CompilePath)"\"
    CXX = icpc
    LINK = icpc
else
    CC = gcc
    CFLAGS += -D_GNU_SOURCE -D __SourceFileName__="\"$(SourceFileName)\"" -D __CompilePath__="\"$(CompilePath)"\"
    CXX = g++
    LINK = g++
endif

ifeq ($(shell uname), Linux)
	CFLAGS += -D__Linux__
	objSourceFile = build/sourceFile.obj
	objTool = objcopy
else
	CFLAGS += -D__MacOS__
	objSourceFile = 
	objTool = 
endif

#debug flag
gdb ?= -O3
#gdb = -g2 -gstabs+
ifneq ($(gdb), -O3)
    CC = gcc
    CXX = g++
    LINK = g++
endif

# -D DIM=xx: for simulation of different dimension.
# -D__triBox__: triclinic box.
# -D__orthBox__: orthogonal box.
CFLAGS += -D DIM=2 -D__triBox__

# BLAS
ifdef MKLROOT
	BLAS_Path = -DMkl
	BLAS_LIB = -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lm
else
	BLAS_INSTALL = ${HOME}/SoftWare/
	BLAS_Path = -DOpenblas -I $(BLAS_INSTALL)/inclue
	BLAS_LIB = -L $(BLAS_INSTALL)/lib -lopenblas -lm
endif

#SparseQR
SPQR_INSTALL = ${HOME}/SoftWare/
SPQR_Path = -DSpqr -I $(SPQR_INSTALL)/include/suitesparse/
SPQR_LIB = -L $(SPQR_INSTALL)/lib64/ -lcholmod -lcamd -lccolamd -lspqr -lsuitesparseconfig -lamd -lcolamd

#voro++
VORO_INSTALL= 
#${HOME}/dengPan/SoftWare/
VORO_Path = 
#-DVoroPlusPlus -I $(VORO_INSTALL)/include
VORO_LIB = 
#-L $(VORO_INSTALL)/lib -lvoro++

# Define linker and its settings
LIB_PATH += $(BLAS_Path) $(SPQR_Path) $(VORO_Path)
LIB += $(BLAS_LIB) $(SPQR_LIB) $(VORO_LIB)

#define size command
SIZE = size

# Path to src files 
vpath %.c .
vpath %.h .
vpath %.cpp .

# Link target
$(EXE): $(OBJ) $(objSourceFile)
	$(LINK) -o $(EXE) $(OBJ) $(objSourceFile) $(LIB)
	$(SIZE) $(EXE)

# Compilation rules
build/%.o: %.c
	$(CC) $(gdb) $(CFLAGS) $(LIB_PATH)  -c $< -o $@
#	@$(CC) $(CFLAGS) $(LIB_PATH) -MM $<  >  $@.dep

build/%.o: %.cpp
	$(CXX) $(gdb) $(CFLAGS) $(LIB_PATH) -c $< -o $@
#	@$(CXX) $(CFLAGS) $(LIB_PATH) -M  $< > $@.dep

$(objSourceFile): $(SourceFileName)
	echo "" > build/sourceFile
	echo "//===Starting of Source Files===" >> build/sourceFile
	for fname in $(SourceFileName); do echo "//===Source File:" $$fname "===" >> build/sourceFile; cat $$fname >> build/sourceFile; done
	$(objTool) -I binary -B i386 -O elf64-x86-64 build/sourceFile $(objSourceFile)
	echo "//===Ending of Source Files===" >> build/sourceFile
	echo "" >> build/sourceFile
	rm build/sourceFile

# Individual dependencies
DEPENDS = $(patsubst %.o,%.o.dep,$(OBJ))
sinclude $(DEPENDS)

# Clean command
clean:
	rm -rf $(OBJ)
	rm -rf $(DEPENDS)
	rm -rf $(objSourceFile)
	rm -rf $(EXE)
