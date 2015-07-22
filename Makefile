# --- excutables to build --- #
BIN1 = EQKSolver
BIN2 = Auxiliary
BIN3 = PredSrcPatterns
BIN4 = MTensor
BINT = Test

BINall = $(BIN1) $(BIN2) $(BIN3) $(BIN4)
all : $(BINall)

# --- compiliers --- #
FC = gfortran
CC = g++

# --- modules --- #
MODULES	:= src_Driver src_RadPattern src_SDContainer src_Synthetic
MOD_DIRS := $(MODULES)
#MOD_DIRS := $(addprefix mod_,$(MODULES))

# --- find source files for cpp and fortran separately --- #
#SRC_DIRS := . $(MOD_DIRS)
#SRCS_C	:= $(foreach sdir,$(SRC_DIRS),$(wildcard $(sdir)/*.cpp))
#SRCS_F	:= $(foreach sdir,$(SRC_DIRS),$(wildcard $(sdir)/*.f))

#vpath %.cpp $(SRC_DIR)
#vpath %.f $(SRC_DIR)

# --- flags --- #
INCLUDES	:= $(addprefix -I,$(MOD_DIRS))
OMPflag = -fopenmp
cflags = -O3 -std=c++11 $(OMPflag) $(INCLUDES)	#-O3
fflags = -e -O2 -ffixed-line-length-132 $(OMPflag)	#-O2
LIBS = -lstdc++ $(OMPflag) -lX11 -lm -rdynamic -lfftw3 -O3

# --- objects --- #
#OBJS_C	:= $(patsubst %.cpp,%.o,$(SRCS_C))
#OBJS_F	:= $(patsubst %.f,%.o,$(SRCS_F))
OBJS_M	:= $(addsuffix .o,$(MOD_DIRS))
OBJS	:= $(OBJS_M)
#OBJS	+= $(patsubst %.cpp,%.o,$(wildcard ./*.cpp))

# --- rules --- #
define make-bin
$(1) : $(OBJS) $(1).o
	$(FC) $$^ -o $$@ $(LIBS)
endef
$(foreach bin,$(BINall),$(eval $(call make-bin,$(bin))))

#$(BIN1) : $(OBJS) $(BIN1).o
#	$(FC) $^ $(LIBS) -o $@

define make-module-obj
$(1).o : $(patsubst %.cpp,%.o,$(filter-out $(1)/surfsyn.cpp,$(filter-out $(1)/Test.cpp,$(wildcard $(1)/*.cpp)))) $(patsubst %.f,%.o,$(wildcard $(1)/*.f))
	ld -r $$^ -o $$@
endef
$(foreach moddir,$(MOD_DIRS),$(eval $(call make-module-obj,$(moddir))))

%.o : %.cpp
	$(CC) $(cflags) -c $< -o $@

%.o : %.f
	$(FC) $(fflags) -c $< -o $@


.PHONY : clean
clean :
	rm -f $(BINall) $(addsuffix .o,$(BINall)) $(OBJS)

.PHONY : clean-mod
clean-mod :
	rm -f $(foreach moddir,$(MOD_DIRS),$(wildcard $(moddir)/*.o))

