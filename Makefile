# --- excutables to build --- #
BIN1 = EQKSolver
BIN2 = Auxiliary
BIN3 = MomentTensor
BIN4 = PredSrcPatterns
BIN5 = MTensor
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
OBJS_M	:= $(addsuffix .o,$(MOD_DIRS))
OBJS	:= $(OBJS_M)
#OBJS	+= $(patsubst %.cpp,%.o,$(wildcard ./*.cpp))

# --- main rules --- #
define make-bin
$(1) : $(OBJS) $(1).o
	$(FC) $$^ -o $$@ $(LIBS)
endef
$(foreach bin,$(BINall),$(eval $(call make-bin,$(bin))))

# --- module object rules --- #
define make-module-obj
$(1).o : $(patsubst %.cpp,%.o,$(filter-out $(wildcard $(1)/*submain.cpp),$(wildcard $(1)/*.cpp))) $(patsubst %.f,%.o,$(wildcard $(1)/*.f))
	ld -r $$^ -o $$@
endef
$(foreach moddir,$(MOD_DIRS),$(eval $(call make-module-obj,$(moddir))))

# --- .cpp rules with dependencies assembled by gcc --- #
define make-cpp
$(shell $(CC) $(cflags) -MM -MT $(patsubst %.cpp,%.o,$(1)) $(1) | tr -d '\\\n' | awk '{print}' )
	$(CC) $(cflags) -c $$< -o $$@
endef
FSRC = $(foreach moddir,$(MOD_DIRS) .,$(wildcard $(moddir)/*.cpp))
$(foreach fcpp,$(FSRC),$(eval $(call make-cpp,$(fcpp))))
#%.o : %.cpp
#	$(CC) $(cflags) -c $< -o $@

# --- .f rules --- #
%.o : %.f
	$(FC) $(fflags) -c $< -o $@


# --- .PHONYs --- #
.PHONY : clean
clean :
	rm -f $(BINall) $(addsuffix .o,$(BINall)) $(OBJS)

.PHONY : clean-mod
ifdef moddir
clean-mod :
	rm -f $(wildcard $(moddir)/*.o) $(moddir).o
else
clean-mod :
	rm -f $(foreach moddir,$(MOD_DIRS),$(wildcard $(moddir)/*.o))
endif

