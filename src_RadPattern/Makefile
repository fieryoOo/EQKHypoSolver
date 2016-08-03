BIN1 = PredSourcePattern

BIN2 = SearchSimilarSourcePattern

#BIN3 = PredAmpSpectrum

BINall = $(BIN1) $(BIN2)

all : $(BINall)

FC = gfortran
CC = g++

MODULES	:= ../src_Driver ../src_RadPattern ../src_SDContainer ../src_Synthetic
MOD_DIRS := $(MODULES)
INCLUDES	:= $(addprefix -I,$(MOD_DIRS))

OMPflag = -fopenmp
cflags = -O3 -std=c++11 $(OMPflag) $(INCLUDES)
fflags = -e -O2 -ffixed-line-length-132
# -Nl30

ldlibs = -lstdc++ -lX11 -lm -fopenmp
# -lmath -lev -lio -lmap

OBJS_Bin := $(addsuffix _submain.o,$(BINall))
OBJS := RadPattern.o rad_pattern4_Love.o rad_pattern4_Rayl.o sourceRad.o surfreadRad.o \
       angles2tensorRad.o phaRad.o unwrapR.o unwrapL.o
#intpol.o unwrap_contin.o
#include $(DSAPMAKE)

# --- main rules --- #
define make-bin
$(1) : $(OBJS) $(1)_submain.o
	$(FC) $$^ -o $$@ $(ldlibs)
endef
$(foreach bin,$(BINall),$(eval $(call make-bin,$(bin))))

%.o : %.f
	$(FC) $(fflags) -c $< -o $@

%.o : %.cpp
	$(CC) $(cflags) -c $< -o $@

clean :
	rm -f $(OBJS) $(BINall) $(OBJS_Bin)
