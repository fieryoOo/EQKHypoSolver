BIN1 = PredSourcePattern

#BIN2 = PredAmpSpectrum

BINall = $(BIN1) #$(BIN2)

all : $(BINall)

fflags = -e $(dbg) -ffixed-line-length-132
# -Nl30
cflags = $(dbg) -std=c++0x #-I$(X11_INCLUDE)

dbg = -g

FC = gfortran

CC = g++

ldlibs = -lstdc++ -lX11 -fopenmp -lm
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
