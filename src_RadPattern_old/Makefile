BIN = PredSourcePattern

fflags = -e $(dbg) -ffixed-line-length-132
# -Nl30
cflags = $(dbg) -std=c++0x #-I$(X11_INCLUDE)

dbg = -g

FC = gfortran

CC = g++

ldlibs = -lstdc++ -lX11 -fopenmp -lm
# -lmath -lev -lio -lmap
OBJS = $(BIN)_submain.o RadPattern.o rad_pattern4_Love.o rad_pattern4_Rayl.o sourceRad.o surfreadRad.o \
       angles2tensorRad.o phaRad.o unwrapR.o unwrapL.o
#intpol.o unwrap_contin.o
#include $(DSAPMAKE)

$(BIN) : $(OBJS)
	$(FC) $^ -o $(BIN) $(ldlibs)

%.o : %.f
	$(FC) $(fflags) -c $< -o $@

%.o : %.cpp
	$(CC) $(cflags) -c $< -o $@

clean :
	rm -f $(OBJS) $(BIN)
