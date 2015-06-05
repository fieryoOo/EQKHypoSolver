
BIN  = Test

# --- compilers --- #
CC = g++

FC = gfortran

DBG =

LIBS = -lstdc++ -lfftw3 -fopenmp

# --- flags --- #
cflags = -std=c++11 -O3 
#$(OMPflag) -lX11 -lm -rdynamic -O3
#fflags =  -Wall -ffixed-line-length-none
fflags = -O2 -Wall -ffixed-line-length-none
#fflags =  -Wall  -O2 -m32 -ffixed-line-length-none

FFLAGS = $(DBG) $(fflags)
CFLAGS = $(DBG) $(cflags)

FOBJS = $(patsubst %.f,%.o,$(wildcard ./*.f)) $(patsubst %.cpp,%.o,$(wildcard ./*.cpp)) $(patsubst %.c,%.o,$(wildcard ./*.c))
#FOBJS = surfsyn.o source.o surfread.o intpol.o pha.o  \
	syn.o FFT.o calspecr.o calspecl.o new_tapwin1.o namer.o \
	force.o angles2tensor.o azd.o wsac_.o rw_sac.o \
	atracer.o rbimod.o read_rect_model.o tracer.o vect.o

$(BIN) : $(FOBJS)
	$(FC) -o $(BIN) $(FOBJS) $(FFLAGS) $(LIBS) $(LDLIBS) 

%.o : %.cpp
	$(CC) $(cflags) -c $< -o $@

%.o : %.f
	$(FC) $(fflags) -c $< -o $@

# ------------- #

INST_DIR = ../../bin

.PHONY : install
install :: $(BIN)
	install -s $(BIN) $(INST_DIR)

.PHONY : clean
clean :
	rm -f $(BIN) core.* $(FOBJS)
