#
BIN  = test

#fflags =  -Wall -ffixed-line-length-none
fflags =  -Wall -O2 -ffixed-line-length-none
#fflags =  -Wall  -O2 -m32 -ffixed-line-length-none

cflags = -std=c++11 -O3

FFLAGS = $(DBG) $(fflags)

LIBS = -lstdc++ -lfftw3

FC = gfortran

CC = gcc

DBG = 

FOBJS = read_rect_model.o read_rect_model_new.o

$(BIN) : $(FOBJS) $(BIN).o
	$(FC)  $(FFLAGS) $^ -o $(BIN) $(LDLIBS) $(LIBS)

%.o : %.cpp
	$(CC) -c -o $@ $< $(cflags)

install :: $(BIN)
	install -s $(BIN) $(INST_DIR)

clean ::
	rm -f $(BIN) core $(FOBJS)
