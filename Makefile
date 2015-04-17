#BIN1 = SearchSource_SA

#BIN2 = SearchSource_SA_MC

#BIN3 = CheckFit

#BIN4 = PredictTraveltimes

#BIN5 = ComputeMisfitsAll

#BIN6 = Auxiliary

BINTEST = Test

OBJS = SDContainer.o Map.o RadPattern.o

VPATH_Rad = RadPattern_src
OBJS_Rad = $(VPATH_Rad)/rad_pattern4_Rayl.o $(VPATH_Rad)/rad_pattern4_Love.o \
	   $(VPATH_Rad)/source.o $(VPATH_Rad)/surfread.o $(VPATH_Rad)/angles2tensor.o \
 	   $(VPATH_Rad)/pha.o $(VPATH_Rad)/unwrapR.o $(VPATH_Rad)/unwrapL.o

OMPflag = -fopenmp

LIBS = -lstdc++ $(OMPflag) -lX11 -lm -rdynamic -g

cflags = -g -std=c++0x $(OMPflag) #-O3

fflags = -e -O2 -ffixed-line-length-132 $(OMPflag) #-O2


FC = gfortran

CC = g++

all : $(BINTEST)

$(BIN1) : $(OBJS) $(OBJS_Rad) $(BIN1).o
	$(FC) $^ $(LIBS) -o $@

$(BIN2) : $(OBJS) $(OBJS_Rad) $(BIN2).o
	$(FC) $^ $(LIBS) -o $@

$(BIN3) : $(OBJS) $(OBJS_Rad) $(BIN3).o
	$(FC) $^ $(LIBS) -o $@

$(BIN4) : $(OBJS) $(OBJS_Rad) $(BIN4).o
	$(FC) $^ $(LIBS) -o $@

$(BIN5) : $(OBJS) $(OBJS_Rad) $(BIN5).o
	$(FC) $^ $(LIBS) -o $@

$(BIN6) : $(OBJS) $(OBJS_Rad) $(BIN6).o
	$(FC) $^ $(LIBS) -o $@

$(BINTEST) : $(OBJS) $(OBJS_Rad) $(BINTEST).o
	$(FC) $^ $(LIBS) -o $@

$(VPATH_Rad)/%.o : $(VPATH_Rad)/%.f
	$(FC) $(fflags) -c $< -o $@

%.o : %.cpp
	$(CC) $(cflags) -c $< -o $@

clean :
	rm -f ${BINTEST} ${BINTEST}.o $(OBJS) $(OBJS_Rad)

