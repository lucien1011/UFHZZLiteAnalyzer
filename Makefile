CXX=c++
CC=gcc
CFLAGS= -std=c++14 -O2 -Wall
INSS=-I./include
INCMSSW=-I${CMSSW_BASE}/src

CFLAGS += `root-config --cflags`
LIBS += `root-config --ldflags --glibs`
LIBS += -L$(ROOFITSYS)/lib -lRooFit -lRooFitCore -lRooStats
LIBS += -L${CMSSW_BASE}/lib/slc6_amd64_gcc630/ -lZZMatrixElementMEKD -lZZMatrixElementMELA -lZZMatrixElementMEMCalculators -lZZMatrixElementPythonWrapper

zz4lOBJ=DarkZ_Ana.o
zxOBJ=ZX_Ana.o
zxDrawOBJ=ZX_Draw.o
upsilonAnaOBJ=Upsilon_Ana.o
kinZfitterOBJ=KinZfitter/KinZfitter.o
helperFunctionOBJ=KinZfitter/HelperFunction.o

.PHONY: clean all main test

all: DarkZ_Ana ZX_Ana ZX_Draw Upsilon_Ana

DarkZ_Ana: DarkZ_Ana.o 
		$(CXX) -o DarkZ_Ana.exe $(zz4lOBJ) $(helperFunctionOBJ) $(kinZfitterOBJ) $(LIBS)

ZX_Ana: ZX_Ana.o
		$(CXX) -o ZX_Ana.exe $(zxOBJ) $(helperFunctionOBJ) $(kinZfitterOBJ) $(LIBS)

ZX_Draw: ZX_Draw.o
		$(CXX) -o ZX_Draw.exe $(zxDrawOBJ) $(helperFunctionOBJ) $(kinZfitterOBJ) $(LIBS)

Upsilon_Ana: Upsilon_Ana.o 
		$(CXX) -o Upsilon_Ana.exe $(upsilonAnaOBJ) $(helperFunctionOBJ) $(kinZfitterOBJ) $(LIBS)

clean:
	@rm *.o *.exe

##############RULES##############
.cc.o:
	$(CXX) $(CFLAGS) $(INSS) ${INCMSSW} -I. -c $<
.cpp.o:
	$(CXX) $(CFLAGS) $(INSS) ${INCMSSW} -I. -c $<

