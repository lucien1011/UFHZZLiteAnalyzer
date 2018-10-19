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
zxWeightOBJ=ZX_Weight.o
upsilonAnaOBJ=Upsilon_Ana.o
jpsiAnaOBJ=Jpsi_Ana.o
wfcAnaOBJ=WFC_Ana.o
kinZfitterOBJ=KinZfitter/KinZfitter.o
helperFunctionOBJ=KinZfitter/HelperFunction.o

.PHONY: clean all main test

all: DarkZ_Ana ZX_Ana ZX_Weight Upsilon_Ana Jpsi_Ana WFC_Ana

DarkZ_Ana: DarkZ_Ana.o 
		$(CXX) -o DarkZ_Ana.exe $(zz4lOBJ) $(helperFunctionOBJ) $(kinZfitterOBJ) $(LIBS)

ZX_Ana: ZX_Ana.o
		$(CXX) -o ZX_Ana.exe $(zxOBJ) $(helperFunctionOBJ) $(kinZfitterOBJ) $(LIBS)

ZX_Weight: ZX_Weight.o
		$(CXX) -o ZX_Weight.exe $(zxWeightOBJ) $(helperFunctionOBJ) $(kinZfitterOBJ) $(LIBS)

Upsilon_Ana: Upsilon_Ana.o 
		$(CXX) -o Upsilon_Ana.exe $(upsilonAnaOBJ) $(helperFunctionOBJ) $(kinZfitterOBJ) $(LIBS)

Jpsi_Ana: Jpsi_Ana.o 
		$(CXX) -o Jpsi_Ana.exe $(jpsiAnaOBJ) $(helperFunctionOBJ) $(kinZfitterOBJ) $(LIBS)

WFC_Ana: WFC_Ana.o 
		$(CXX) -o WFC_Ana.exe $(wfcAnaOBJ) $(helperFunctionOBJ) $(kinZfitterOBJ) $(LIBS)

clean:
	@rm *.o *.exe

##############RULES##############
.cc.o:
	$(CXX) $(CFLAGS) $(INSS) ${INCMSSW} -I. -c $<
.cpp.o:
	$(CXX) $(CFLAGS) $(INSS) ${INCMSSW} -I. -c $<

