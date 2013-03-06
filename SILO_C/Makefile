MPPL =      $(HOME)/hym/MPP/MPPL
INFLAGS =   -r8 -o
DBG1 =       -g -Mbounds -Minfo -v
DBG_FLAGS =  $(DBG1) $(INFLAGS)
#OPT_FLAGS =  -O3 -mcpu=opteron -mtune=opteron
OPT_FLAGS =  -fast
HYBM_FLAGS = $(OPT_FLAGS) $(INFLAGS)
MPI_f90 = ftn
F90 = pgf90

#==============================================================================#

CXX = pgCC

#==============================================================================#

test: 
	$(CXX) -o tester.exe src_drivers/tester.cpp

#==============================================================================#
        
SRCDRV = src_drivers
SRCPKG = src_package
PDIR = /p/hym/cmyers
SILO = $(SILO_LIB) -lsilo -I$(SILO_INC)
INC = $(SILO) -I$(SRCPKG) -I$(SRCDRV)

BF  = Basic_Functions
AR  = ASCII_Read
AW  = ASCII_Write
SR  = SILO_Read
SW  = SILO_Write
HDO = HYM_DataObj
SCO = SILO_CycObj
InterF  = Interp_Functions
IntegF  = Integ_Functions
SOBJ = $(BF).o $(AW).o $(HDO).o $(SCO).o
POBJ = $(BF).o $(SR).o $(InterF).o $(IntegF).o

F3D = HYM_SILO
FP  = Probe_SILO
SMD = SILO_Mode_Data
RCC = SILO_mode_data_RCC
SM2 = SILO_mode_data_v2
GJV = Get_Jmax_vmax_n0
GFH = Get_FR_height_P

$(BF).o: $(SRCPKG)/$(BF).cpp 
	$(CXX) $(INC) -c $(SRCPKG)/$(BF).cpp
$(AW).o: $(SRCPKG)/$(AW).cpp
	$(CXX) $(INC) -c $(SRCPKG)/$(AW).cpp
$(SR).o: $(SRCPKG)/$(SR).cpp
	$(CXX) $(INC) -c $(SRCPKG)/$(SR).cpp
$(SW).o: $(SRCPKG)/$(SW).cpp
	$(CXX) $(INC) -c $(SRCPKG)/$(SW).cpp
$(HDO).o: $(SRCPKG)/$(HDO).cpp
	$(CXX) $(INC) -c $(SRCPKG)/$(HDO).cpp
$(SCO).o: $(SRCPKG)/$(SCO).cpp
	$(CXX) $(INC) -c $(SRCPKG)/$(SCO).cpp
$(InterF).o: $(SRCPKG)/$(InterF).cpp
	$(CXX) $(INC) -c $(SRCPKG)/$(InterF).cpp
$(IntegF).o: $(SRCPKG)/$(IntegF).cpp
	$(CXX) $(INC) -c $(SRCPKG)/$(IntegF).cpp

$(F3D).o: $(SRCDRV)/$(F3D).hpp $(SRCDRV)/$(F3D).cpp
	$(CXX) $(INC) -c $(SRCDRV)/$(F3D).cpp
$(FP).o: $(SRCDRV)/$(FP).cpp
	$(CXX) $(INC) -c $(SRCDRV)/$(FP).cpp
$(RCC).o: $(SRCDRV)/$(RCC).cpp
	$(CXX) $(INC) -c $(SRCDRV)/$(RCC).cpp
$(SM2).o: $(SRCDRV)/$(SM2).cpp
	$(CXX) $(INC) -c $(SRCDRV)/$(SM2).cpp

silo:  $(SOBJ) $(F3D).o
	$(CXX) -g -O0 $(F3D).o $(SOBJ) $(INC) -o ../$(F3D).exe

probe: $(POBJ) $(FP).o
	$(CXX) $(FP).o $(POBJ) $(INC) -o $(FP).exe
RCC: $(POBJ) $(RCC).o
	$(CXX) $(RCC).o $(POBJ) $(INC) -o $(RCC).exe
#	./$(RCC).exe $(PDIR)/RunData_CounterH/2011_07_25_ct_RCC/SILO/ $(PDIR)/RunData_CounterH/2011_07_25_ct_RCC/Extracted_Data/ 265
modes2: $(POBJ)
	$(CXX) $(POBJ) $(INC) -o $(SM2).exe $(SRCDRV)/$(SM2).cpp
#	./$(SM2).exe $(PDIR)/RunData_CounterH/2011_01_21_ct_HR/SILO/ $(PDIR)/RunData_CounterH/2011_01_21_ct_HR/SILO_mode_data_profiles/ 0
getmax: $(POBJ)
	$(CXX) $(POBJ) $(INC) -o $(GJV)_w_mins.exe $(SRCDRV)/$(GJV).cpp
getheight: $(POBJ)
	$(CXX) $(POBJ) $(INC) -o $(GFH).exe $(SRCDRV)/$(GFH).cpp

clean:
	rm *.o *.exe

