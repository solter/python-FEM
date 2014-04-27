#generates all the plots for problem 2
import FEM

FEM.FEMcalc("setups/fin401A.setup",True)
FEM.FEMcalc("setups/fin401B.setup",True)
FEM.FEMcalc("setups/fin801A.setup",True)
FEM.FEMcalc("setups/fin801B.setup",True)
FEM.FEMcalc("setups/fin801A2.setup",True)
FEM.FEMcalc("setups/fin801B2.setup",True)

#generates all the plots for problem 2
import FVM

FVM.FVMcalc("setups/fin2A1.setup",True)
FVM.FVMcalc("setups/fin2A2.setup",True)
FVM.FVMcalc("setups/fin2B1.setup",True)
FVM.FVMcalc("setups/fin2B2.setup",True)
FVM.FVMcalc("setups/fin21.setup",True)
FVM.FVMcalc("setups/fin22.setup",True)
FVM.FVMcalc("setups/fin2B_t100.setup",True)
FVM.FVMcalc("setups/fin2B2_t100.setup",True)
