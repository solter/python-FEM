#generates all the plots for problem 2
import FEM

print "calculating fin401A"
FEM.FEMcalc("setups/fin401A.setup",True)
print "calculating fin401B"
FEM.FEMcalc("setups/fin401B.setup",True)
print "calculating fin801A"
FEM.FEMcalc("setups/fin801A.setup",True)
print "calculating fin801B"
FEM.FEMcalc("setups/fin801B.setup",True)
print "calculating fin801A2"
FEM.FEMcalc("setups/fin801A2.setup",True)
print "calculating fin801B2"
FEM.FEMcalc("setups/fin801B2.setup",True)

#generates all the plots for problem 2
import FVM

print "calculating fin2A1"
FVM.FVMcalc("setups/fin2A1.setup",True)
print "calculating fin2A2"
FVM.FVMcalc("setups/fin2A2.setup",True)
print "calculating fin2B1"
FVM.FVMcalc("setups/fin2B1.setup",True)
print "calculating fin2B2"
FVM.FVMcalc("setups/fin2B2.setup",True)
print "calculating fin21"
FVM.FVMcalc("setups/fin21.setup",True)
print "calculating fin22"
FVM.FVMcalc("setups/fin22.setup",True)
print "calculating fin2B_t100"
FVM.FVMcalc("setups/fin2B_t100.setup",True)
print "calculating fin2B2_t100"
FVM.FVMcalc("setups/fin2B2_t100.setup",True)
