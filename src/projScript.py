from scipy.integrate import quad
import math
import FEM
import matplotlib.pyplot as plt

#define the l2 error as || f1 - f2 ||_L2
def l2(f1,f2):
  tmp = quad(lambda x: (f1((x,)) - f2((x,)))*(f1((x,)) - f2((x,))), 0, 1)
  return math.sqrt(tmp)

#go through and calculate all the FE solutions
fems = []
for probNum in range(1,3):
  fems.append([])
  for meshSz in range(0,5):
    if(probNum == 2 && meshSz == 3):
      continue
    fems[probNum-1].append([])
    for poly in range(1,4):
      if(probNum == 2 && poly == 3):
        continue
      fems[probNum-1][meshSz].append(
        FEM.FEMcalc("setups/prob%d_%d_%d.setup"%(probNum,5*2**meshSz,poly))
      )

#TODO: implement accuracy & convergence tests

ex1sol = lambda x: math.sin(x[0]) - math.sin(1)*x[0]
err10 = l2(ex1sol,prob140)
err20 = l2(ex1sol,prob140)
err40 = l2(ex1sol,prob140)
err80 = l2(ex1sol,prob180)

plt.figure()
pltSoln(prob120,[ex1sol,'exact soln'])
plt.title('40 grid points')
plt.savefig('p1_40.png')

plt.figure()
pltSoln(prob140,[ex1sol,'exact soln'])
plt.title('80 grid points')
plt.savefig('p1_80.png')


plt.figure()
pltSoln(prob180,[ex1sol,'exact soln'])
plt.title('40x40 grid points')
plt.savefig('p2_40.png')


with open("proj1Outs.txt",'w') as f:
  s = "Project 1 problems\n1\n=========\n"
  f.write(s)
  s = "l2 error for 40 mesh points:\n%f\n\n"%err40
  f.write(s)
  s = "l2 error for 80 mesh points:\n%f\n\n"%err80
  f.write(s)
  s = "order accuracy: %f\n\n"%math.log(err40/err80,2)
  f.write(s)
  s = "Figures saved as p1_40.png, p1_80.png and p2_40.png"
  f.write(s)
