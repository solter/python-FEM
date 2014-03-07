from scipy.integrate import quad
import math
import FEM
import matplotlib.pyplot as plt

def l2(f1,f2):
  tmp = quad(lambda x: (f1((x,)) - f2((x,)))*(f1((x,)) - f2((x,))), 0, 1)
  return math.sqrt(tmp)

prob140 = FEM.FEMcalc("prob140.setup")
prob180 = FEM.FEMcalc("prob180.setup")
prob240 = FEM.FEMcalc("prob240.setup")

ex1sol = lambda x: math.sin(x[0]) - math.sin(1)*x[0]
err40 = l2(ex1sol,prob140)
err80 = l2(ex1sol,prob180)

plt.figure()
prob140.pltSoln([ex1sol,'exact soln'])
plt.title('40 grid points')
plt.savefig('p1_40.png')

plt.figure()
prob140.pltSoln([ex1sol,'exact soln'])
plt.title('80 grid points')
plt.savefig('p1_80.png')


plt.figure()
prob140.pltSoln([ex1sol,'exact soln'])
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


