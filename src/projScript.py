from scipy.integrate import quad
from scipy.integrate import dblquad
import numpy as np
import math
import FEM
import matplotlib.pyplot as plt
#for animating the 3d plots
import sys
sys.path.append('/home/solter/Documents/Scripts')
import animPlot as anim

#define the l2 error as || f1 - f2 ||_L2
def l2d(f1,f2,dim):
  tmp = float('nan') 
  if(dim == 1):
    tmp = quad(lambda x: (f1((x,)) - f2((x,)))*(f1((x,)) - f2((x,))), 0, 1,epsrel=5e-3)[0]
  if(dim == 2):
    return .1
    tmp = dblquad(lambda x,y: (f1((x,y)) - f2((x,y)))*(f1((x,y)) - f2((x,y))), 0, 1, lambda x: 0, lambda y: 1,epsrel=5e-3)[0]
  return math.sqrt(tmp)

#go through and calculate all the FE solutions
femsol = [[],[]]
meshSz = [[],[]]
poly = [[],[]]
meshSz[0] = [5*2**i for i in range(5)]
meshSz[1] = [5*2**i for i in range(3,4)]
#meshSz[1] = [5*2**i for i in range(4)]
poly[0] = range(1,4)
poly[1] = range(1,3)
#poly[1] = range(1,3)
for probNum in range(len(femsol)):
  for meshNum in range(len(meshSz[probNum])):
    femsol[probNum].append([])
    for polyNum in range(len(poly[probNum])):
      namestr = "prob%d_%d_%d.setup"%(probNum+1,meshSz[probNum][meshNum],poly[probNum][polyNum])
      print "Problem " + namestr + " ..."
      femsol[probNum][meshNum].append(
        FEM.FEMcalc("setups/" + namestr)
      )
      print " solved"

print "calculating l2 errors"
#go through and calculate l2 errors and convergence rates
exsol = [None,None]
exsol[0] = lambda x: math.sin(x[0]) - math.sin(1)*x[0]
exsol[1] = lambda x,y: femsol[1][-1][-1]#use best solution as exact solution
l2err = [[],[]]
conRate = [[],[]]
for probNum in range(len(femsol)):
  for meshNum in range(len(meshSz[probNum])):
    l2err[probNum].append([])
    conRate[probNum].append([])
    for polyNum in range(len(poly[probNum])):
      l2err[probNum][meshNum].append(
        l2d(femsol[probNum][meshNum][polyNum],exsol[probNum],probNum+1)
      )
      if(meshNum > 0):
        conRate[probNum][meshNum].append(
          math.log(
            l2err[probNum][meshNum-1][polyNum] / l2err[probNum][meshNum][polyNum]
          ,2)
        )
      else:
        conRate[probNum][meshNum].append([float('nan')])

print "printing plots"
#create plots -> the 3d plots will be gif images which
#rotate the figure for better viewing
angles = []
fignames = ["",""]
for probNum in range(len(femsol)):
  if(probNum == 1):
      angles = np.linspace(0,360,21)[:-1] # Take 20 angles between 0 and 360

  for meshNum in range(len(meshSz[probNum])):
    for polyNum in range(len(poly[probNum])):
      ax = FEM.pltSoln(femsol[probNum][meshNum][polyNum])
      nameStr = "outputs/prob%d_%d_%d"%(probNum+1,meshSz[probNum][meshNum],poly[probNum][polyNum])
      #save the plots
      if(probNum == 0):
        plt.savefig(nameStr + ".png") 
        fignames[0] += nameStr + ".png\n"
      elif(probNum == 1):
        plt.savefig(nameStr + ".png") 
        anim.rotanimate(ax, angles,nameStr + '.gif',delay=20)
        fignames[0] += nameStr + ".png,  " + nameStr + ".gif\n" 
print 

#Writing output
with open("outputs/proj1Results.txt",'w') as f:
  f.write("MA5629\nProject 1\nPeter Solfest\n\n")
  f.write("=================================\n\n")
  f.write("h represents the following:\n")
  f.write("1D -> the interval length\n")
  f.write("2D -> the smallest triangle edge length\n\n")
  f.write("<#> l2 err represents the l2 error of the solution\
  using polynomial basis of degree <#>\n")
  f.write("CR represents the rate of convergence\n")
  f.write("For the source code, see the following github repo:\n")
  f.write("https://github.com/solter/python-FEM\n")
  

  for probNum in range(len(femsol)):
    f.write(("\n\nProblem %d\n=============\n\n"%(probNum+1)))
    if(probNum == 0):
      f.write(" M  |    h     |")
      for polyNum in range(len(poly[probNum])):
        f.write((" %d l2 err |   CR    |"%poly[probNum][polyNum]))
      f.write("\n-----------")
      for polyNum in range(len(poly[probNum])):
        f.write("---------------------------")

      for meshNum in range(len(meshSz[probNum])):
        h = 0
        if(probNum == 0):
          h = 1./meshSz[probNum][meshNum]
        elif(probNum == 1):
          h = 1./(2+meshSz[probNum][meshNum])
          h *= math.sqrt(2)
        
        f.write(("\n %2d | %1.2e |"%(meshSz[probNum][meshNum]h)))
        for polyNum in range(len(poly[probNum])):
          f.write((" %1.2e |"%l2err[probNum][meshNum][polyNum]))
          if(meshNum == 0):
            f.write("         |")
          else:
            f.write(("   %.2f   |"%conRate[probNum][meshNum][polyNum]))
    
    f.write("The following are images of the solutions:")
    f.write(fignames[probNum])
