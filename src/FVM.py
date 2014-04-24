import numpy as np
import numpy.linalg as la
import scipy as sp
import scipy.sparse as sparse
import scipy.integrate as integrate
import scipy.sparse.linalg as spla
import string
import importlib
import math
import matplotlib.pyplot as plt
import time
from FEM import mesh
from FEM import simQuad
from FEM import polynom
from FEM import timeStep

class FVMcalc(object):
  
  def __init__(self, fname="FVM.setup"):
    self.domain = None #mesh of domain
    self.bdry = {0:None} #dictionary defining what to do with B.V.'s
    self.soln = None #The average cell values of the solution.
    self.recOrder = 1 #The reconstruction order
    self.numFlux = 1 #The numerical flux to use. See numFlux
    self.timeOrder = 1 #The accuracy of the time integration.
    self.initCond = None #the initial condition
    self.funcName = None #the name of the function file for the init condition
    self.endTime = 1. #the time to run simulation to
    self.outputFile = None

    self.parseInputFile(fname)
    self.getStarted()
    self.solnSteps()

  def parseInputFile(self, fname):
    """See example .setup file for how this parses
    stuff -> note that lines beginning with # are comments"""
    with open(fname,'r') as fl:
      secType = -1 #section types are 
      #0-MESH:, 1-FVM_FORM:, 
      for line in fl:#for every line in the file
        rdln = line.split()
        if(len(rdln) == 0):
          #skip empty lines
          continue
        if(rdln[0][0] == "#"):
          #skip comments
          continue
        
        #check for section type
        if(rdln[0] == "MESH:"):
          secType = 0
          continue
        if(rdln[0] == "FVM_FORM:"):
          secType = 1
          continue

        #parse the sections
        if(secType == 0):#MESH:
          if(rdln[0] == "file="):
            self.domain = mesh(rdln[1])

          elif(rdln[0].isdigit()):
            self.bdry[int(rdln[0])] = rdln[1]

        elif(secType ==1):#FVM_FORM:
          if(rdln[0] == "recon_order="):
            try:
              self.recOrder = int(rdln[1])
            except ValueError:
              raise NameError("reconstruction order must be an integer")
          elif(rdln[0] == "numFlux="):
            if(rdln[1] == "GLF"):
              self.numFlux = 1
            elif(rdln[1] == "Godunov"):
              self.numFlux = 0
            else:
              raise NameError(rdln[1] + " numFlux not supported")
          
          elif(rdln[0] == "timeOrder=")
            try:
              self.recOrder = int(rdln[1])
            except ValueError:
              raise NameError("time order must be an integer")

          elif(rdln[0] == "initCondition="):
            funcmod = importlib.import_module(rdln[1])
            funcmod = reload(funcmod)
            self.funcName = rdln[1] + " , " + rdln[2]
            self.initCond = getattr(funcmod, rdln[2])
          
          elif(rdln[0] == "runTo="):
            try:
              self.endTime = float(rdln[1])
            except ValueError:
              raise NameError("Endtime must be an integer")

  def getStarted(self):
    #initialize solution to empty list
    self.soln = np.zeros(len(self.domain.poly))

    #generate first step
    for i in range(len(self.domain.poly)):
      self.soln = simpQuad(initCond,
        self.domain.verts[self.domain.poly[i][0]],
        self.domain.verts[self.domain.poly[i][1]],
        )

    #start recording to output file
    self.outputFile = open(
      time.strftime("out/%Y_%M_%d_%H_%M_FVM.dat",time.localtime()),
      'w')
    self.outputFile.write(
      "#FVM for Berger's eqn output with the following parameters:\n
      #Reconstruction order = %d\n
      #Numerical Flux Scheme %d (0 - godunov, 1 - GLF)\n
      #Time integration order %d\n
      #Initial condition %s\n,
      #Endtime = %f\n#\nxs: " 
      % (self.recOrder, self.numFlux, self.timeOrder, 
      self.funcName, self.endTime)
    )

    for i in range(len(self.domain.poly))
      xmid = 
        .5*self.domain.verts[self.domain.poly[i][0]] +
        .5*self.domain.verts[self.domain.poly[i][1]]
      self.outputFile.write(" %f "%xmid)
   
    self.outputFile.write("\n\nt\t\t| vals\n%f | %s"%(0.0, self.soln))

  def solnSteps(self):
    """Calculate steps and store them to output file
    """
    t = 0.0
    while(t < self.endTime):#step through problem until final solution
      #calculate dt, assume uniform partition
      dt = self.domain.verts[self.poly[0][1]] - 
        self.domain.verts[self.poly[0][0]]
      dt /= max(abs(self.soln))
      if(self.timeorder > 1):
        #for stability in more accurate problems
        dt *= .5
      self.soln = timeStep(self.soln, updateF, self.timeorder, dt)
      t += dt
      self.outputFile.write(
        "\n\nt\t\t| vals\n%f | %s"%(t, self.soln))

    self.outputFile.close();

  def ENO(self, j, soln):
    """Does an ENO reconstruction
    j -> the initial cell in the stencil
    soln -> u vector to reconstruct values with

    returns a list of length 2, the first
    element is the reconstruction at j-1/2,
    the right second is the reconstruction at j+1/2

    Note that this currently assumes periodic boundary conditions,
    and that the mesh intervals are ordered
    """

    #the trivial case of first order
    if(self.recOrder == 1):
      return [soln[stencil[0]], soln[stencil[0]]]

    #initialize stencil
    stMin = j
    stMax = j
    stencil = [j]
    n = len(soln)
    
    #form stencil
    while(len(stencil) < self.recOrder):
      if( 
        diff(stencil + [(stMax + 1)%n]) >
        diff(stencil + [(stMin - 1)%n])
      ):#if right element larger difference
        stencil += [(stMin - 1)%n]
        stMin -= 1
      else:
        stencil == [(stMax + 1)%n]
        stMax += 1
    
    #calculate polynomial, via Ac = b, where c is the coeficient list
    A = np.zeros((self.recOrder,self.recOrder),dtype=float)
    b = np.zeros(self.recOrder)
    
    for i in range(self.recOrder):
      xl = self.domain.verts[self.domain.poly[stencil[i]][0]]
      xr = self.domain.verts[self.domain.poly[stencil[i]][1]]
      #check if this needs to wrap around:
      if(stencil[i] - stencil[0] > self.recOrder - 1):#push back
        xl -= self.domain.verts[self.domain.poly[0][0]]
        xr -= self.domain.verts[self.domain.poly[0][0]]
      elif(stencil[0] - stencil[i] > self.recOrder - 1):#push forward
        xl += self.domain.verts[self.domain.poly[-1][1]]
        xr += self.domain.verts[self.domain.poly[-1][1]]

      b[i] = soln[stencil[i]] * (xr - xl)
      for j in range(self.recOrder):
        A[i,j] = (xr**(j+1) - xl**(j+1))/(j+1)
    
    coefMat = la.solve(A,b)
    xl = self.domain.verts[self.domain.poly[stencil[0]][0]]
    xr = self.domain.verts[self.domain.poly[stencil[0]][1]]
    ul = coefMat * [[xl**i] for i in range(self.recOrder)]
    ur = coefMat * [[xr**i] for i in range(self.recOrder)]
    return [ul, ur]

  def diff(self, sten):
    """Computes D[stencil] without scaling by interval length
    """
    if(len(sten) == 1):
      return self.soln[sten[0]]
    else:
      return diff(sten[1:]) - diff(sten[:-1]) 

  def numFlux(self, ul, ur, alpha = 0):
    """The numerical flux.
    self.numFlux = 0 -> Gudonov scheme
    self.numFlux = 1 -> Lax-Friedrich scheme

    ul -> the left value
    ur -> the right value
    alpha -> only used for LF, should be the max|f'(u)|,
      which for Bergers equation is max|u|
    """
    if(self.numFlux == 0):
      if(ul < ur):
        if(ul > 0):
          return .5*ul*ul
        elif(ur < 0):
          return .5*ur*ur
        else:#if ul < 0 and ur > 0
          return 0
      elif(ur > ul):
        if(ur + ul > 0):
          return .5*ul*ul
        else:
          return .5*ur*ur
      else:
        return .5*ur*ul
    elif(self.numFlux == 1):
      return .5* (
        .5*ul**2 + .5*ur**2 - alpha*(ur - ul)
      )

    else:
      raise NameError("your flux is not supported")

  def updateF(self,u):
    """this calculates u_t = F(u)
    """

    #get info at far left
    u = np.array(u)
    unew = np.zeros(u.size)
    alpha = max(abs(u))
    u0 = ENO(0,u)
    up = u0[:]
    v0 = [self.domain.verts[self.domain.polys[0][0]],
      self.domain.verts[self.domain.polys[0][1]]]
    vp = v0[:]

    #for every interval
    for i in range(1,len(u)):
      #get subsequent interval info
      un = ENO(i,u)
      vn = [self.domain.verts[self.domain.polys[i][0]],
        self.domain.verts[self.domain.polys[i][1]]]
      #calculate numeric flux between i and i-1
      flx = numFlux(up[1],un[0],alpha)
      #update new u vals
      unew[i] += flx/(vn[1] - vn[0])
      unew[i-1] -= flx/(vp[1] - vp[0])
      #move new variables into previous variables slot
      up = un[:]
      vp = vn[:]

    #take care of interface between first and last value
    #calculate numeric flux between n and 0
    flx = numFlux(up[1],u0[0],alpha)
    #update new u vals
    unew[0] += flx/(v0[1] - v0[0])
    unew[-1] -= flx/(vp[1] - vp[0])
     
    return unew

#TODO: Plot solutions
