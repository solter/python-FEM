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
from mpl_toolkits.mplot3d import Axes3D


class polynom(object):
  """A class to deal with 1D and 2D polynomials.
  """

  def __init__(self, dim, deg = 1, coef = np.array([])):
    """Initializes a polynomial
    
    PARAMETERS: 
    dim : the dimension of the space (>2 may be buggy)
    deg : the expected degree of the polynomial
    coef : the coeffiecient matrix. If empty, will be initialized to zeros
    """
    self.dim = dim
    self.deg = deg
    self.coef = np.zeros((self.deg+1)*np.ones(self.dim),dtype=float)
    if(coef.size > 0):
      if(coef.ndim == dim):
        self.coef = coef
      else:
        self.setFlatCoef(coef)
        
  def addElem(self, idx, coef):
    """Add nonzero polynomial element
    
    PARAMETERS:
    idx : index of parameter. Must be a tuple of same length
      as dimensions of polynomial.
    coef : value of the coeficient

    example:
    for a 2-D polynomial, add the element 10*x1*x2^2:
    addElem((1, 2), 10)
    """
   
    self.deg = max([sum(idx), self.deg])
    self.dim = self.dim
    
    #check that the appropriate number of dimensions is given
    if(len(idx) != self.dim):
      raise NameError("Index and polynomial dimension mismatch")

    #expand the coef matrix if needed
    if(len(self.coef) <= self.deg):
      pred = len(self.coef)
      padIdx = [pred]*self.dim
      for i in range(self.dim):
        if(i > 0):
          padIdx[i-1] = self.deg + 1
        padIdx[i] = self.deg + 1 - pred
        self.coef = np.append(self.coef,np.zeros(padIdx),i)
  
    self.coef[idx] += coef

  def __str__(self):
    """Generates a string representation"""
    toret =  self.__genStr((),"")
    if(len(toret) > 2):
      return toret[:-2]
    else:
      return "0"

  def __genStr(self, idx, toret):
    """Helper function for __str__"""
    if self.dim == len(idx):
      #if coef at idx is 0, don't print it
      if(self.coef[idx] == 0):
        return toret
      #if coef at idx isn't 0, add it to return string
      else:
        toret += str(self.coef[idx])
        for i in range(len(idx)):
          if (idx[i] != 0):
            toret += "*x%d^%d" %(i, idx[i])
        return toret + " + "
    else:#recurse
      for i in range(self.coef.shape[len(idx)]):
        tmp = idx + (i,)
        toret = self.__genStr(tmp, toret);
      return toret

  def __add__(p1, p2): 
    #transform scalar into polynomial
    if(np.isscalar(p2)):
      tmp = p2
      p2 = polynom(self.dim,0)
      p2.addElem(tuple(np.zeros(self.dim,dtype=int)),tmp)
  
    dim = p1.dim
    if(dim != p2.dim):
      raise NameError("polynomial dimension mismatch")
    if(p1.deg < p2.deg):
      tmp = p1
      p1 = p2
      p2 = tmp

    #make temporary polynomial which has p2 data
    #and add a 0 element of the right degree so they can be added
    toret = polynom(p1.dim,p1.deg,p1.coef.copy())
    if(toret.dim == 1):
      toret.coef[0:p2.deg+1] += p2.coef
    elif(toret.dim == 2):
      toret.coef[0:p2.deg+1,0:p2.deg+1] += p2.coef
    else:
      raise NameError("%d-d Polynomials not supported"%self.dim)

    return toret

  #variations on addition
  def __radd__(self,other):
    return self.__add__(self,other)

  def __sub__(self,p2):
    return self.__add__(-1*p2)

  def __rsub__(self,p2):
    self = -1*self
    return self + p2

  def __mul__(p1, p2):
    #transform scalar into polynomial
    if(np.isscalar(p2)):
      tmp = p2
      p2 = polynom(p1.dim,0)
      p2.addElem(tuple(np.zeros(p1.dim,dtype=int)),tmp)
   
    #initialize polynomial
    toret = polynom(p1.dim, p2.deg + p1.deg)
   
    #make sure p1 has larger degree
    if(p1.deg < p2.deg):
      tmp = p2
      p2 = p1
      p1 = tmp

    if(toret.dim == 1):
      nzs = zip(p2.coef.nonzero()[0])
      for idx in nzs:
        toret.coef[idx[0]:idx[0] + p1.deg + 1] += p2.coef[idx] * p1.coef  
    elif(toret.dim == 2):
      nzs = zip(p2.coef.nonzero()[0],p2.coef.nonzero()[1])   
      for idx in nzs:
        toret.coef[idx[0]: idx[0] + p1.deg + 1, idx[1]: idx[1] + p1.deg + 1] += p2.coef[idx] * p1.coef  
    else:
      raise NameError("%d-d Polynomials not supported"%self.dim)
    
    return toret

  #variation on multiplication
  def __rmul__(p1,p2):
    return p1.__mul__(p2)

  #copy to not get have same stuff referenced
  def copy(self):
    return polynom(self.dim, self.deg, self.coef.copy())

  def diff(self, var = 0, norm = True):
    """Generates the derivative w.r.t. x<var> of the polynomial
    
    Parameters:
    var : the variable to take the derivative w.r.t
    norm : If False this will find the antiderivative
    """
    
    if(self.dim == 1):
      if(var != 0):
        raise NameError("Only valid variable is x0")
        
      if(norm):#if a derivative
        toret = polynom(self.dim,self.deg - 1, np.delete(self.coef,0,0))
        for n in range(toret.deg+1):
          toret.coef[n] *= n+1
      else:#the antiderivative
        toret = polynom(self.dim, self.deg + 1, np.insert(self.coef,0,np.zeros(1),0))
        for n in range(1,toret.deg+1):
          toret.coef[n] /= n
      
    elif(self.dim == 2):
      if(not (var == 0 or var == 1)):
        raise NameError("Only valid variables are x0 and x1")
      
      if(norm):#if a derivative
        toret = polynom(self.dim,self.deg - 1, np.delete(np.delete(self.coef,0,var),-1,1 - var) ) 
        for n in range(toret.deg+1):
          if(var==0):
            toret.coef[n,:] *= n+1
          else:
            toret.coef[:,n] *= n+1
      else:#the antiderivate
        if(var == 0):
          idx1 = 0
          idx2 = self.deg+1
        else:
          idx1 = self.deg+1
          idx2 = 0
        tmpCoef =  np.insert(self.coef,idx1, np.zeros(self.deg+1),0)
        tmpCoef = np.insert(tmpCoef,idx2,np.zeros(self.deg+2),1)
        toret = polynom(self.dim, self.deg + 1, tmpCoef)
        for n in range(1,toret.deg+1):
          if(var==0):
            toret.coef[n,:] /= n
          else:
            toret.coef[:,n] /= n
    else:
      raise NameError("%d-d Polynomials not supported"%self.dim)

    return toret

  def __call__(self,x):
    """Evaluates the polynomial at x

    for 2-D:
    if one of the values in x is 't',
    this variable is left untouched and a 1d polynomial is returned

    Example:
      p1 = x0 + 2*x1
      p1((1,'t')) = 1 + 2*x
      p1(('t',1)) = x + 2
    """
    
    if(len(x) != self.dim):
      raise NameError("Input need the right dimensions")

    toret = 0
    if(self.dim==1):
      #for all nonzero elements
      for idx in zip(self.coef.nonzero()[0]):
        toret += self.coef[idx] * x[0]**idx[0]

    elif(self.dim==2):
      #for all nonzero elements
      if(x[0] == 't' or x[1] == 't'):
        toret = polynom(1,self.deg)

      for idx in zip(self.coef.nonzero()[0],self.coef.nonzero()[1]):
        if(x[0] == 't' and x[1] == 't'):
          toret.addElem((idx[0] + idx[1],),self.coef[idx])
        elif(x[0] == 't'):
          toret.addElem((idx[0],),self.coef[idx] * x[1]**idx[1])
        elif(x[1] == 't'):
          toret.addElem((idx[1],),self.coef[idx] * x[0]**idx[0])
        else:
          toret += self.coef[idx] * x[0]**idx[0] * x[1]**idx[1]
    else:
      raise NameError("%d-d polynomials not supported"%self.dim)
    
    return toret

  def integrate(self,x0,x1):
    """integrate polynomial from x0 to x1
    
    Parameters:
    x0,x1 -> tuples representing the endpoints of integration.
    elements must be either numbers or 't'
    if empty uses example 1.1 and 2.2
      Examples:
        1D: 
        1.1 integrate((0,),(1,)) = int_0^1 f dx (returns a number)
        1.2 integrate((0,),('t',)) = int_0^t f dx (returns a polynomial)

        2D:
        1.1 integrate((0,0),(1,2)) = int_0^1 int_0^2 f dx0 dx1 (returns a number)
        1.2 integrate((0,0),(1,'t')) = int_0^1 int_0^x1 f dx0 dx1 (returns a number)
        1.3 integrate((0,0),('t','t')) = int_0^t int_0^x1 f dx0 dx1 (returns a 1d polynomial)
    """
   
    if(len(x0) == 0):
      if(self.dim == 1):
        x0 = (0,)
        x1 = (1,)
      elif(self.dim == 2):
        x0 = (0,0)
        x1 = (1,'t')

    if(self.dim == 1):
      p = self.diff(0,False)
      if(x1[0] == x0[0]):
        return 0
      elif(x1[0] == 't'):
        return p - p((x0[0],))
      elif(x0[0] == 't'):
        return p((x1[0],)) - p
      else:
        return p((x1[0],)) - p((x0[0],))
    elif(self.dim == 2):
      if(x1[0] == x0[0] or x1[1] == x0[1]):
        return 0
      else:
        p1 = self.diff(0,False)
        p1 = p1((x1[1],'t')) - p1((x0[1],'t'))
        return p1.integrate((x0[0],),(x1[0],))
    else:
      raise NameError("%d-d Polynomials not supported"%self.dim)

  def getFlatCoef(self, x, diffOrder=0, diffVariable=()):
    """Gets the coeficients of polynomial when evaluated at x

    returns a 1-d array of values
    Parameters:
    x -> n-dimensional point to evaluate at
    diffOrder -> If this polynomial represents a diffOrder^th 
      derivative,the output array will be coefficients of original
      polynomial
    diffVariable -> if diffOrder != 0, this is a tuple of length
      poly.dim. Where diffVariable[i] represents the number of
      derivatives taken w.r.t. the i^th variable.
      By default all derivatives are assumed to be taken w.r.t
      the 0^th variable
    """
    
    if(len(diffVariable) == 0):
        diffVariable = diffVariable + (diffOrder,)
        for i in range(1,self.dim):
          diffVariable = diffVariable + (0,)
    
    if(diffOrder != sum(diffVariable)):
      raise NameError("diffVariable and diffOrder mismatch")

    if(len(diffVariable) != self.dim):
      raise NameError("diffVariable must be of dimension %d"%self.dim)

    if(len(x) != self.dim):
      raise NameError("x must be of dimension %d"%self.dim)

    if(self.dim==1):
      toret = np.zeros(self.deg+1 + diffOrder)
      #for all nonzero elements
      for idx in zip(self.coef.nonzero()[0]):
        toret[idx[0] + diffOrder] = self.coef[idx] * x[0]**idx[0]

    elif(self.dim==2):
      #zeros of length n*(n+1)/2, where n = self.deg + 1 + diffOrder
      toret = np.zeros(((self.deg + 1 + diffOrder)*(self.deg + 1 + diffOrder + 1))/2)
      for idx in zip(self.coef.nonzero()[0],self.coef.nonzero()[1]):
        diag = idx[0] + diffVariable[0] + idx[1] + diffVariable[1]
        flatIdx = sum(range(1,diag + 1)) + idx[1] + diffVariable[1]
        toret[flatIdx] += self.coef[idx] * x[0]**idx[0] * x[1]**idx[1]
    
    else:
      raise NameError("%d-d polynomials not supported"%self.dim)
    
    return toret

  def setFlatCoef(self,flatCoef):
    """Takes a 1-D array and uses it to establish coeficients

    NOTE: the polynomial must be initialized with the appropriate
    dimension and degree before calling this
    Parameters:
    coef -> the coefficients in a flat matrix.
    If coef is empty, will create a polynomial with 1's as all coefficients
    Examples:
    1D: [1 1 0 -4] <-> 1 + x -4x^3
    2D (deg 3): [1 1 0 -4 5 0 -.5 0 0 1] <-> 1 + x0 -4 x0^2 +5 x0 x1 -.5 x0^3 + x1^3
    """
    
    if(self.dim == 1):
      if(len(flatCoef) == 0):
        flatCoef = np.ones(self.deg + 1)
      
      if(len(flatCoef) != self.deg + 1):
        raise NameError("flatCoef incorrect length")

      self.coef = np.array(flatCoef)
      return

    elif(self.dim == 2):
      totlen = ((self.deg + 1)*(self.deg + 2))/2
      if(len(flatCoef) == 0):
        flatCoef = np.ones(totlen)

      if(len(flatCoef) != totlen):
        raise NameError("flatCoef incorrect length")

      for diag in range(self.deg+1):
        for col in range(diag+1):
          idx = (diag - col,col)
          flatIdx = sum(range(1,diag + 1)) + col
          self.coef[idx] = flatCoef[flatIdx]
      
    else:
      raise NameError("%d-d Polynomials not supported"%self.dim)

class mesh(object):
  """A mesh for the grid

  members:
  verts -> list of vertices. 
    verts[i] contains the coordinates for this vertex
  v_is_bndry -> v_is_bndry[i] = True if verts[i] is on the boundary
  poly -> list of elements. 
    poly[i] contains a list of verts which define its borders
  v_res -> list of elements which the vertex's are in.
    v_res[i] is a list of all elements which verts[i] is a part of
  """

  def __init__(self,fname):
    self.readMesh(fname)

  def readMesh(self,fname):
    """Read in nodes and values from a tetgen files
    
    This is set up to deal with 1-D and 2-D .node and .ele files
    """
    getFormat = True
    #parse the node file
    with open(fname + ".node",'r') as f:
      skip = 0
      for line in f:#for each line in the file
        #if not a comment
        if(line.split()[0][0] != "#"):
          if(getFormat):
            tmp = [int(s) for s in line.split()[0:5]]
            #generate matrix of empty 0's
            self.verts = np.zeros((tmp[0], tmp[1]))
            #make a vector of False
            self.v_is_bndry = [0 for i in range(tmp[0])]
            self.v_res = [[] for i in range(tmp[0])]
            skip = tmp[2]# num of attributes
            if(tmp[3] == 0):
              raise NameError("node file must have boundary markers")
            getFormat = False
          else:
            s = line.split()
            idx = int(s.pop(0)) - 1#the node number
            #populate the node values
            for i in range(self.verts.shape[1]):
              self.verts[(idx,i)] = float(s.pop(0))
            #if it's a boundary node, mark it as such
            if(int(s[skip]) != 0):
              self.v_is_bndry[idx] = int(s[skip])
    
    getFormat = True
    #parse the ele file
    with open(fname + ".ele") as f:
      for line in f:#for every line in the file
        #if not comment
        if(line.split()[0][0] != "#"):
          if(getFormat):
             tmp = [int(s) for s in line.split()[0:2]]
             #generate matrix of empty 0's
             self.poly = np.zeros((tmp[0],tmp[1]),dtype=int)
             getFormat = False
          else:
            s = line.split()
            idx = int(s.pop(0)) - 1#the poly number
            #populate the poly values
            for i in range(self.poly.shape[1]):
              nd = int(s.pop(0)) - 1
              self.poly[(idx,i)] = nd#add node to this element
              self.v_res[nd].append(idx)#add this element nodes residence list

  def refineMesh(self):
    """refines the mesh"""
    raise NameError("This function is currently broken")
    Overts = self.verts.copy()
    Ov_is_bndry = self.v_is_bndry.copy()
    Opoly = self.poly.copy()
    Ov_res = self.v_res.copy()

    self.verts = Overts.tolist()
    self.poly = np.zeros((2**self.dim * Opoly.shape[0],Opoly.shape[1]))
    self.v_res = []
    for i in len(Opoly):
      if(self.dim == 1):
        #THIS FAILS
        nidx = len(self.verts)
        self.vert.append((.5*Overts[Opoly[i][0]] + .5*Overts[Opoly[i][1]]).tolist())
        self.poly[2*i] = [Opoly[i][0], nidx]
        self.poly[2*i + 1] = [nidx, Opoly[i][1]]
      elif(self.dim == 2):
        #THIS FAILS
        nidx = len(self.verts)
        self.vert.append((.5*Overts[Opoly[i][0]] + .5*Overts[Opoly[i][1]]).tolist())
        self.vert.append((.5*Overts[Opoly[i][1]] + .5*Overts[Opoly[i][2]]).tolist())
        self.vert.append((.5*Overts[Opoly[i][2]] + .5*Overts[Opoly[i][0]]).tolist())
        self.poly[4*i] = [Opoly[i][0], nidx, nIdx+2]
        self.poly[4*i + 1] = [nidx, Opoly[i][1],nIdx+1]
        self.poly[4*i + 2] = [nidx+1, Opoly[i][2], nidx+2]
        self.poly[4*i + 3] = [nidx, nidx+1, nidx+2]

class FEMcalc(object):
  """The base class which calculates FEM for a problem
  """
  """
  #class attributes
  dim = 0 #dimension of the problem
  polyDeg = 0 #degree of polynomial approximation
  reg = 0 #regularity of basis functions
  domain = None #mesh of domain
  bdry = {0:None} #dictionary defining what to do with B.V.'s
  blin = [[],[]] #list of operations the bilinear form applies
  lin = [[],[]] #linear form
  funcmod = None #module with user defined functions
  refPolys = [] #polynomials on the reference element
  elemStiffMat = np.array([]) #element stiffness matrix
  dofs = [] #indices of corresponding basis function
  globStiffMat = None#global stiffness matrix. TODO: sparse
  rhs = None #Right hand side of linear equation
  soln = None #The coefficients of the solution.
  """

  def __init__(self,fname="FEM.setup"):
    self.dim = 0 #dimension of the problem
    self.polyDeg = 0 #degree of polynomial approximation
    self.reg = 0 #regularity of basis functions
    self.domain = None #mesh of domain
    self.bdry = {0:None} #dictionary defining what to do with B.V.'s
    self.blin = [[],[]] #list of operations the bilinear form applies
    self.lin = [[],[]] #linear form
    self.funcmod = None #module with user defined functions
    self.refPolys = [] #polynomials on the reference element
    self.elemStiffMat = np.array([]) #element stiffness matrix
    self.dofs = [] #indices of corresponding basis function
    self.globStiffMat = None#global stiffness matrix. 
    self.rhs = None #Right hand side of linear equation
    self.soln = None #The coefficients of the solution.
    self.parseInputFile(fname)
    self.genRefBasis()
    self.genElemStiffMat()
    self.genGlobalStiffMat()
    self.genRHS()
    self.findSolution()
   
  def refMesh(self):
    self.domain.refineMesh()
    self.genGlobalStiffMat()
    self.genRHS()
    self.findSolution()

  def parseInputFile(self, fname):
    """See example .setup file for how this parses
    stuff -> note that lines beginning with # are comments"""
    with open(fname,'r') as fl:
      secType = -1 #section types are 
      #0-MESH:, 1-SPACE:, 2-FEM_FORM: 3-FUNCTION_FILE
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
        if(rdln[0] == "SPACE:"):
          secType = 1
          continue
        if(rdln[0] == "FEM_FORM:"):
          secType = 2
          continue

        #parse the sections
        if(secType == 0):#MESH:
          if(rdln[0] == "file="):
            self.domain = mesh(rdln[1])
            self.dim = self.domain.verts.shape[1]

          elif(rdln[0].isdigit()):
            self.bdry[int(rdln[0])] = rdln[1]

        elif(secType ==1):#SPACE:
          if(rdln[0] == "PD"):
            try:
              self.polyDeg = int(rdln[1])
            except ValueError:
              raise NameError("polynomial degree must be an integer")
          elif(rdln[0] == "REG"):
            try:
              self.reg = int(rdln[1])
            except ValueError:
              raise NameError("regularity must be an integer")

        elif(secType == 2):#FEM_FORM:
          if(rdln[0] == "a(u,v)="):
            rdln.pop(0)
            pos = 0
            inu = -1
            while (rdln[0] != "|"):
              nelem = rdln.pop(0)
              if(nelem == "("):
                self.blin[0].append([])
                self.blin[1].append([])
                inu = 0
              elif(nelem == "grad"):
                self.blin[inu][pos].append('grad')
                rdln.pop(0)
              elif(nelem == ","):
                inu  = 1 
              elif(nelem == "u" or nelem == "v"):
                self.blin[inu][pos].append('eval')
                self.blin[inu][pos].append(nelem)
              elif(nelem == "+" and inu >= 0):
                self.blin[inu][pos].append('plus')
              elif(nelem == "-"):
                self.blin[inu][pos].append('sub')
              elif(nelem == ")"):
                inu = -1

          elif(rdln[0] == "L(v)="):
            rdln.pop(0)
            pos = 0
            inu = -1
            while (rdln[0] != "|"):
              nelem = rdln.pop(0)
              if(nelem == "("):
                self.lin[0].append([])
                self.lin[1].append([])
                inu = 0
              elif(nelem == "grad"):
                self.lin[inu][pos].append('grad')
                rdln.pop(0)
              elif(nelem == ","):
                inu  = 1 
              elif(nelem == "v"):
                self.lin[inu][pos].append('eval')
                self.lin[inu][pos].append('v')
              elif(nelem == "+" and inu >= 0):
                self.lin[inu][pos].append('plus')
              elif(nelem == "-"):
                self.lin[inu][pos].append('sub')
              elif(nelem == ")"):
                inu = -1
              elif(inu == 0):
                self.lin[inu][pos].append('eval')
                self.lin[inu][pos].append(nelem)

          elif(rdln[0] == "funcfile="):
            self.funcmod = importlib.import_module(rdln[1])
            self.funcmod = reload(self.funcmod)

    #check that regularity is possible with polynomial degree
    if(self.polyDeg < 1):
      self.polyDeg = 1

    if(self.dim == 1):
      if(self.polyDeg + 1 < 2*self.reg + 2):
        self.polyDeg = 2*self.reg + 1
        print ("Polynomial degree was insufficient to achieve\n \
          requested regularity. Polynomial degree reset to %d\n"%self.polyDeg)
    elif(self.dim == 2):
      if((self.polyDeg + 1)*(self.polyDeg+2) < 3*(self.reg + 1)*(self.reg + 2)):
        self.polyDeg = math.ceil(math.sqrt(3*(self.reg+1)*(self.reg+2)-1.5) - 1.5)
        print ("Polynomial degree was insufficient to achieve\n \
          requested regularity. Polynomial degree reset to %d\n"%self.polyDeg)

  def genRefBasis(self):
    """Generates basis functions on the reference element

    Note that in 1D, the reference element is [0,1],
    and in 2D the reference element is the triangle
    defined by the points (0,0), (0,1), (1,1)
    
    The degrees of freedom are below
    default 1D DOF:
    enough derivatives to satisfy regularity on boundary,
    and regulary spaced points between
    The returned list will correspond to the following elements being 1:
    [f[0] f'[0] ... f^(r)[0] f[1/n] f[2/n] ... f[n-1/n] f^(r)[1] ... f'[1] f[1]]
    default 2D DOF (supports only regularity = 0):
    uses equi-distant points along the border, as well as the corners.
    if needed, uses the midpoint as well
    The returned list will correspond to the following elements being 1:
    [f[0,0] f[0,1/n] ... f[0,1] f[1/n,1] ... f[1,1] f[1-1/n,1-1/n]... f[1/n,1/n] ]
    """
    
    if(self.dim == 0):
      raise NameError("Must define dimension before getting RefBasis")

    ptmp = polynom(self.dim, self.polyDeg)
    ptmp.setFlatCoef(())#generate dense polynomial of appropriate degree
    nDOF = ptmp.getFlatCoef(self.dim*(0,)).shape#this is a tuple
    coefMat = np.zeros(nDOF + nDOF,dtype=float)#nDOF by nDOF zeros
    nDOF = nDOF[0]
    if(self.dim == 1):
      rownum = 0
      pd = ptmp.copy()
      for i in range(self.reg+1):
        coefMat[rownum] = pd.getFlatCoef((0,))
        rownum += 1
        coefMat[-rownum] = pd.getFlatCoef((1,))
        pd = pd.diff()
      
      for rn in range(rownum ,nDOF - rownum):
        pt = (float(rn +1 - rownum)/float(nDOF - rownum),)
        coefMat[rn] = ptmp.getFlatCoef(pt)
      
    elif(self.dim == 2):
      if(self.reg > 0):
        raise NameError("Regularity degree %d not supported for 2d"%self.reg)
      
      #Generate all the DOF points
      rownum = -1
      nperSide = nDOF/3
      for i in range(nperSide):
        rownum += 1
        coefMat[rownum] = ptmp.getFlatCoef((0,0 + float(i)/nperSide))
      for i in range(nperSide):
        rownum += 1
        coefMat[rownum] = ptmp.getFlatCoef((0 + float(i)/nperSide,1))
      for i in range(nperSide):
        rownum += 1
        coefMat[rownum] = ptmp.getFlatCoef((1 - float(i)/nperSide,1 - float(i)/nperSide))
      
      
      if(nDOF%3 == 1):
        rownum += 1
        coefMat[rownum] = ptmp.getFlatCoef((x,y))
      
    else:
      raise NameError("%d-d polynomials not supported"%self.dim)
     
    #will solve and put coefficients along the rows
    coefMat = la.solve(coefMat,np.eye(nDOF)).transpose()
    self.refPolys = [
      polynom(self.dim,self.polyDeg,coefMat[i]) 
      for i in range(nDOF)]
    self.elemStiffMat = np.zeros(2*(nDOF,))

  def genElemStiffMat(self):
    """Generates the element stiffness matrix
    """
    
    if(len(self.elemStiffMat) == 0):
      raise NameError("Must generate polys on ref element before calling genElemStiffMat")

    rownum = -1
    #for each polynomial, on right hand side of inner product
    for polyU in self.refPolys:
      rownum += 1
      colnum = -1
      #for each term in self.blin form
      polyLeft = []
      for i in range(len(self.blin[0])):
        #-----Left hand side--------
        #if its a list of stuff
        if(type(self.blin[0][i]) == list):
          #for each element in this list
          for j in range(len(self.blin[0][i])):
            if(self.blin[0][i][j] == "grad"):
              #calculate the gradient
              for dim in range(self.dim):
                polyLeft.append(polyU.diff(dim))
            else:
              raise NameError("%s not currently supported"%self.blin[0][i][j])
        #else if its a string
        else:
          raise NameError("%s not currently supported"%self.blin[0][i])
        
        #for each polynomial on left hand side of inner product
        for polyV in self.refPolys: 
          colnum += 1    
          #----Right hand side-------
          #if its a list of stuff
          polyRight = []
          if(type(self.blin[1][i]) == list):
            #for each element in this list
            for j in range(len(self.blin[1][i])):
              if(self.blin[1][i][j] == "grad"):
                #calculate the gradient            
                for dim in range(self.dim):
                  polyRight.append(polyV.diff(dim))
              else:
                raise NameError("%s not currently supported"%self.blin[0][i][j])
          else:
            raise NameError("%s not currently supported"%self.blin[0][i][j])
            
          #compute entry in stiffness matrix
          if(len(polyLeft) != len(polyRight)):
            print polyLeft
            print polyRight
            raise NameError("Inner product has different dimensions on left and right")
          
          p2int = polynom(self.dim,0)
          for i2 in range(len(polyLeft)):
            p2int = p2int + polyLeft[i2]*polyRight[i2]
          
          self.elemStiffMat[rownum,colnum] = p2int.integrate((),())

  def genGlobalStiffMat(self):
    """Generates the global stiffness matrix
    """
    globDOF = 0
    #find the global degrees of freedom, and arrange them appropriately
    for eleNum in range(len(self.domain.poly)):
      eleDOF = []
      eleVerts = self.domain.poly[eleNum]
      for i in range(len(self.refPolys)):
        #figure out if this DOF corresponds with a vertex
        vertNum = -1
        if(i == 0):
          vertNum = 0
        elif(self.dim == 1 and i == len(self.refPolys)-1):
          vertNum = 1
        elif(self.dim == 2):
          if(i==(len(self.refPolys)/3)):
            vertNum = 1
          if(i==2*(len(self.refPolys)/3)):
            vertNum = 2
        #if a vertex
        if(vertNum >= 0):
          dofEle = min(self.domain.v_res[eleVerts[vertNum]])
          if(self.bdry[ self.domain.v_is_bndry[ eleVerts[vertNum] ] ] == "HD"):
            eleDOF.append(-1)
          elif(self.domain.v_is_bndry[ int(eleVerts[vertNum]) ] != 0):
            raise NameError("boundry condition not supported")
          elif(dofEle == eleNum):
            eleDOF.append(globDOF)
            globDOF += 1
          else:
            vidx = self.domain.poly[dofEle].tolist().index(eleVerts[vertNum])
            if(self.dim == 1):
              eleDOF.append(self.dofs[dofEle][-vidx])
            else:
              eleDOF.append(self.dofs[dofEle][vidx * (len(self.refPolys)/3)])
        #if edge
        elif(self.dim == 2 and i < len(self.refPolys) - 1):
          tmp = i/(len(self.refPolys)/3)
          vn1 = eleVerts[tmp]
          vn2 = eleVerts[(tmp+1)%3]
          #list of all triangles the edge is in
          triIn = [k for k in self.domain.v_res[vn1] if k in self.domain.v_res[vn2]]
          if(len(triIn) == 1):#the edge is on the boundary
            eleDOF.append(-1)
          elif(min(triIn) == eleNum):
            eleDOF.append(globDOF)
            globDOF += 1
          else:
            dofEle = min(triIn)
            vidx = self.domain.poly[dofEle].tolist().index(vn1)
            eleDOF.append(self.dofs[dofEle][vidx * (len(self.refPolys)/3) + i%(len(self.refPolys)/3)])
          
        #if interior
        else:
          eleDOF.append(globDOF)
          globDOF += 1
      
      self.dofs.append(eleDOF)

    self.globStiffMat = sparse.lil_matrix((globDOF,globDOF))

    if(self.reg != 0):
      raise NameError("%d-order regularity not supported",self.reg)

    #for each mesh interval/triangle
    for eleIdx in range(len(self.domain.poly)):
      scale = 1. 
      vs = self.domain.poly[eleIdx]
      if(self.dim == 1):
        x2 = self.domain.verts[vs[1]][0]
        x1 = self.domain.verts[vs[0]][0]
        #1/h
        scale = abs((x2-x1))
        scale = scale/scale**2 #2 gradients and an integral
        
      elif(self.dim == 2):
        x3 = self.domain.verts[vs[2]][0]
        x2 = self.domain.verts[vs[1]][0]
        x1 = self.domain.verts[vs[0]][0]
        y3 = self.domain.verts[vs[2]][1]
        y2 = self.domain.verts[vs[1]][1]
        y1 = self.domain.verts[vs[0]][1]
        #1/area of triangle = abs(jacobian)
        scale = (x3 - x2 )*( y2 - y1) + (x2 - x1)*(y3 - y2)
        scale = abs(scale)/scale**2
        
      for i1 in range(len(self.dofs[eleIdx])):
        rowIdx = self.dofs[eleIdx][i1]
        if(rowIdx < 0):
          continue
        for i2 in range(i1, len(self.dofs[eleIdx])):
          colIdx = self.dofs[eleIdx][i2]
          if(colIdx < 0):
            continue
          
          val = scale * self.elemStiffMat[i1,i2]
          self.globStiffMat[rowIdx,colIdx] += val
          if(rowIdx != colIdx):
            self.globStiffMat[colIdx,rowIdx] += val

    self.globStiffMat = self.globStiffMat.tocsr()

  def genRHS(self):
    """Generates the RHS for solve
    """
    self.rhs = np.zeros(self.globStiffMat.shape[0])
    #calculate L(v)
    for i in range(len(self.lin[0])):
      #-----Left hand side--------
      #if its a list of stuff
      lhsFunc = None
      if(type(self.lin[0][i]) == list):
        #for each element in this list
        skip = False
        for j in range(len(self.lin[0][i])):
          if(skip):
            skip = False
            continue
          if(self.lin[0][i][j] == "eval"):
            #save the function
            skip = True
            lhsFunc = getattr(self.funcmod,self.lin[0][i][j+1])
          else:
            raise NameError("%s not currently supported"%self.lin[0][i][j])
      #else if its a string
      else:
        raise NameError("%s not currently supported"%self.lin[0][i])
      
      #----Right hand side-------
      for eleIdx in range(len(self.domain.poly)):
        #scaling factor
        scale = 1. 
        #transformation into reference element
        trans = None
        vs = self.domain.poly[eleIdx]
        if(self.dim == 1):
          x2 = self.domain.verts[vs[1]][0]
          x1 = self.domain.verts[vs[0]][0]
          #1/h
          scale = abs((x2-x1))
          #trans = lambda x: (float(x[0]) - x1)/(x2 - x1)
          trans = lambda x: (x[0]*(x2-x1) + x1,) 
          
        elif(self.dim == 2):
          x3 = self.domain.verts[vs[2]][0]
          x2 = self.domain.verts[vs[1]][0]
          x1 = self.domain.verts[vs[0]][0]
          y3 = self.domain.verts[vs[2]][1]
          y2 = self.domain.verts[vs[1]][1]
          y1 = self.domain.verts[vs[0]][1]
          #1/area of triangle = abs(jacobian)
          scale = (x3 - x2 )*( y2 - y1) + (x2 - x1)*(y3 - y2)
          scale = abs(scale)
          #trans = lambda x: ( 
          #  (x1*x[1] - x2*x[1] - x[0]*y1 + x2*y1 + x[0]*y2 - x1*y2)/
          #    (x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3)
          #  , 
          #  ((x3-x2)*(x[1] - y1) + (x[0]-x1)*(y2-y3))/
          #    (x3*(y2-y1) + x2*(y1-y3) + x1*(y3-y2))
          #  )
          trans = lambda x: ( 
            (x3 - x2)*x[0] + (x2 - x1)*x[1] + x1
            , 
            (y3-y2)*x[0] + (y2 - y1)*x[1] + y1
            )

        for i1 in range(len(self.dofs[eleIdx])):
          rowIdx = self.dofs[eleIdx][i1]
          if(rowIdx < 0):
            continue
          
          #if its a list of stuff
          rhsFunc = None
          if(type(self.lin[1][i]) == list):
            #for each element in this list
            skipv = False
            for j in range(len(self.lin[1][i])):
              if(skipv):
                skipv = False
                continue
              if(self.lin[1][i][j] == "eval"):
                #save the function
                skipv = True
                rhsFunc = self.refPolys[i1]
              else:
                raise NameError("%s not currently supported"%self.lin[1][i][j])
          else:
            raise NameError("%s not currently supported"%self.lin[1][i][j])
        
          if(self.dim == 1):
            tmp = scale *integrate.quad(
              lambda x: lhsFunc(trans((x,))) * rhsFunc((x,)), 0,1
            )[0]
            self.rhs[rowIdx] += tmp
          elif(self.dim ==2):
            self.rhs[rowIdx] += scale *integrate.dblquad(
              lambda x,y: lhsFunc(trans((x,y))) * rhsFunc((x,y)), 0,1,lambda x: 0, lambda x: x
            )[0]

  def findSolution(self):
    """Finds the solution once everything has been created
    """
    self.soln = spla.spsolve(self.globStiffMat,self.rhs)
 
  def __call__(self,x):
    """Returns the value of the solution at x

    If x is not inside the mesh, will return nan
    """
    if(x in self.domain.verts):
      vidx = self.domain.verts.tolist().index(list(x))
      pidx = min(self.domain.v_res[vidx])
      vpidx = self.domain.poly[pidx].tolist().index(vidx)
      if(self.dim==1):
        solnIdx = self.dofs[pidx][-vpidx]
      elif(self.dim==2):
        solnIdx = self.dofs[pidx][vpidx * (len(self.refPolys)/3)]
      #if its on the boundary
      if(solnIdx < 0):
        return 0.
      else:
        return self.soln[solnIdx]
    else:
      #determine which element it lies in
      eleNum = 0
      while(not self.__inside__(eleNum,x)):
        eleNum += 1
      
      if(eleNum > len(self.domain.poly)):
        return float('nan')
      
      val = 0.
      trans = None
      vs = self.domain.poly[eleNum]
      if(self.dim == 1):
        x2 = self.verts[vs[1]][0]
        x1 = self.verts[vs[0]][0]
        #1/h
        trans = lambda x: (float(x[0]) - x1)/(x2 - x1)
      elif(self.dim == 2):
        x3 = self.verts[vs[2]][0]
        x2 = self.verts[vs[1]][0]
        x1 = self.verts[vs[0]][0]
        y3 = self.verts[vs[2]][1]
        y2 = self.verts[vs[1]][1]
        y1 = self.verts[vs[0]][1]
        #1/area of triangle = abs(jacobian)
        trans = lambda x: ( 
          (x1*x[1] - x2*x[1] - x[0]*y1 + x2*y1 + x[0]*y2 - x1*y2)/
            (x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3)
          , 
          ((x3-x2)*(x[1] - y1) + (x[0]-x1)*(y2-y3))/
            (x3*(y2-y1) + x2*(y1-y3) + x1*(y3-y2))
            )
          
      for i in range(len(self.refPolys)):
        solnIdx = self.dofs[eleNum][i]
        if(solnIdx < 0):
          continue
        print self.soln[solnIdx]
        print self.refPolys(trans(x))
        val += self.soln[solnIdx] * self.refPolys(trans(x))
      print val
      return val

  def __inside__(self,eleNum,x):
    vs = self.domain.poly[eleNum]
    pts = np.zeros((self.dim+1,self.dim))
    x = np.array(x)
    for i in range(self.dim+1):
      pts[i] = self.domain.verts[vs[i]]

    if(self.dim == 1):
      #this implies that x was greater than one point and less than the other point
      return (pts[0][0] - x)*(pts[1][0] - x) < 0
    elif(self.dim == 2):
      v0 = pts[2] - pts[0]
      v1 = pts[1] - pts[0]
      v2 = x - pts[0]
      #use barycentric coordinate method
      dot00 = np.dot(v0, v0)
      dot01 = np.dot(v0, v1)
      dot02 = np.dot(pts, v2)
      dot11 = np.dot(v1, v1)
      dot12 = np.dot(v1, v2)
      
      invDenom = 1. / (dot00 * dot11 - dot01 * dot01)
      u = (dot11 * dot02 - dot01 * dot12) * invDenom
      v = (dot00 * dot12 - dot01 * dot02) * invDenom
      
      return ((u >= 0) and (v >= 0) and (u + v < 1))

  def pltSoln(self,xargs = []):
    """plots the solution.
    
    Parameters:
    xargs -> 
    a list containing a function handle and string.
    This will be plotted along with the solution
    """
    if(self.dim == 1):
      x = self.domain.verts
      y = np.zeros(len(x))
      if(len(xargs) == 2):
        y1 = np.zeros(len(x))
      for i in range(len(x)):
        y[i] = self(x[i])
        if(len(xargs) == 2):
          y1[i] = xargs[0](x[i])
      
      plt.figure()
      plt.plot(x,y,label = 'num. soln.')
      if(len(xargs) == 2):
        plt.plot(x,y1,label = xargs[1])
        plt.legend(loc=2)

    elif(self.dim == 2):

      fig = plt.figure()
      ax = fig.add_subplot(111, projection='3d')
      
      x = self.domain.verts[:,0]
      y = self.domain.verts[:,1]
      z = np.zeros(len(x))
      for i in range(len(z)):
        z[i] = self((x[i],y[i]))
      ax.plot_trisurf(x, y, z, linewidth=0.2)

    plt.show()

