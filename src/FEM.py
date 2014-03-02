import numpy as np
import scipy as sp

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
    if(coef.size > 0):
      self.coef = coef
    else:
      self.coef = np.zeros((self.deg+1)*np.ones(self.dim),dtype=float)

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
    return self.__add__(self,-1*p2)

  def __rsub__(self,p2):
    return self.__sub__(p2,self)

  def __mul__(p1, p2):
    #transform scalar into polynomial
    if(np.isscalar(p2)):
      tmp = p2
      p2 = polynom(self.dim,0)
      p2.addElem(tuple(np.zeros(self.dim,dtype=int)),tmp)
   
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
    return __mul__(p1,p2)

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
        toret = polynom(self.dim, self.deg + 1, np.insert(\
            np.insert(self.coef,0-var, np.zeros(self.deg+1),0),\
            -1+var,np.zeros(self.deg + 2),1))
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
      if(x[0] == 't' and x[1] == 't'):
        raise NameError("One value must be a number")
      elif(x[0] == 't' or x[1] == 't'):
        toret = polynom(1,self.deg)

      for idx in zip(self.coef.nonzero()[0],self.coef.nonzero()[1]):
        if(x[0] == 't'):
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
      Examples:
        1D: 
        integrate((0,),(1,)) = int_0^1 dx (returns a number)
        integrate((0,),('t',)) = int_0^t dx (returns a polynomial)

        2D:
        integrate((0,0),(1,2)) = int_0^1 int_0^2 dx0 dx1 (returns a number)
        integrate((0,0),(1,'t')) = int_0^1 int_0^x1 dx0 dx1 (returns a number)
        integrate((0,0),('t','t')) = int_0^t int_0^x1 dx0 dx1 (returns a 1d polynomial)
    """
    
    if(self.dim == 1):
      p = self.diff(0,False)
      if(x1[0] == x2[0]):
        return 0
      elif(x1[0] == 't'):
        return p - self(x0[0])
      elif(x0[0] == 't'):
        return self(x1[0]) - p
      else:
        return self(x1[0]) - self(x0[0])
    elif(self.dim == 2):
      pass
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
    with open(fname + ".node") as f:
      skip = 0
      for line in f:#for each line in the file
        #if not a comment
        if(line.split()[0][0] != "#"):
          if(getFormat):
            tmp = [int(s) for s in line.split()[0:5]]
            #generate matrix of empty 0's
            self.verts = np.zeros((tmp[0], tmp[1]))
            #make a vector of False
            self.v_is_bndry = [False for i in range(tmp[0])]
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
              self.v_is_bndry[idx] = True
    
    getFormat = True
    #parse the ele file
    with open(fname + ".ele") as f:
      for line in f:#for every line in the file
        #if not comment
        if(line.split()[0][0] != "#"):
          if(getFormat):
             tmp = [int(s) for s in line.split()[0:2]]
             #generate matrix of empty 0's
             self.poly = np.zeros((tmp[0],tmp[1]))
             getFormat = False
          else:
            s = line.split()
            idx = int(s.pop(0)) - 1#the poly number
            #populate the poly values
            for i in range(self.poly.shape[1]):
              nd = int(s.pop(0)) - 1
              self.poly[(idx,i)] = nd#add node to this element
              self.v_res[nd].append(idx)#add this element nodes residence list
    
