import numpy as np
import scipy as sp
#
#class polynomial(object):
#  """A class to deal with multidimensional polynomials.
#  """
#
#  def __init__(self, dim, deg = 1, coef = np.array([])):
#    """Initializes a polynomial
#    
#    PARAMETERS: 
#    dim : the dimension of the space (>2 may be buggy)
#    deg : the expected degree of the polynomial
#    coef : the coeffiecient matrix. If empty, will be initialized to zeros
#    """
#    self.dim = dim
#    self.deg = deg
#    if(coef.size > 0):
#      self.coef = coef
#    else:
#      self.coef = np.zeros((self.deg+1)*np.ones(self.dim))
#
#  def addElem(self, idx, coef):
#    """Add nonzero polynomial element
#    
#    PARAMETERS:
#    idx : index of parameter. Must be a tuple of same length
#      as dimensions of polynomial.
#    coef : value of the coeficient
#
#    example:
#    for a 2-D polynomial, add the element 10*x1*x2^2:
#    addElem((1, 2), 10)
#    """
#   
#    self.deg = max([sum(idx), self.deg])
#    self.dim = self.dim
#    
#    #check that the appropriate number of dimensions is given
#    if(len(idx) != self.dim):
#      raise NameError("Index and polynomial dimension mismatch")
#
#    #expand the coef matrix if needed
#    if(len(self.coef) <= self.deg):
#      pred = len(self.coef)
#      self.coef = np.append(self.coef,np.zeros([self.deg +1 - pred, pred]),0)
#      self.coef = np.append(self.coef,np.zeros([self.deg+1, self.deg+1 - pred]),1)
#  
#    self.coef[idx] += coef
#
#  def __str__(self):
#    """Generates a string representation"""
#    toret =  self.__genStr((),"")
#    if(len(toret) > 2):
#      return toret[:-2]
#    else:
#      return "0"
#
#  def __genStr(self, idx, toret):
#    """Helper function for __str__"""
#    if self.dim == len(idx):
#      #if coef at idx is 0, don't print it
#      if(self.coef[idx] == 0):
#        return toret
#      #if coef at idx isn't 0, add it to return string
#      else:
#        toret += str(self.coef[idx])
#        for i in range(len(idx)):
#          if (idx[i] != 0):
#            toret += "*x%d^%d" %(i, idx[i])
#        return toret + " + "
#    else:#recurse
#      for i in range(self.coef.shape[len(idx)]):
#        tmp = idx + (i,)
#        toret = self.__genStr(tmp, toret);
#      return toret
#
#  def __add__(p1, p2):
#    dim = p1.dim
#    if(dim != p2.dim):
#      raise NameError("polynomial dimension mismatch")
#    if(p1.deg < p2.deg):
#      tmp = p1
#      p1 = p2
#      p2 = tmp
#
#    toret = polynomial(dim,p1.deg,p1.coef.copy())
#    toret.coef[0:p2.deg + 1, 0:p2.deg + 1] += p2.coef
#    return toret
#

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
    
