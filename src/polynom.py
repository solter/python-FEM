import numpy as np
import scipy as sp

class polynomial(object):
  """A class to deal with multidimensional polynomials.
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
      self.coef = np.zeros((self.deg+1)*np.ones(self.dim))

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
      self.coef = np.append(self.coef,np.zeros([self.deg +1 - pred, pred]),0)
      self.coef = np.append(self.coef,np.zeros([self.deg+1, self.deg+1 - pred]),1)
  
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
    dim = p1.dim
    if(dim != p2.dim):
      raise NameError("polynomial dimension mismatch")
    if(p1.deg < p2.deg):
      tmp = p1
      p1 = p2
      p2 = tmp

    toret = polynomial(dim,p1.deg,p1.coef.copy())
    toret.coef[0:p2.deg + 1, 0:p2.deg + 1] += p2.coef
    return toret



