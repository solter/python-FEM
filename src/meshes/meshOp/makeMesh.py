import numpy as np

def make1dMesh(n, xl = 0, xr = 1, name = "1deq"):

  with open("%s%d.node"%(name,n),'w') as f:
    s = "#Uniform 1d mesh on [%g,%g] with %d partitions\n"%(xl,xr,n)
    f.write(s)
    s = "%d 1 0 1\n"%(n+1)
    f.write(s)
    xs = np.linspace(xl, xr, n+1)
    for i in range(len(xs)):
      s = "%d %f %d\n"%(i+1, xs[i], xs[i] == xl or xs[i] == xr)
      f.write(s)

  with open("%s%d.ele"%(name,n),'w') as f:
    s = "#Uniform 1d mesh on [%g,%g] with %d partitions\n"%(xl,xr,n)
    f.write(s)
    s = "%d 2 0\n"%(n)
    f.write(s)
    for i in range(1,n+1):
      s = "%d %d %d\n"%(i,i,i+1)
      f.write(s)

def make2dMesh(n):
  with open("2deq%d.node"%n,'w') as f:
    s = "#Uniform 1d mesh on [0,1]x[0,1] with %d internal nodes per x\n"%n
    f.write(s)
    s = "%d 2 0 1\n"%((n+2)**2)
    f.write(s)
    xs = np.linspace(0,1,n+2)
    for i in range(len(xs)):
      for j in range(len(xs)):
        s = "%d %f %f %d\n"%(i*len(xs) + j + 1, xs[j], xs[i], xs[i] == 0 or xs[i] == 1 or xs[j] == 0 or xs[j] == 1)
        f.write(s)

  with open("2deq%d.ele"%n,'w') as f:
    s = "#Uniform 1d mesh on [0,1]x[0,1] with %d internal nodes per x\n"%n
    f.write(s)
    s = "%d 3 0\n"%(2*(n+1)**2)
    f.write(s)
    polnum = 0
    for i in range(1,n+2):
      for j in range(1,n+2):
        polnum += 1
        s = "%d %d %d %d\n"%(polnum,(i-1)*(n+2) + j,(i-1)*(n+2) + j+1,i*(n+2)+j)
        f.write(s)
        polnum += 1
        s = "%d %d %d %d\n"%(polnum,(i-1)*(n+2) + j+1,i*(n+2)+j+1,i*(n+2)+j)
        f.write(s)

