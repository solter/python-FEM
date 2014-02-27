import polynom as p
reload(p)

p1 = p.polynomial(2);
p2 = p.polynomial(2);
p1.addElem((0,0),1)
p1.addElem((1,1),4)
p2.addElem((0,0),-2)
p2.addElem((3,1),2)

