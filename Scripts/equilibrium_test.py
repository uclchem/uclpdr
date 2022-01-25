import numpy as np


a=[[4,2,3],[3,-5,2],[-2,3,8]]
a=np.asarray(a)
b=np.asarray([8,-14,27])

#need strictly upper triangle (no diag)
u=np.triu(a,k=1)
#but sum of lower triangle and diag
dl=np.tril(a)
#which we only use by inverse
dl_inv=np.linalg.inv(dl)

p=-dl_inv.dot(u)
q=dl_inv.dot(b)

x=np.asarray([2.0,4.0,2.375])
for i in range(10):
    x=p.dot(x)+q
    print(x)