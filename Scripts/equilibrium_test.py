import numpy as np


# a=[[4,2,3],[3,-5,2],[-2,3,8]]
# a=np.asarray(a)
# b=np.asarray([8,-14,27])
# #need strictly upper triangle (no diag)
# u=np.triu(a,k=1)
# #but sum of lower triangle and diag
# dl=np.tril(a)
# #which we only use by inverse
# dl_inv=np.linalg.inv(dl)

# p=-dl_inv.dot(u)
# q=dl_inv.dot(b)
# x=np.zeros(3)
# for i in range(10):
#     x=p.dot(x)+q
# print(x)


#Robertson Problem
#https://arxiv.org/pdf/2103.15341.pdf
b=np.zeros(3)
r=[0.04,3e7,1e4]
x=np.asarray([0.9,0.06,0.12])
for i in range(10):
    a=np.asarray(
        [[-r[0],0,r[2]*x[1]],
        [r[0],-r[1]*x[1],-r[2]*x[2]],
        [0,r[1]*x[1]*x[1],0]])
    print(a)
    #need strictly upper triangle (no diag)
    u=np.triu(a,k=1)
    #but sum of lower triangle and diag
    dl=np.tril(a)
    print(dl)
    #which we only use by inverse
    dl_inv=np.linalg.inv(dl)

    p=-dl_inv.dot(u)
    q=dl_inv.dot(b)

    x=p.dot(x)+q
    print(x)
