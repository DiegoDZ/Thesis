#Take advantage of the periodic boundary conditions in the creation of C(t) matrix.
# !!! Only for system 'fluid' with periodic boundary conditions.

import numpy as np

Ct     = np.loadtxt('Ct-PBC-2500steps-noStat.dat')
nSteps = np.shape(Ct)[0]

nNodes      = 60              
nBlocks     = 1
nVar        = int(np.sqrt(nBlocks))
dim         = nVar * nNodes

def pbc(C):
    Cstat = np.zeros((nNodes,nNodes))
    for i in range(nNodes):
        for j in range(nNodes):
            for k in range(nNodes):
                Cstat[i,j] += C[(i+k)%nNodes,(j+k)%nNodes]
    return Cstat/nNodes

#Change format: vector-> matrix
def reshape_vm(A):
    B = A.reshape(nBlocks,nNodes*nNodes).reshape(nVar,nVar,nNodes,nNodes).transpose(0,2,1,3).reshape(dim,dim)
    return B

#Change format: matrix-> vector
def reshape_mv(A):
    B = A.reshape(nVar,nNodes,nVar,nNodes).swapaxes(1,2).ravel()
    return B

CtStat   = np.zeros((nSteps, nNodes ** 2))
for i in range(nSteps):
    print i 
    C = reshape_vm(Ct[i,:])
    CtStat[i,:] = reshape_mv(((pbc(Ct)) + pbc(C.T)) / 2)

np.savetxt('Ct-PBC-2500steps.dat', CtStat)
