#------------------------------------------------------------------------------
#                            buildCt.py
#------------------------------------------------------------------------------
# This script builds the matrices of correlations only with a bunch of
# correlations computing with LAMMPS. 
# The matrices of correlations built are:
# <gx(t)gx>, <Sxz(t)Sxz>, <Sxz(t)Fx>, <Fx(t)Sxz> and <Fx(t)Fx>
#------------------------------------------------------------------------------
#                         Author   : @DiegoDZ
#                         Date     : January  2018
#                         Modified : September 2018
#                         Run      : python buildCt.py  (Run with python 2.7)
#------------------------------------------------------------------------------

import numpy as np

#----------------------------------
#Inputs
#----------------------------------
nSteps   = 500                              #nCols = nNodes ** 2
dt       = 0.004
nNodes   = 34                               #Number of fluid nodes 
nBlocks  = 1                                #Do not confuse with the number of blocks we will use to build the matrix
                                                         #of correlations C(t). This number of blocks refers to the number of different
                                                         #correlation files (i.e. <gxgx(t)>, <gxrho(t)> ...).
                                                         #In this case the number of blocks is only one because the matrix of
                                                         #correlations will be <gxgx(t)>.
sBlocks       = int(np.sqrt(nBlocks))
dim           = sBlocks * nNodes
#We define the number of nodes for the blocks calculated with LAMMPS. Note the different between the number of nodes of C(t) and the number of nodes of block1 and block2.

#----------------------------------
#Subrutines
#----------------------------------

#Change format: vector-> matrix  (For matrices of correlations)
def reshape_vm(A):
    B = A.reshape(nBlocks,nNodes*nNodes).reshape(sBlocks,sBlocks,nNodes,nNodes).transpose(0,2,1,3).reshape(dim,dim)
    return B

#Change format: matrix-> vector  (For matrices of correlations)
def reshape_mv(A):
    B = A.reshape(sBlocks,nNodes,sBlocks,nNodes).swapaxes(1,2).ravel()
    return B

def step(t):
    step = int(round(t/dt))
    return step


#-------------------------------------------MATRICES--------------------------------
    
Ct      = np.loadtxt('Ct-x2Box-WALLS-untilt1-AVG.dat')
SxzFxt  = np.loadtxt('corr-SxzFx-x2Box-WALLS-untilt1-AVG.dat')
#SxzFxt  = np.loadtxt('corr-SxzFx-x2Box-WALLS-untilt1.dat.1')



CtStat       = np.zeros((nSteps, nNodes ** 2))
corr_SxzSxzt = np.zeros((nSteps, nNodes ** 2))
corr_SxzFxt  = np.zeros((nSteps, nNodes ** 2))
corr_FxSxzt  = np.zeros((nSteps, nNodes ** 2))
corr_FxFxt   = np.zeros((nSteps, nNodes ** 2))

for k in range(nSteps):
    print 'step='+str(k)
    C           = np.zeros((nNodes, nNodes))

    CSim = 0.5 * (Ct[k,:].reshape(nNodes-1, nNodes-1)+(Ct[k,:].reshape(nNodes-1, nNodes-1)).T)
    CStat = 0.5 * (CSim + np.rot90(CSim,2))
    C[1:,1:] = CStat
    CtStat[k,:] = reshape_mv(C)
    #C[1:,1:] = Ct[k,:].reshape(nNodes-1,nNodes-1)
    #CtStat[k,:] = 0.5 * reshape_mv(C + C.T) 
    #CtStat[k,:] = reshape_mv(C)
    
    A                 = SxzFxt[k,:].reshape(2*(nNodes-1), 2*(nNodes-1))
    B                 = np.zeros((2*nNodes, 2*nNodes))
    B[1:nNodes,1:nNodes]    = A[:nNodes-1, :nNodes-1]
    B[1:nNodes,nNodes+1:]   = A[:nNodes-1, nNodes-1:]
    B[nNodes+1:, 1:nNodes]  = A[nNodes-1:, :nNodes-1]
    B[nNodes+1:, nNodes+1:] = A[nNodes-1:, nNodes-1:]
   
    #stat.
    B[:nNodes, :nNodes] = 0.5 * (B[:nNodes, :nNodes] + np.rot90(B[:nNodes, :nNodes], 2)) 

    corr_SxzSxzt[k,:] = reshape_mv(B[:nNodes, :nNodes])
    corr_SxzFxt[k,:]  = reshape_mv(B[:nNodes, nNodes:])
    corr_FxSxzt[k,:]  = reshape_mv(B[nNodes:, :nNodes])
    corr_FxFxt[k,:]   = reshape_mv(B[nNodes:, nNodes:])

np.savetxt('Ct-x2Box-WALLS.dat',     CtStat)
np.savetxt('SxzSxz-x2Box-WALLS.dat', corr_SxzSxzt)
np.savetxt('SxzFx-x2Box-WALLS.dat',  corr_SxzFxt)
np.savetxt('FxSxz-x2Box-WALLS.dat',  corr_FxSxzt)
np.savetxt('FxFx-x2Box-WALLS.dat',   corr_FxFxt)


#np.savetxt('SxzSxz-t0-new', reshape_vm(corr_SxzSxzt[0,:]))
#np.savetxt('SxzFx-t0-new', reshape_vm(corr_SxzFxt[0,:]))
#np.savetxt('FxSxz-t0', reshape_vm(corr_FxSxzt[0,:]))
#np.savetxt('FxFx-t0', reshape_vm(corr_FxFxt[0,:]))

#EOF

