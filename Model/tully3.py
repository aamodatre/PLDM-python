import numpy as np

class parameters():
    """ Defining runtime parameters for the Model """
    NSteps = 600 #int(2*10**6)
    NTraj = 10      # Number of trajectories per job
    dtN = 2         # Sampling frequency
    dtE = dtN/20    
    NStates = 2     # Two site tully model
    M = 2000        # Mass in atomic units
    initState = 0   # Starting in the lower state
    nskip = 5      

def Hel(R):
    """ Defining the Electronic Hamiltonian
     as a function of nuclear coordinates """

    """ Potential energy in the diabatic representation """
    # Constant Parameters based on Tully 1990
    
    A = 6e-4
    B = 0.10
    C = 0.90
    
    VMat = np.zeros((2,2))
    
    VMat[0,0] = A
    VMat[1,1] = -A
    if R < 0:
        VMat[1,0] = B * np.exp(C * R)
    else:
        VMat[1,0] = B * (2 - np.exp(-C * R))
    VMat[0,1] = VMat[1,0]

    return VMat

def dHel0(R):
    return 0

def dHel(R):
    """ Defining the gradient matrix"""

    # Constant Parameters based on Tully 1990
    
    A = 6e-4
    B = 0.10
    C = 0.90
    
    dVMat = np.zeros((2,2,1))
    dVMat[0,0,0] = 0
    dVMat[1,1,0] = 0

    if R < 0:
        dVMat[1,0,0] = B*C*np.exp(C*R)
    else:
        dVMat[1,0,0] = B*C*np.exp(-C*R)
    
    dVMat[0,1,0] = dVMat[1,0,0]
    
    return dVMat

def initR():
    """ Initialising trajectories """
    R0 = -9.0
    P0 = 30
    alpha = 1.0
    sigR = 1.0/np.sqrt(2.0*alpha)
    sigP = np.sqrt(alpha/2.0)
    R = np.random.normal()*sigR + R0
    P = np.random.normal()*sigP + P0
    return np.array([R]), np.array([P])
