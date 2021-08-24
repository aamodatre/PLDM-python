#/usr/bin/python3
""" Modelled for Linear Absorption Spectrum """
"""Ref.: https://doi.org/10.1063/1.4884945"""

import numpy as np


class parameters():
   ESteps = 20
   NSteps = 1000 #int(2*10**6)
   NTraj = 2000
   dtN = 1
   dtE = dtN/ESteps
   NStates = 2
   M = 1728.26
   initStatef = 0
   initStateb = 1
   nskip = 1
   ωg = 0.013669
   ωe = 0.99*ωg
   D = 0.2
   epsilon = 5

def Hel(R):
    ωg = parameters.ωg
    ωe = parameters.ωe
    D = parameters.D
    m = parameters.M
    epsilon = 5

    VMat = np.zeros((2,2))
    VMat[0,0] = 0.5*m*(ωg**2)*(R**2)
    VMat[1,1] = 0.5*m*(ωe**2)*((R-D)**2) + epsilon
    V0 = (VMat[0,0]+VMat[1,1])/2.
    VMat[0,0] = VMat[0,0] - V0
    VMat[1,1] = VMat[1,1] - V0
    return VMat

def dHel(R):
    ωg = parameters.ωg
    ωe = parameters.ωe
    D = parameters.D
    m = parameters.M

    dVMat = np.zeros((2,2,1))
    # dVMat[0,0,0] = m*(ωg**2)*(R)
    dVMat[1,1,0] = m*(ωe**2)*((R-D)) - m*(ωg**2)*(R)
    
    return dVMat

def dHel0(R):
    ωg = parameters.ωg
    m = parameters.M
    return m*(ωg**2)*(R)

def initR():
    R0 = 0.0
    P0 = 0.0
    ωg = parameters.ωg
    # β = 10**6
    
    # sigP = np.sqrt( ωg / ( 2 * np.tanh( 0.5*β*ωg ) ) )
    
    sigP = np.sqrt( ωg / 2 ) # In the limit T = 0 K
    sigR = sigP/ωg
 
    R = np.random.normal()*sigR + R0
    P = np.random.normal()*sigP + P0

    return np.array([R]), np.array([P])