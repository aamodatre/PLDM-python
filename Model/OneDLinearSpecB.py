#/usr/bin/python3
""" Modelled for Linear Absorption Spectrum - Subotnik JCTC"""
"""Ref.: DOI: 10.1021/acs.jctc.5b00510"""
import numpy as np


class parameters():
   NSteps = 600 #int(2*10**6)
   NTraj = 10
   dtN = 2
   dtE = dtN/20
   NStates = 2
   M = 1728.26
   initStatef = 0
   initStateb = 1
   nskip = 1
   

def bathparams(N):
    ω = np.zeros(N+1, dtype = np.float64)
    ω[0] = 0.002278 # 500 cm-1
    ωc = 2*ω[0] 
    m = 1728.26
    eta = (ω[0]*m)**(-1)
    c = np.zeros(N+1)
    
    for i in range(1, N+1):
        ω[i] = -ωc*np.log((i-(3/2))/(N-1))
        
        c[i] = ω[i]*np.sqrt((2*eta*m*ωc)/(np.pi*(N-1)))    


def Hel(R):

    VMat = np.zeros((2,2))

    VMat[0,0] = 
    VMat[1,0] = C * np.exp(-D * R**2)
    VMat[0,1] = VMat[1,0]
    VMat[1,1] = -A * np.exp(-B * R**2) + E0
