#/usr/bin/python3

"""Model for 2 Level Coupled-Dimer"""
""" Model -> Ref.:https://doi.org/10.1063/1.5006824"""
""" Parameters -> Ref.:https://doi.org/10.1063/1.5031788 """

import numpy as np
from scipy import constants as sc

class parameters():
   ESteps = 20
   NSteps = 2481 # For Spectroscopy simulation
   # NSteps = 8269 # For Populations simulation
   NTraj = 350
   dtN = 5.
   dtE = dtN/ESteps
   nskip = 1
   NStates = 4
   # NStates = 3
   
   # cj and ωj is an array of dimensions (100,1)
   cj = np.loadtxt("/scratch/aatre3/NAMD/PLDM-python/Discretization/cj.txt")
   ωj = np.loadtxt("/scratch/aatre3/NAMD/PLDM-python/Discretization/ωj.txt")
   ω = np.zeros((2*len(ωj)))
   c = np.zeros((2*len(cj)))
   ω[0:100] = ωj
   ω[100:200] = ωj
   c[0:100] = cj
   c[100:200] = cj
   ndof = int(len(c))

   initStatef = 2
   initStateb = 0

   M = 1
   Eh2c = 219474.6313632   # Hartree -> cm inverse
   ωc = 18.0/Eh2c   # 18 cm-1 in Hartree

def Hel(R):
   
   N = parameters.NStates
   cj = parameters.c
   ωj = parameters.ω
   ωc = parameters.ωc
   m = parameters.M
   Eh2c = parameters.Eh2c

   # On-Site Energies:
   ε = np.array([50,-50], dtype = np.float64)/Eh2c
   ε_bar = 10000.0/Eh2c  # in cm-1

   VConst= np.zeros((N, N), dtype=np.float64)
   VConst[0,0] = 0                           # State 00
   VConst[1,1] = ε_bar + ε[0]                # State 01
   VConst[2,2] = ε_bar + ε[1]                # State 10
   VConst[3,3] = 2*ε_bar + ε[0] + ε[1]       # State 11
   # Diabatic Couplings
   VConst[1,2] = 100.0/Eh2c                    # J-10
   VConst[2,1] = 100.0/Eh2c                    # J-01
   
   VCouple = np.zeros((N,N), dtype = np.float64)
 
   VCouple[1,1] = np.sum(cj[0:100]*R[0:100])
   VCouple[2,2] = np.sum(cj[100:200]*R[100:200])
   VCouple[3,3] = np.sum(cj*R)

   VMat = VConst + VCouple #+ VBath

   # V0 = (VMat[0,0]+VMat[1,1])/2.
   # VMat[0,0] = VMat[0,0] - V0
   # VMat[1,1] = VMat[1,1] - V0
   return VMat

def  dHel(R):
   N = parameters.NStates
   cj = parameters.c
   ωj = parameters.ω
   ndof =  parameters.ndof
   dVij = np.zeros((N,N,ndof), dtype = np.float64)

   # state 01
   dVij[1,1,0:100] = cj[0:100]
   # state 10
   dVij[2,2,100:] = cj[100:]
   # state 11
   dVij[3,3,:] = cj

   return dVij

def dHel0(R):
   ωj = parameters.ω
   return (ωj**2.0) * R

def initR():
   R0 = 0.0
   P0 = 0.0

   ω = parameters.ω
   T = 300 # Kelvin
   # Eh2J = 4.3597447222071e-18
   # β = (T*sc.k/Eh2J)**(-1) # 300 K == 1052.8
   β = 1052.583416013556
   sigP = np.sqrt( ω / ( 2 * np.tanh( 0.5*β*ω ) ) )
   sigR = sigP/ω
   
   # R = np.random.normal()*sigR + R0
   # P = np.random.normal()*sigP + P0

   """ Changed R and P sampling """
   R = np.random.normal(R0, sigR, size = parameters.ndof)
   P = np.random.normal(P0, sigP, size = parameters.ndof)

   return R, P