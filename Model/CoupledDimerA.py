#/usr/bin/python3

"""Model for 2 Level Coupled-Dimer"""
""" Model -> Ref.:https://doi.org/10.1063/1.5006824"""
""" Parameters -> Ref.:https://doi.org/10.1063/1.5031788 """

import numpy as np
from scipy import constants as sc

class parameters():
   ESteps = 20
   NSteps = 1000 #int(2*10**6)
   NTraj = 1000
   dtN = 1
   dtE = dtN/ESteps
   nskip = 1
   NStates = 4
   
   # cj and ωj is an array of dimensions (400,1)
   cj = np.loadtxt("./../Discretization/cj.txt")
   ωj = np.loadtxt("./../Discretization/ωj.txt")
   ndof = len(cj)

   """ Temporary derivation """
   # Tr[µ(t) µ(0) p(0)]
   # |0>, |1> , |2> , |3>
   # µ(0) = µ01(|0><1| + |1><0|)  + µ02(|0><2| + |2><0|) 
   # p(0) = |0><0| . pB
   # = TrB[   (TrE[µ(t) µ(0) |0><0|] -  TrE[µ(t) |0><0|  µ(0) ]).pB]
   # TrE[µ(t) (µ01|1><0| +  µ02|2><0|) ] - TrE[µ(t) (µ01|0><1| +  µ02|0><2|) ]

   initStatef = 2
   initStateb = 2

   Eh2c = 219474.6313632   # Hartree -> cm inverse
   
   M = 1          # """Needs to be confirmed"""
   ωc = 18/Eh2c   # 18 cm-1 in Hartree
   λ = 50/Eh2c    # 50 cm-1 in Hartree

def Hel(R):
   
   N = parameters.NStates
   cj = parameters.cj
   ωj = parameters.ωj
   ωc = parameters.ωc
   λ = parameters.λ
   m = parameters.M
   Eh2c = parameters.Eh2c

   # On-Site Energies:
   ε = np.array([0,-50,+50, 0], dtype = np.float64)
   ε_bar = 10000/Eh2c  # in cm-1

   VConst= np.zeros((N, N), dtype=np.float64)
   for i in range(len(VConst)):
      VConst[i,i] = ε_bar + (ε[i])/Eh2c
   VConst[-1,-1] = 2*ε_bar
   # Diabatic Couplings
   VConst[1,2] = VConst[2,1] = 100/Eh2c
   
   VCouple = np.zeros((N,N), dtype = np.float64)
 
   VCouple[1,1] = np.sum(cj[100:200]*R[100:200])
   VCouple[2,2] = np.sum(cj[200:300]*R[200:300])
   VCouple[3,3] = np.sum(cj[300:400]*R[300:400])

   VMat = VConst + VCouple #+ VBath

   # V0 = (VMat[0,0]+VMat[1,1])/2.
   # VMat[0,0] = VMat[0,0] - V0
   # VMat[1,1] = VMat[1,1] - V0
   return VMat

def  dHel(R):
   N = parameters.NStates
   cj = parameters.cj
   ωj = parameters.ωj
   ndof =  parameters.ndof
   dVij = np.zeros((N,N,ndof), dtype = np.float64)

   # state 00
   dVij[0,0,0:ndof] = 0 
   # state 01
   dVij[1,1,100:200] = cj[100:200]
   # state 10
   dVij[2,2,200:300] = cj[200:300]
   # state 11
   dVij[3,3,300:400] = cj[300:400]

   return dVij

def dHel0(R):
   ωj = parameters.ωj
   return (ωj**2.0) * R

def initR():
   R0 = 0.0
   P0 = 0.0

   ω = parameters.ωj
   T = 300 # Kelvin
   Eh2J = 4.3597447222071e-18
   β = (T*sc.k/Eh2J)**(-1) # 300 K == 1052.8
    
   sigP = np.sqrt( ω / ( 2 * np.tanh( 0.5*β*ω ) ) )
   sigR = sigP/ω
   
   # R = np.random.normal()*sigR + R0
   # P = np.random.normal()*sigP + P0

   """ Changed R and P sampling """
   R = np.random.normal(R0, sigR, size = parameters.ndof)
   P = np.random.normal(P0, sigP, size = parameters.ndof)

   return np.array([R]), np.array([P])