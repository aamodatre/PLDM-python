#/usr/bin/python3
"""Model for 2D Bilinear Spectra"""
"""Ref.:https://doi.org/10.1063/1.5031788"""
import numpy as np


def model(M=2):

   """Units in the following table are in cm-1"""
   #              |  M0 |  M1 |
   eps = np.array([500.0, 400.0])
   delta = np.array([100.0, 100.0])


   return eps[M], delta[M]

# From Spin-Boson
class parameters():
   NSteps = 200 #int(2*10**6)
   NTraj = 200
   dtN = 0.01       
   dtE = dtN/20
   NStates = 2
   M = 1    
   initStateb = 0
   initStatef = 1
   nskip = 10
   eps, delta = model(1) # model3


""" Questions
1. What is the NSteps variable. 
2. Set appropriate value for dtN
"""