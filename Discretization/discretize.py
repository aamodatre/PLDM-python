import numpy as np
from numpy import pi as π
import matplotlib.pyplot as plt

def J(ω):
    ωc = 0.000082014 # 18 cm-1 in Hartree
    λ = 0.00022782 # 50 cm-1 in Hartree
    return 2*λ*((ω/ωc)/ (1+(ω/ωc)**2))

def getParameters(N, Jω):
    
    # Spectral Density Parameters
    ωc = 0.000082014 # 18 cm-1 in Hartree
    λ = 0.00022782 # 50 cm-1 in Hartree

    # Frequencies considered
    ω = np.linspace(0.00000001,250*ωc,30000) # upper limit - reconsider
    dω = ω[1] - ω[0]
    Jω = Jω(ω)

    Fω = np.zeros((len(ω)))
    
    # Integral as per ipynb file:
    for iω in range(len(ω)):
        Fω[iω] = (1.0/π) * np.sum( Jω[:iω]/ω[:iω] ) * dω
    
    λs = Fω[-1] 

    print("\nReorganization Energy : {}\n".format(λs))

    ωj = np.zeros((N))
    cj = np.zeros((N))
    
    for i in range(N):
        j = i+1
        ωj[i] = ω[np.argmin(np.abs(Fω - ((j-0.5)/N) *  λs))] 
        cj[i] = 2.0 * ωj[i] * (λs/(2.0 * N))**0.5 

    # dj = 0.5 * cj**2.0 / ωj**3
    
    fig = plt.figure(figsize = (10,10))
    ax = fig.add_subplot(1,1,1)
    ax.set_title("${J(\omega)}/{\omega}$", fontsize = 15)
    ax.set_xlabel("$\omega$", fontsize = 15)
    ax.plot(ω,Jω/ω)
    plt.show()
    return ωj, cj

# Degrees of Freedom
N = 100

ωj, cj = getParameters(N,J)

ω = np.zeros((400))
c = np.zeros((400))

ω[:100] = ωj
ω[100:200] = ωj
ω[200:300] = ωj
ω[300:400] = ωj

c[:100] = 0
c[100:200] =  cj
c[200:300] =  cj
c[300:400] =  2 * cj

np.savetxt("ωj.txt", ω)
np.savetxt("cj.txt", c)