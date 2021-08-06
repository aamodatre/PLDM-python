#!/software/anaconda3/2020.07/bin/python
#SBATCH -p action 
#SBATCH -o output.log
#SBATCH --mem-per-cpu=1GB
#SBATCH -t 6:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
import sys, os
sys.path.append(os.popen("pwd").read().replace("\n",""))
sys.path.append(os.popen("pwd").read().replace("\n","")+"/Model")
#-------------------------
import pldm as method
import CoupledDimerA as model
# stype = 'sampled'
stype = "focused"
#-------------------------
from multiprocessing import Pool
import time 
import numpy as np

t0 = time.time()
#----------------
trajs = model.parameters.NTraj
#----------------
try:
    fold = sys.argv[1]
except:
    fold = "."
#------------------------------------------------------------------------------------------
#--------------------------- SBATCH -------------------------------------------------------
sbatch = [i for i in open('parallel.py',"r").readlines() if i[:10].find("#SBATCH") != -1 ]
cpu = int(sbatch[-1].split("=")[-1].replace("\n","")) 
nodes = int(sbatch[-2].split()[-1].replace("\n",""))
print (f"nodes : {nodes} | cpu : {cpu}")
procs = cpu * nodes  
ntraj = procs * trajs
print (f"Total trajectories {ntraj}")
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------

t1 = time.time()
with Pool(procs) as p:

    NSteps = model.parameters.NSteps
    NTraj = model.parameters.NTraj
    NStates = model.parameters.NStates

    #------ Arguments for each CPU------------------
    args = []
    for j in range(procs):
        par = model.parameters() 
        par.ID   = j
        par.SEED   = np.random.randint(0,100000000)
        
        #---- methods in model ------

        par.dHel = model.dHel
        par.dHel0 = model.dHel0
        par.initR = model.initR
        par.Hel   = model.Hel
        par.stype = stype

        args.append(par)
    #-----------------------------------------------
    print(f"Running : {par.NTraj*procs}  trajectories in {procs} cpu" )

    #------------------- parallelization -----------------------------------
    rho_ensemble  = p.map(method.runTraj, args)
    #-----------------------------------------------------------------------

#------------------- Gather --------------------------------------------
rho_sum = np.zeros(rho_ensemble[0].shape, dtype = rho_ensemble[0].dtype)
for i in range(procs):
    for t in range(rho_ensemble[0].shape[-1]):
        rho_sum[:,:,t] += rho_ensemble[i][:,:,t]


# PiiFile = open(f"{fold}/Pii.txt","w") 
# NTraj = model.parameters().NTraj
# for t in range(rho_ensemble[0].shape[-1]):
#     PiiFile.write(f"{t * model.parameters.nskip * model.parameters.dtE/2} \t")
#     for i in range(NStates):
#         PiiFile.write(str(rho_sum[i,i,t].real / (  procs * NTraj ) ) + "\t")
#     PiiFile.write("\n")
# PiiFile.close()

PijFile = open(f"{fold}/Pij.txt","w") 
NTraj = model.parameters().NTraj
for t in range(rho_ensemble[0].shape[-1]):
    PijFile.write(f"{t * model.parameters.nskip * model.parameters.dtE/2} \t")
    for i in range(NStates):
        for j in range(NStates):
                PijFile.write(str(rho_sum[i,j,t].real / (  procs * NTraj ) ) + "\t")
                PijFile.write(str(rho_sum[i,j,t].imag / (  procs * NTraj ) ) + "\t")
    PijFile.write("\n")
PijFile.close()

t2 = time.time()-t1
print(f"Total Time: {t2}")
print(f"Time per trajectory: {t2/ntraj}")
