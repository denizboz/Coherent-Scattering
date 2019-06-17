#This is the simulation of an electron incoming from the left to a nuclear region with n spins.
#For a range of different reflection amplitudes, it calculates the probability for the electron
#to depart the region from the left, and saves it to a generated text file.

import numpy as np
from matplotlib import pyplot as plt
import itertools as itt
import time

def totalref(nuclei,lenn,electron,npos,epos,branch,rcount,tcount,sumup):
    if((electron[0]==0 and electron[1]==-1) or electron[0]==lenn):
        if(sum(nuclei)==lenn/2+1):
            sumup.append([list(nuclei),rcount,tcount])   #save the output and r&t counts.
    else:
        if(electron[1]==1):                 #electron is up, inside somewhere
            if(nuclei[electron[0]]==1):     #if spins are parallel
                electron[0]+=1              #let the electron pass
                totalref(nuclei,lenn,electron,npos,epos,branch,rcount,tcount,sumup)
            elif(nuclei[electron[0]]==0):         #if spins are anti-parallel
                for i in range(lenn):
                    npos[branch][i]=nuclei[i]     #save the nucspin config
                for i in range(2):
                    epos[branch][i]=electron[i]   #save the electron config
                branch+=1
                electron[0]+=1      #electron transmits
                tcount += 1
                totalref(nuclei,lenn,electron,npos,epos,branch,rcount,tcount,sumup)
                tcount -= 1
                branch-=1
                for i in range(lenn):
                    nuclei[i]=npos[branch][i]   #turn back to first previous branch
                for i in range(2):
                    electron[i]=epos[branch][i]
                branch+=1
                nuclei[electron[0]]=1
                electron[1]=-1          #electron reflects, spin flips
                rcount += 1
                totalref(nuclei,lenn,electron,npos,epos,branch,rcount,tcount,sumup)
                rcount -= 1
                branch-=1
        elif(electron[1]==-1):              #electron is down, inside somewhere
            if(nuclei[electron[0]-1]==0):   #if spins are parallel
                electron[0]-=1              #let the electron pass
                totalref(nuclei,lenn,electron,npos,epos,branch,rcount,tcount,sumup)
            elif(nuclei[electron[0]-1]==1): #if the spins are anti-parallel
                for i in range(lenn):
                    npos[branch][i]=nuclei[i]     #save the nucspin config
                for i in range(2):
                    epos[branch][i]=electron[i]   #save the electron config
                branch+=1
                electron[0]-=1      #electron transmits
                tcount += 1
                totalref(nuclei,lenn,electron,npos,epos,branch,rcount,tcount,sumup)
                tcount -= 1
                branch-=1
                for i in range(lenn):
                    nuclei[i]=npos[branch][i]     #turn back to first previous branch
                for i in range(2):
                    electron[i]=epos[branch][i]
                branch+=1
                nuclei[electron[0]-1]=0
                electron[1]=1          #electron reflects, spin flips
                rcount += 1
                totalref(nuclei,lenn,electron,npos,epos,branch,rcount,tcount,sumup)
                rcount -= 1
                branch-=1
    return sumup            #return the nx2 matrix of outputs and rt multip's.

def SpinGenerator(n):
    index = []
    x = np.array(list(itt.product([0,1],repeat=n)))
    for i in range(len(x)):
        if(sum(x[i])!=n/2):
            index.append(i)
    x = np.delete(x,index,0)
    return x

def Bigtotal(n,r):
    t = (1-r**2)**0.5*1j
    start_time = time.clock()
    perms = SpinGenerator(n)
    size = len(perms)
    alltotal = 0
    for i in range(size):    #for each input
        electron = [0,1]
        branch = 0
        rcount = 0
        tcount = 0
        npos = np.zeros((50*n,n),dtype=int)
        epos = np.zeros((50*n,2),dtype=int)
        sumup = []
        sum_matrix = totalref(perms[i],n,electron,npos,epos,branch,rcount,tcount,sumup)
        sum_matrix = np.array(sum_matrix)
        calcMatrix = np.zeros((len(sum_matrix),2),dtype=np.ndarray)
        for i in range(len(sum_matrix)):
            calcMatrix[i,0] = sum_matrix[i,0]
            calcMatrix[i,1] = r**sum_matrix[i,1]*t**sum_matrix[i,2]
        for i in range(len(calcMatrix)):                #for the same outputs, sum the rt multip's and
            for j in range(i+1,len(calcMatrix)):        #gather them in one output
                if(calcMatrix[i,0]==calcMatrix[j,0]):
                    calcMatrix[i,1] += calcMatrix[j,1]
                    calcMatrix[j] = 0
        index = []
        for i in range(len(calcMatrix)):                #delete the residues (empty lines)
            for j in range(1):
                if(calcMatrix[i,j]==0) and (calcMatrix[i,j+1]==0):
                    index.append(i)
        calcMatrix = np.delete(calcMatrix,index,0)
        summ = 0
        for i in range(len(calcMatrix)):
            summ += abs(calcMatrix[i,1])**2
        alltotal += summ
    totalp = alltotal/size
    return totalp

rs = np.linspace(0,1,101)

start_time = time.clock()

n = 6       #number of nuclear spins

Pr = []
for r in rs:
    AA = Bigtotal(n,r)
    Pr.append(AA)
np.savetxt('RTnucspin' + str(n) + '.txt', Pr)

time_elapsed = time.clock() - start_time
print("Computation time:",time_elapsed)
