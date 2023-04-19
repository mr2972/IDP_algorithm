#!/bin/python3
import sys
import numpy as np
import pandas as pd
import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use('classic')
matplotlib.use('pdf')
mpl.rcParams['xtick.major.size'] = 8
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['ytick.major.size'] = 8
mpl.rcParams['ytick.major.width'] = 2
mpl.rcParams['axes.linewidth']= 2
#picture=input("Enter Disprot ID: ")
#BINSIZE=input("Input BINSIZE= ")
#BINSIZE=int(BINSIZE)
BINSIZE=10

data= sys.argv[1]
#print(data)
aa=['K','R','H','E','D','G','A','V','C','I','L','M','F','Y','W','P','S','T','N','Q']
amino_acids=[]
for line in data:
    if line.startswith('>')!=True:
        #print(line)
        for alphabet in line:
            #print(alphabet)
            if alphabet.upper() in aa:
                entry=alphabet.upper()
                amino_acids.append(entry)
#print(amino_acids)
total_aa=len(amino_acids)     ##total number of amino acids
namex=['K','R','E','D']
namep=['K','R']
namen=['E','D']
namenonpolar=['V','I','L','M','F','Y','W',]
namepolar=['S','T','N','Q']
others=['C','A','H']
gly_pro=['G','P']
##*****************************************************
#print(amino_acids)
i=0
charged_aa=[]
polar_aa=[]
hydrophobic_aa=[]
netcharge_aa=[]
others_aa=[]
gly_pro_aa=[]
while i <= total_aa-1:
    j=i
    positivex=0
    negativex=0
    polarx=0
    nonpolarx=0
    total_chargedx=0
    net_chargedx=0
    othersx=0
    gly_prox=0
    #print(i,j)
    while j<i+BINSIZE:
        if j<= total_aa-1:
            if amino_acids[j] in namex:
                total_chargedx+=1
                if amino_acids[j] in namep:
                    positivex+=1
                if amino_acids[j] in namen:
                    negativex+=1
            elif amino_acids[j] in namepolar:
                polarx+=1
            elif amino_acids[j] in namenonpolar:
                nonpolarx+=1
            elif amino_acids[j] in gly_pro:
                gly_prox+=1
            elif amino_acids[j] in others:
                othersx+=1
        j+=1
    if j>total_aa:
        BX=total_aa%BINSIZE
    else:
        BX=BINSIZE
    charged_aa.append(total_chargedx/BX)
    polar_aa.append(polarx/BX)
    hydrophobic_aa.append(nonpolarx/BX)
    net_chargedx=positivex-negativex
    netcharge_aa.append(net_chargedx/BX)
    others_aa.append(othersx/BX)
    gly_pro_aa.append(gly_prox/BX)
    i+=BINSIZE
length=len(charged_aa)
#print(length,charged_aa,total_aa)
lenh=0
start=0
end=BINSIZE
size_plot=total_aa//BINSIZE
#result=pd.read_csv("table.csv", index_col=0)
#result = pd.DataFrame(result.loc[result['ID'] == picture])
##Cell2
#print("Start,END,hydrophobic_aa[i],polar_aa[i],charged_aa[i],others_aa[i],gly_pro_aa[i]")
dict_={"START":[],"END":[],"Hydrophobe":[],"Polar":[],"Charged":[],"Others":[],"G-P":[]}
i=0
while i<len(charged_aa):
    #print((i)*BINSIZE,(i+1)*BINSIZE,hydrophobic_aa[i]*100,polar_aa[i]*100,charged_aa[i]*100,others_aa[i]*100,gly_pro_aa[i]*100)
    dict_["START"].append((i)*BINSIZE)
    dict_["END"].append((i+1)*BINSIZE)
    dict_["Hydrophobe"].append(hydrophobic_aa[i]*100)
    dict_["Polar"].append(polar_aa[i]*100)
    dict_["Charged"].append(charged_aa[i]*100)
    dict_["Others"].append(others_aa[i]*100)
    dict_["G-P"].append(gly_pro_aa[i]*100)
    i+=1
dict_=pd.DataFrame(dict_)
dict_

##Cell3
length=dict_.loc[:,"START"]
length=len(length)
i=0
ID_region1={}
while i<length:
    if dict_.loc[i,"Hydrophobe"]<25:
        ID_region1[(i)*10]=1 # is IDP
    else:
        ID_region1[(i)*10]=0 #is not IDP
    i+=1
#print(ID_region1)
#**
keys1=[]
for keys in ID_region1.keys():
    keys1.append(keys)
i=0
ID_region2={}
while i<len(keys1):
    if i==0:
        if ID_region1[(i)*10]==1 and ID_region1[(i+1)*10]==1:
            ID_region2[(i)*10]=1 # is IDP
    elif i<len(keys1)-1:
        if ID_region1[(i)*10]==1 and ID_region1[(i+1)*10]==1:
            ID_region2[(i)*10]=1 # is IDP
        elif ID_region1[(i-1)*10]==1 and ID_region1[(i)*10]==1:
            ID_region2[(i)*10]=1 # is IDP
        if ID_region1[(i)*10]==0 and ID_region1[(i+1)*10]==1 and ID_region1[(i-1)*10]==1:
            ID_region2[(i)*10]=1 # is IDP
        value1=(i+2)*10
        if value1 in ID_region1.keys():
            if ID_region1[(i)*10]==1 and ID_region1[(i+1)*10]==0 and ID_region1[(i+2)*10]==1:
                ID_region2[(i)*10]=1 # is IDP
        value2=(i-2)*10
        if value2 in ID_region1.keys():
            if ID_region1[(i)*10]==1 and ID_region1[(i-1)*10]==0 and ID_region1[(i-2)*10]==1:
                ID_region2[(i)*10]=1 # is IDP
    elif i==len(keys1):
        if ID_region1[(i-1)*10]==1 and ID_region1[(i)*10]==1:
            ID_region2[(i)*10]=1 # is IDP
    i+=1


disorder_value=[]
i=0
while i<total_aa:
    disorder_value.append(0)
    i+=1

Disorder_region=[]
for keys in ID_region2.keys():
    start=int(keys)
    end=start+10
    r=str(start)+"-"+str(end)
    Disorder_region.append(r)
    while start<=int(keys)+10:
        disorder_value[start]=1
        start+=1
print("******************************************")
print("Total amino acids= %d"%total_aa)
print("DISORDER REGION: %s "%(Disorder_region))
print("******************************************")  
plt.figure(facecolor='white')
plt.plot(np.arange(1,total_aa+1,1),disorder_value)
plt.xlabel("Residue Number", fontsize=14, fontweight="bold")
plt.ylabel("Disorderness", fontsize=14, fontweight="bold")
plt.xlim(1,total_aa)
plt.ylim(-0.1,1.1)
plt.savefig("disorder.png")
#plt.show()
