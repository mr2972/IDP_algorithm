#!/bin/python3
#@title
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.style.use('classic')
plt.rcParams['xtick.major.size'] = 8
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.size'] = 8
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['axes.linewidth'] = 2


amino_acids = input("Enter the FASTA sequence: ")
amino_acids = list(amino_acids)
total_aa = len(amino_acids)  # total number of amino acids
namex = ['K', 'R', 'E', 'D']
namep = ['K', 'R']
namen = ['E', 'D']
namenonpolar = ['V', 'I', 'L', 'M', 'F', 'Y', 'W']
namepolar = ['S', 'T', 'N', 'Q']
others = ['C', 'A', 'H']
gly_pro = ['G', 'P']

print(f'******')
print(f'The Default Binsize is 10 and cutoff is 25')
print(f'In case you want to change it write the number in respective prompts')
print(f'******')
BINSIZE = input("Enter Binsize: ") or 10
BINSIZE = int(BINSIZE)

cutoff = input("Enter Cutoff: ") or 25
cutoff = int(cutoff)

i = 0
charged_aa = []
polar_aa = []
hydrophobic_aa = []
netcharge_aa = []
others_aa = []
gly_pro_aa = []

while i <= total_aa - 1:
    j = i
    positivex = 0
    negativex = 0
    polarx = 0
    nonpolarx = 0
    total_chargedx = 0
    net_chargedx = 0
    othersx = 0
    gly_prox = 0

    while j < i + BINSIZE:
        if j <= total_aa - 1:
            if amino_acids[j] in namex:
                total_chargedx += 1
                if amino_acids[j] in namep:
                    positivex += 1
                if amino_acids[j] in namen:
                    negativex += 1
            elif amino_acids[j] in namepolar:
                polarx += 1
            elif amino_acids[j] in namenonpolar:
                nonpolarx += 1
            elif amino_acids[j] in gly_pro:
                gly_prox += 1
            elif amino_acids[j] in others:
                othersx += 1
        j += 1

    if j > total_aa:
        BX = total_aa % BINSIZE
    else:
        BX = BINSIZE

    charged_aa.append(total_chargedx / BX)
    polar_aa.append(polarx / BX)
    hydrophobic_aa.append(nonpolarx / BX)
    net_chargedx = positivex - negativex
    netcharge_aa.append(net_chargedx / BX)
    others_aa.append(othersx / BX)
    gly_pro_aa.append(gly_prox / BX)
    i += BINSIZE

length = len(charged_aa)
lenh = 0
start = 0
end = BINSIZE
size_plot = total_aa // BINSIZE

# Cell2
dict_ = {
    "START": [],
    "END": [],
    "Hydrophobe": [],
    "Polar": [],
    "Charged": [],
    "Others": [],
    "G-P": []
}

i = 0
while i < len(charged_aa):
    dict_["START"].append(i * BINSIZE)
    dict_["END"].append((i + 1) * BINSIZE)
    dict_["Hydrophobe"].append(hydrophobic_aa[i] * 100)
    dict_["Polar"].append(polar_aa[i] * 100)
    dict_["Charged"].append(charged_aa[i] * 100)
    dict_["Others"].append(others_aa[i] * 100)
    dict_["G-P"].append(gly_pro_aa[i] * 100)
    i += 1

dict_ = pd.DataFrame(dict_)

# Cell3
length = dict_.loc[:, "START"].values[:]
length = len(length)
i = 0
ID_region1 = {}

while i < length:
    if int(dict_.loc[i, "Hydrophobe"]) < int(cutoff):
        ID_region1[(i) * BINSIZE] = 1  # is IDP
    else:
        ID_region1[(i) * BINSIZE] = 0  # is not IDP
    i += 1

keys1 = []
for keys in ID_region1.keys():
    keys1.append(keys)

i = 0
ID_region2 = {}
while i < len(keys1):
    if i == 0:
        if ID_region1[(i) * BINSIZE] == 1:
            ID_region2[(i) * BINSIZE] = 1  # is IDP
        elif ID_region1[(i) * BINSIZE] == 0:
            if ID_region1[(i) * BINSIZE] == 0 and ID_region1[(i + 1) * BINSIZE] == 1:
                ID_region2[(i) * BINSIZE] = 1  # is IDP
                ID_region1[0] = 1  # Updating for first bin in ID_region1
    elif i < len(keys1) - 1:
        if (
            ID_region1[(i) * BINSIZE] == 1
            and ID_region1[(i + 1) * BINSIZE] == 1
        ):
            ID_region2[(i) * BINSIZE] = 1  # is IDP
        elif (
            ID_region1[(i - 1) * BINSIZE] == 1
            and ID_region1[(i) * BINSIZE] == 1
        ):
            ID_region2[(i) * BINSIZE] = 1  # is IDP
        if (
            ID_region1[(i) * BINSIZE] == 0
            and ID_region1[(i + 1) * BINSIZE] == 1
            and ID_region1[(i - 1) * BINSIZE] == 1
        ):
            ID_region2[(i) * BINSIZE] = 1  # is IDP

        value1 = (i + 2) * BINSIZE
        if value1 in ID_region1.keys():
            if (
                ID_region1[(i) * BINSIZE] == 1
                and ID_region1[(i + 1) * BINSIZE] == 0
                and ID_region1[(i + 2) * BINSIZE] == 1
            ):
                ID_region2[(i) * BINSIZE] = 1  # is IDP

        value2 = (i - 2) * BINSIZE
        if value2 in ID_region1.keys():
            if (
                ID_region1[(i) * BINSIZE] == 1
                and ID_region1[(i - 1) * BINSIZE] == 0
                and ID_region1[(i - 2) * BINSIZE] == 1
            ):
                ID_region2[(i) * BINSIZE] = 1  # is IDP

    elif i == len(keys1) - 1:
        if ID_region1[(i) * BINSIZE] == 1:
            ID_region2[(i) * BINSIZE] = 1  # is IDP
        elif (
            ID_region1[(i) * BINSIZE] == 0
            and ID_region1[(i - 1) * BINSIZE] == 1
        ):
            ID_region2[(i) * BINSIZE] = 1  # is IDP

    i += 1

disorder_value = []
i = 0
while i < total_aa:
    disorder_value.append(0)
    i += 1

Disorder_region = []
for keys in ID_region2.keys():
    start = int(keys)
    end = start + BINSIZE
    r = str(start) + "-" + str(end)
    Disorder_region.append(r)

    span = int(keys) + BINSIZE
    if span > total_aa - 1:
        span = total_aa - 1

    while start <= span:
        disorder_value[start] = 1
        start += 1

print("******************************************")
print("Total amino acids= %d" % total_aa)
print("DISORDER REGION: %s " % (Disorder_region))
print("******************************************")

plt.figure(figsize=(12,7),facecolor='white')
plt.plot(np.arange(1, total_aa + 1, 1), disorder_value)
plt.xlabel("Residue Number", fontsize=14, fontweight="bold")
plt.ylabel("Disorderness", fontsize=14, fontweight="bold")
plt.xlim(1, total_aa)
plt.ylim(-0.1, 1.1)
plt.yticks(np.arange(0,2,1), fontsize=18)
plt.savefig("disorder.png")
plt.show()
