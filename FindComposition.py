#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 12:40:49 2022

@author: andrew
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

numbin = 20     #number of bins/slices of sphere
ext =""
#file="CoreShell5050_relaxedAuNi.stru"
file="10.stru"

data = pd.read_csv(file, skiprows=3)
data = data.drop(data.columns[3:], axis=1)
data[["atoms", "x"]] = data[data.columns[0]].str.split(expand=True)
data = data.drop(data.columns[0], axis=1)
data = data.set_axis(["y", "z", "atom", "x"], axis=1, inplace=False)
data["x"] = data["x"].astype(float)
data["y"] = data["y"].astype(float)
data["z"] = data["z"].astype(float)
data["distance"] = np.sqrt((data["x"]*data["x"] + data["y"]*data["y"] + data["z"]*data["z"])*(4.0789*4.0789))
print("Max distance in angstroms "+str(data["distance"].max()))
print("Min distance in angstroms "+str(data["distance"].min()))
realBins = np.linspace(0,data["distance"].max(),numbin)
data["bin"] = pd.cut(data["distance"], bins=realBins, right=False, precision=2,include_lowest='True')

data2 = data.groupby("bin")["atom"].value_counts(normalize=True).unstack()
data2 = data2.fillna(0)

fig, ax = plt.subplots(figsize=(7,5))
for key in data2.keys():
    data3 = data2[key].reset_index()
    ax.scatter(range(len(data3['bin'])),data3[key])

#ax = data.plot(figsize=(7,5), kind='bar')
x = []
for c in data3.bin:
    x.append(str(c))

ax.set_title("Composition of Relaxed AuNi", fontsize=20)
ax.set_ylabel("Fraction", fontsize=16)
ax.tick_params(axis='y', labelsize=12)
ax.set_xticks(np.arange(len(x)))
ax.set_xticklabels(x)
ax.set_xlabel("Distance [A]", fontsize=16)
ax.tick_params(axis='x', labelsize=10, labelrotation = 70)
ax.legend()
ax.grid(True)
plt.tight_layout()
plt.savefig(ext+"composition.pdf")



