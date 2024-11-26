# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 14:26:21 2024

@author: bwehlin
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import collections  as mc

def plotDomain(xPath, yPath, ax=None, pointsize=0.1, linewidth=0.4):
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 3))
        
    X = np.loadtxt(xPath, delimiter=' ', comments=None, skiprows=0)
    Y = np.loadtxt(yPath, delimiter=' ',  comments=None, skiprows=0)

    print("X:", X)
    print("Y:", Y)

    
    nI = X.shape[0] - 2
    nJ = X.shape[1] - 2
    
    hlines = [[(X[i, j], Y[i, j]), (X[i+1, j], Y[i+1, j])] for i in range(nI+1) for j in range(nJ+2)]
    vlines = [[(X[i, j], Y[i, j]), (X[i, j+1], Y[i, j+1])] for i in range(nI+2) for j in range(nJ+1)]
    
    # X = X.reshape(-1, 1)
    # Y = Y.reshape(-1, 1)
    X = load_clean_grid(xPath)
    Y = load_clean_grid(yPath)

    
    ax.scatter(X, Y, s=pointsize)

    hlc = mc.LineCollection(hlines, linewidths=linewidth)
    vlc = mc.LineCollection(vlines, linewidths=linewidth)
    ax.add_collection(hlc)
    ax.add_collection(vlc)

    ax.axis('equal')
    
    return ax

xPath = 'xGrid.txt'
yPath = 'yGrid.txt'

with open(xPath, 'r') as file:
    print("xGrid.txt:")
    for i, line in enumerate(file):
        print(f"Row {i}: {repr(line)}")

with open(yPath, 'r') as file:
    print("\nyGrid.txt:")
    for i, line in enumerate(file):
        print(f"Row {i}: {repr(line)}")

def clean_file(filepath):
    with open(filepath, 'r') as file:
        cleaned_lines = [line.strip() for line in file if line.strip()]  # Remove empty lines and strip spaces
    with open(filepath, 'w') as file:
        file.write('\n'.join(cleaned_lines))

clean_file('xGrid.txt')
clean_file('yGrid.txt')

def load_clean_grid(filepath):
    with open(filepath, 'r') as file:
        return np.array([
            list(map(float, line.split()))  # Convert each space-separated value to float
            for line in file if line.strip()  # Skip empty lines
        ])

ax = plotDomain(xPath, yPath)

plt.show()
