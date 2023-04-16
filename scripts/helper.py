import os
import numpy as np

def readfile(fname):
    script_dir = os.path.dirname(os.path.realpath(__file__))
    repo_dir = script_dir[:script_dir.rfind("/")]
    fpath = repo_dir + "/build/data/" + fname

    lines = []
    with open(fpath) as f:
        lines = f.readlines()

    dim = len(lines[0].split())
    npoints = len(lines)

    data = np.zeros((dim,npoints))

    i = 0
    for line in lines:
        e = 0
        for val in line.split():
            data[e][i] = val
            e += 1
        i += 1

    return data

def writeToFile(fname, data):
    script_dir = os.path.dirname(os.path.realpath(__file__))
    repo_dir = script_dir[:script_dir.rfind("/")]
    fpath = repo_dir + "/build/data/" + fname

    if os.path.isfile(fpath):
        check = input("Overwrite contents of " + fname + " ? [y/n] ")
        if check.lower() != "y":
            print(fname + " not overwritten")
            return

    f = open(fpath, "w")
    for line in data.T:
        for val in line:
            f.write(str(val) + " ")
        f.write("\n")


def plot_dir():
    return os.path.dirname(os.path.realpath(__file__)) + "/plots/"

def sotaiPhases(ax, textPos=0):
    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()
    ax.axvline(2.1, color='black', linestyle='dashed', linewidth=0.8)
    ax.axvline(3, color='black', linestyle='dashed', linewidth=0.8)
    ax.axvline(3.5, color='black', linestyle='dashed', linewidth=0.8)
    ax.axvline(6, color='black', linestyle='dashed', linewidth=0.8)
    if(textPos == 0):
        textPos = 0.9*ymax
    else:
        textPos *= ymax
    ax.text((xmin + 2.1)/2, textPos, r'I', ha = 'center', fontsize=8)
    ax.text((2.1 + 3)/2, textPos, r'II', ha = 'center', fontsize=8)
    ax.text((3 + 3.5)/2, textPos, r'III', ha = 'center', fontsize=8)
    ax.text((3.5 + 6)/2, textPos, r'IV', ha = 'center', fontsize=8)
    ax.text((6 + xmax)/2, textPos, r'V', ha = 'center', fontsize=8)

def totaiPhases(ax, textPos=0):
    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()
    ax.axvline(2.55, color='black', linestyle='dashed', linewidth=0.8)
    ax.axvline(3.54, color='black', linestyle='dashed', linewidth=0.8)
    if(textPos == 0):
        textPos = 0.9*ymax
    else:
        textPos *= ymax
    ax.text((xmin + 2.55)/2, textPos, r'I', ha = 'center', fontsize=8)
    ax.text((2.55 + 3.54)/2, textPos, r'II', ha = 'center', fontsize=8)
    ax.text((3.54 + xmax)/2, textPos, r'III', ha = 'center', fontsize=8)
