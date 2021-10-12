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

def plot_dir():
    return os.path.dirname(os.path.realpath(__file__)) + "/plots/"
