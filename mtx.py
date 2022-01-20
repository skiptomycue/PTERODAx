import json
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm

path = 'COVX/mtx/'
plt.rcParams.update({'figure.max_open_warning': 60})
groups = 44

def cvxRoll(oldpath, path):

    for root, dirs, files in os.walk(oldpath, topdown=False):

        for name in sorted(files):

            roll = open(oldpath + name, 'r')

            lines = roll.readlines()

            k = -1
            l = 0

            srot = []

            while l + k + 1 < len(lines):

                if lines[l][2] in ['v', 'd']:
                    srot.append(lines[l].replace('\n', ''))
                    k += 1

                else:
                    srot[k] += str(lines[l]).replace('\n', '')

                l += 1

            textfile = open(path + name, "w")
            for element in srot:
                textfile.write(element + "\n")
            textfile.close()
#cvxRoll('COVX/OLD/', path)

def buildMTX(lines, l, z):

    words = lines[l].lower().split()

    ctrl = [int(a) for a in lines[l + 1].lower().split()[1:(groups * 2 + 1)]]
    w = np.array(ctrl).reshape(int(len(ctrl) / 2), 2).transpose()[0].tolist()
    b = np.array(ctrl).reshape(int(len(ctrl) / 2), 2).transpose()[1].tolist()

    vals = [float(a) for a in lines[l + 2].lower().split()[1:]]

    cvx[z][words[2]] = np.zeros((groups, groups)).tolist()

    k = 0

    if sum(w) == len(vals):

        for i in range(groups):
            for j in range(groups):
                if j >= (i - b[i] + 1) and j < (w[i] + i - b[i] + 1):
                    cvx[z][words[2]][i][j] = vals[k]
                    k += 1

        plotCovx(cvx[z][words[2]], [z, words[2]])

    else:

        print(z+'-'+words[2]+' is corrupt')
def checkMTX(lines, l, z):

    words = lines[l].lower().split()

    if words[1] in za_spesa and words[2] in MT_spesa and words[1] == words[3] and words[2] == words[4]:

        buildMTX(lines, l, z)

    elif words[1] == '9228' and words[2] in MT_spesa and words[1] == words[3] and words[2] == words[4]:

        buildMTX(lines, l, '92235')

    elif words[1] == '9237' and words[2] in MT_spesa and words[1] == words[3] and words[2] == words[4]:

        buildMTX(lines, l, '92238')
def plotCovx(A, name):

    fig, axs = plt.subplots()

    im = axs.imshow(A, cmap='RdYlGn', norm=SymLogNorm(linthresh=1E-10, base=10))
    fig.colorbar(im, orientation='vertical')
    axs.set(xlabel='Energy groups', title='Covariance matrix for '+name[0]+' MT='+name[1])

    fig.savefig('COVX/PLOT/'+name[0]+'-'+name[1]+'.png')
def shop():

    for z in za_spesa:

            roll = open(path+'cvx.mat'+z, 'r')
            lines = roll.readlines()
            cvx[z]={}

            for l in range(len(lines)):

                if '7d' in lines[l][0:3]:

                    checkMTX(lines, l, z)

cvx = {}

za_spesa = ['92235', '92238', '94239', '01001', '08016']
MT_spesa = ['452', '18', '102', '2', '4']
MT_names = ['nubar', 'fission', 'capture', 'elastic', 'inelastic']

shop()

with open('COVX/mtx.json', 'w') as fp:
    json.dump(cvx, fp)
