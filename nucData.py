import numpy as np
import json
import math
import serpentTools as ST
from serpentTools.settings import rc
from periodictable import elements
import matplotlib.pyplot as plt
import serpent

model    = 'UO2-pin'                                # INPUT MODEL
energy   =  2                                       # INPUT ENERGY GROUPS
PASSI    =  20                                      # INPUT STEP NUMBER
fpSwitch =  False                                   # SWITCH TO FULL NUCLIDE CHART

### INITS ###

file     = model+'/'+str(energy)+'_groups'
res = ST.read(file + '/REP_res.m')
dep = ST.read(file + '/REP_dep.m')
keff = ST.read(file+'/REP_res.m').resdata['absKeff'][:,0]
time =  dep.days
giorni = dep.days[-1]
T = dep.days[-1]*24*3600
MXT=ST.MicroXSTuple

### TIME STEPS ###

tempo_het  = np.linspace(0,5,10).tolist() + np.linspace(5,100,20).tolist()[1:] + np.linspace(100,giorni,100).tolist()[1:]
tempo_homo = np.linspace(0, giorni, PASSI)
tempo = tempo_homo
steps = len(tempo)
nodo = min([int(len(dep.days)/2),int(len(tempo)/2)])

### SERPENT ###

REACTIONS=serpent.REACTIONS
vol=serpent.volumes[model]
VOL=[sum(a) for a in vol]
REG=serpent.regions[model]
UNI=serpent.unis[model]
DET=serpent.detectors[model]
mat=serpent.materials[model]

fuelId=REG.index('Fuel')
chi  = res.universes[UNI[fuelId],0.0, 0, 0.0].infExp['infChit']
grid = res.universes[UNI[fuelId],0.0, 0, 0.0].groups[::-1]
ene  = len(chi)


### FUNCTIONS ###

def getMM(zai):

    z=zai[:-4]
    a=zai[-4:-1]

    if int(a) in elements[int(z)].isotopes:

        MM=elements[int(z)][int(a)].mass

    elif zai == 60000:

        MM = 12

    else:

        MM = int(a)

    return MM

def getName(zai):

    z=zai[:-4]
    a=zai[-4:-1]

    if int(a) in elements[int(z)].isotopes:

        Name=elements[int(z)][int(a)].name

    elif zai == 60000:

        Name = 'carbon'

    else:

        Name = 'sarcazzo'

    Name += '-'+str(a)
    Name=Name[0].upper()+Name[1:]

    return Name

def xsPlot(sigma, MT, zai, t):

    x = grid[:-1]
    id= ZAI.index(zai)
    name = nuc[id].name
    fig, axs = plt.subplots()

    y = [sigma[MT][e][t][id]*1E+24/x[e] for e in range(ene)]
    axs.step(x, y, where = 'pre')
    tit = name+ ' MT '+MT+' cross section\n'
    axs.set(xlabel='Energy [MeV]', ylabel='cross section [barn/MeV]', title=tit)

    axs.set_yscale('log')
    axs.set_xscale('log')

    fig.savefig(model+'xs/'+zai+'_'+MT+'_xs.png')

def rescaleTime(xs, xsTime, newTime):

    newXs = []

    for t in newTime:

        j = 0

        while xsTime[j]+1 < t and j < (len(xsTime)-2):

            j+=1

        width = time[j+1] - time[j]

        #newXs.append(xs[j]*(xsTime[j+1]-t)/width+xs[j+1]*(t-xsTime[j])/width)
        newXs.append(xs[nodo])

    return newXs

def readMDEP():

    print('\nScanning nuclear data')

    mdep = []

    for i in range(len(time)):

        mdep.append(ST.read(file+'/REP_mdx'+str(i)+'.m'))

    return mdep

def getZai(mdep):

    out = []

    if fpSwitch == True:

        for u in UNI:

            zai = []

            for key in mdep[0].xsVal[u].keys():

                if key[0] not in [60000, 300000] and key[1] == 101:

                    zai.append(key[0])

            if u == UNI[fuelId]:

                zai.extend(serpent.sama)

            aiz=sorted(zai)

            fuori=[str(z) for z in aiz]

            out.append(fuori)

    else:

        out=serpent.zais[model]

    ZAI=sum(out,[])

    return out, ZAI

### OBJECTS ###

class xs:

    def __init__(self, **kwargs):

        keys = ['fis', 'cap', 'sca', 'nu']

        for key in keys:

            self.__setattr__(key, np.zeros(ene))

            if key in kwargs:
                self.__setattr__(key, kwargs[key])

        keys = ['gamma', 'lam', 'v']

        for key in keys:

            self.__setattr__(key, 0)

            if key in kwargs:
                self.__setattr__(key, kwargs[key])                      # cla                       # xs

def xsDef(**kwargs):

    keys = ['fis', 'cap', 'sca', 'nu']

    xs = {}

    for key in keys:

        xs[key] = np.zeros(ene)

        if key in kwargs:
            xs[key] = kwargs[key]

    keys = ['gam', 'lam', 'v']

    for key in keys:

        xs[key] = 0

        if key in kwargs:
            xs[key] = kwargs[key]

    return xs

class Nuclide:

    def __init__(self, name, idReg, zai, mat, vol, xs, nuc, **kwargs):

        if zai in xs[UNI[idReg]].keys() :

            for key in xs[UNI[idReg]][zai].keys():

                xs[UNI[idReg]][zai][key] = (np.array(xs[UNI[idReg]][zai][key]) / vol * VOL[idReg]).tolist()

            self.xs = {**xs[UNI[idReg]][zai], **xsDef(**kwargs)}

        else :

            self.xs = xsDef(**kwargs)


        if name[-1] in ['0','1']:

            self.name = getName(zai)

        else:

            self.name = name

        self.id = len(nuc)
        self.reg = idReg
        self.zai = zai
        self.mat = mat

        self.dens = 0
        self.comp = 1E-40
        self.vol  = 0
        self.at   = 0

        if int(zai[-1]) == 0:

            for i in range(len(mat)):

                if dep.materials[mat[i]].getValues('days','mdens', zai=int(zai))[0][0] > 0:

                    self.dens += dep.materials[mat[i]].getValues('days','mdens', zai=int(zai))[0][0]
                    self.vol  += vol[i]
                    self.at   += dep.materials[mat[i]].getValues('days','mdens', zai=int(zai))[0][0] / getMM(zai) * vol[i] * 6.022E+23
                    self.comp += (dep.materials[mat[i]].getValues('days','adens', zai=int(zai))[0][0]+1E-40) * vol[i] / VOL[idReg] * 6.022E+23

### MAIN ###

def buildScatt():

    scatt=[]
    for i in range(len(REG)):
        scatt.append(res.universes[UNI[i],0.0, 0, 0.0].infExp['infS0'].reshape(ene,ene)*VOL[i])

    return scatt

def buildAlbe(det):

    dire=['in','out']

    R=np.ones((2,2,ene))

    zero = np.zeros(ene)

    for i in [0,1]:
        for j in [0,1]:

            R[i][j]=abs(det.detectors[str(dire[i])+str(dire[j])].tallies[::-1])


    R_out=(abs(det.detectors['bordo'].tallies[::-1]))

    RR=np.array([[R[0,1], -R[0,0], zero],[-R[0,1], (R[0,0]+R[1,1]), -R[1,0]],[zero, -R[1,1], R[1,0]+R_out]])

    F=[]

    reg=serpent.detectors[model]

    for i in range(len(reg)):

        #F.append((det.detectors[reg[i]].tallies[::-1]+np.ones(ene)*1E-5)/VOL[i])
        F.append((det.detectors[reg[i]].tallies[::-1])/VOL[i])

    S=(RR/F).transpose()

    SS=np.array([S[j].transpose() for j in range(ene)])

    F=np.array(F).transpose().ravel()

    return SS, F

def cycleAlbe(time):

    Albe=[]
    Flux=[]


    for t in range(len(time)):

        det = ST.read(file+'/REP_det'+str(int(t))+'.m')
        Albedo, Flusso, = buildAlbe(det)
        Albe.append(Albedo)
        Flux.append(Flusso)

    return Albe, Flux,

def buildXS(xsTime, newTime, ZAI, mdep):

    xs={}

    for u in UNI:

        xs[u]={}

        for z in ZAI:

            xs[u][z]={}

            for r in REACTIONS:

                if MXT(zai=int(z),mt=int(r),metastable=0) in mdep[0].xsVal[u].keys() or MXT(zai=int(z),mt=int(r),metastable=1) in mdep[0].xsVal[u].keys() :

                    xs[u][z][r]=[]

                    for i in range(len(mdep)):

                        xs[u][z][r].append(mdep[i].getXS(universe=u, isotope=int(z), reaction=int(r))[0]*1E-24)

                    xs[u][z][r]=rescaleTime(xs[u][z][r], xsTime, newTime)

                else:

                    xs[u][z][r]=np.zeros((len(newTime),ene)).tolist()
                    #xs[u][z][r]=np.zeros((len(xsTime),2)).tolist()

        #xs[u]['10020']['101'] = ( np.array(xs[u]['80160']['101'])/2 + np.array(xs[u]['10020']['101'])).tolist()


    return xs

def buildNuc(zai, xs):

    nuc = []
    ZAI = []
    MAT = []

    nuU5 = np.array([2.65] * int((ene / 2)) + [2.43] * int((ene / 2)))
    nuPu = np.array([3.16] * int((ene / 2)) + [2.86] * int((ene / 2)))

    vU5 = 200.7E+6 * 1.6E-19  # J/fiss
    vPu = 207.0E+6 * 1.6E-19

    for i in range(len(UNI)):

        if i == fuelId:

            for z in zai[i]:

                nu = nuU5
                v  = vU5

                if z == '942390':
                    nu = nuPu
                    v  = vPu

                elem = Nuclide(z, i, z, mat[i], vol[i], xs, nuc, nu=nu, v=v)
                nuc.append(elem)
                ZAI.append(elem.zai)
                MAT.append(elem.mat)
        else:

            for z in zai[i]:
                elem = Nuclide(z, i, z, mat[i], vol[i], xs, nuc)
                nuc.append(elem)
                ZAI.append(elem.zai)
                MAT.append(elem.mat)


    return nuc, ZAI, MAT

def buildSig(nuc):

    sig = {}

    for key in REACTIONS+['sca']:

        sig[key]=[]

        for nuclide in nuc:

            if key == 'sca':
                sig[key].append((np.array(nuclide.xs[key]) / np.array(nuclide.comp)).tolist())
                #sig[key].append((np.array(nuclide.xs[key]) / np.array(nuclide.comp)).tolist()*steps)

            else:
                sig[key].append(nuclide.xs[key])

        sig[key] = np.array(sig[key]).transpose()

    for key in ['nu', 'gam', 'lam', 'v']:

        sig[key]=[]

        for nuclide in nuc:

            sig[key].append(nuclide.xs[key])
            #sig[key].append([nuclide.xs[key]]*steps)

        sig[key] = np.array(sig[key], dtype=object).transpose()

    sig['removal'] = np.array(sig['101'], dtype=object) + np.array(sig['18'], dtype=object)

    return sig

def buildWhere(nuc):

        where = np.zeros((len(REG),len(nuc)))

        for a in nuc:

            where[a.reg][a.id]=1

        return where

def versorXS(At,sig):

        xs = sig.copy()

        for key in xs.keys():

            xs[key] = np.zeros((ene,steps,len(At)))

        return xs


def main():

    scatt = buildScatt()

    ### MICRO DEPLETION ###

    Albe, Flux = cycleAlbe(time)
    mdep       = readMDEP()
    zai, ZAI   = getZai(mdep)

    xs   = buildXS(time, tempo, ZAI, mdep)
    Albe = rescaleTime(Albe, time, tempo)

    ### NUCLIDES ###

    nuc, ZAI, MAT = buildNuc(zai, xs)
    sig = buildSig(nuc)

    ### VECTORS ###

    At = [a.at for a in nuc]
    where = buildWhere(nuc)
    xs = versorXS(At,sig)

    print('Done!\n')

    return scatt, Albe, Flux, nuc, sig, zai, ZAI, MAT, where, At, xs

scatt, Albe, Flux, nuc, sig, zai, ZAI, MAT, where, At, xs = main()



