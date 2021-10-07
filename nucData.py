import numpy as np
import json
import math
import serpentTools as ST
from serpentTools.settings import rc
from periodictable import elements
import matplotlib.pyplot as plt

file     = '2_groups'            # INPUT ENERGY GROUPS
fpSwitch =  False                # INPUT FULL NUCLIDE CHART

res = ST.read(file + '/REP_res.m')
dep = ST.read(file + '/REP_dep.m')
time =  dep.days
giorni = dep.days[-1]
T = dep.days[-1]*24*3600
MXT=ST.MicroXSTuple

#tempo = np.linspace(0,5,10).tolist() + np.linspace(5,100,15).tolist()[1:] + np.linspace(100,giorni,25).tolist()[1:]
tempo = np.linspace(0,5,10).tolist() + np.linspace(5,100,20).tolist()[1:] + np.linspace(100,giorni,100).tolist()[1:]
#tempo = np.linspace(0, giorni, 1000)

steps = len(tempo)
nodo = int(len(dep.days)/2)

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

volZr = 2.16817E-01
volHe = 2.84600E-02

volClad = volZr + volHe
volFuel = 6.92643E-01
volCool = 1.14346E+00

VOL=[volFuel, volClad, volCool]

######################
### SERPENT MANUAL ###
######################

chi  = res.universes['20',0.0, 0, 0.0].infExp['infChit']
grid = res.universes['20',0.0, 0, 0.0].groups[::-1]
ene  = len(chi)

### scattering ###

scaFuel = res.universes['20',0.0, 0, 0.0].infExp['infS0'].reshape(ene,ene)*VOL[0]
scaClad = res.universes['2' ,0.0, 0, 0.0].infExp['infS0'].reshape(ene,ene)*VOL[1]
scaCool = res.universes['50',0.0, 0, 0.0].infExp['infS0'].reshape(ene,ene)*VOL[2]

scatt = [scaFuel, scaClad, scaCool]

#scaFuel = [7.40651E-04, 2.06149E-03]
#scaClad = [9.56363E-04, 4.30693E-04]
#scaCool = [2.87179E-02, 2.69137E-03]

#scatt = np.array([scaFuel, scaClad, scaCool]).transpose()

### nubar ###

nuU5 = np.array([2.65]*int((ene/2)) + [2.43]*int((ene/2)))
nuPu = np.array([3.16]*int((ene/2)) + [2.86]*int((ene/2)))

vU5 = 200.7E+6 * 1.6E-19      # J/fiss
vPu = 207.0E+6 * 1.6E-19

#####################
### REGIONS CLASS ###
#####################

reg = []


class Region:

    def __init__(self, name, id, vol):
        self.id = id
        self.vol = vol
        self.name = name


def regDef(name, id, vol):
    while len(reg) < (id + 1):
        reg.append(0)
    cla = Region(name, id, vol)
    reg[cla.id] = cla
    return cla

fuel = regDef('Fuel', 0, volFuel)
clad = regDef('Cladding', 1, volClad)
cool = regDef('Coolant', 2, volCool)


##################
### XS OBJECTS ###
##################

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

####################
### SERPENT AUTO ###
####################



#corrBoro = np.concatenate((np.ones(int(steps/50*5+1))*2.1,np.ones(steps-int(steps/50*5))*1.6) , axis=0)

###############
### ALBEDOS ###
###############

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

    reg=['fuel', 'clad', 'cool']

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



Albe, Flux = cycleAlbe(time)
keff = ST.read(file+'/REP_res.m').resdata['absKeff'][:,0]


#######################
### MICRO DEPLETION ###
#######################

def readMDEP():

    print('\nParsing nuclear data')

    mdep = []

    for i in range(len(time)):

        mdep.append(ST.read(file+'/REP_mdx'+str(i)+'.m'))

    return mdep

mdep=readMDEP()

UNI = ['20', '2', '50']
REACTIONS = ['101', '102', '103', '107', '16', '17', '18' ]
#sama = [621481, 621501, 621511, 631491, 641521, 631511, 611491, 631491, 621491, 621501, 631501]
sama =[]
def getzaiFuel():

    zaiFuel=[]

    for key in mdep[0].xsVal['20'].keys():

        if key[0] not in [60000, 300000] and key[1] == 101:
        #if key[0] not in [60000, 300000] and key[1] == 101:

            zaiFuel.append(key[0])

    zaiFuel.extend(sama)

    zai=sorted(zaiFuel)

    out=[str(z) for z in zai]

    return out
def getzaiElse():
    zaiElse = []

    for key in mdep[0].xsVal['2'].keys():

        #if key[0] not in [130270, 10010, 10020, 80160] and key[1] == 101:
        if key[1] == 101:

            if mdep[0].getXS(universe='2', isotope=key[0], reaction=101)[0][1] != 0:

                zaiElse.append(key[0])

    zai=sorted(zaiElse)

    out=[str(z) for z in zai]

    return out

zaiOld=['531350','541350', '601490', '611490', '621490', '922340', '922350', '922380', '922390', '932390', '942390']
zaiFuel=zaiOld
zaiElse=['400900', '400910', '400920', '400940', '400960', '80160']

if fpSwitch == True:

    zaiFuel=getzaiFuel()
    zaiElse=getzaiElse()

ZAI=zaiFuel+zaiElse+['20040']+['80160']+['10010']

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
def buildXS(xsTime, newTime):

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



xs = buildXS(time, tempo)
Albe = rescaleTime(Albe, time, tempo)


################
### NUCLIDES ###
################

######################
### NUCLIDES CLASS ###
######################

nuc=[]
ZAI=[]
MAT=[]

class Nuclide:

    def __init__(self, name, idReg, zai, mat, vol, **kwargs):

        if zai in xs[UNI[idReg]].keys() :

            for key in xs[UNI[idReg]][zai].keys():

                xs[UNI[idReg]][zai][key] = (np.array(xs[UNI[idReg]][zai][key]) / vol * reg[idReg].vol).tolist()

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

        if int(zai[-1]) != 0:
            self.dens = 0
            self.comp = 0

        else:
            self.dens = dep.materials[mat].getValues('days','mdens', zai=int(zai))[0][0]
            self.comp = (dep.materials[mat].getValues('days','adens', zai=int(zai))[0][0]+1E-40) * vol / reg[idReg].vol * 6.022E+23

        if 'xmat' in kwargs.keys():
            self.dens+=dep.materials[kwargs['xmat']].getValues('days', 'mdens', zai=int(zai))[0][0]
            self.xmat = kwargs['xmat']

        else:
            self.xmat = False

        self.vol = vol
        self.at = self.dens / getMM(zai) * vol * 6.022E+23


def nucDef(name, idReg, zai, mat, vol, **kwargs):

    cla = Nuclide(name, idReg, zai, mat, vol, **kwargs)

    nuc.append(cla)
    ZAI.append(cla.zai)
    MAT.append(cla.mat)

    return cla


def buildFuel():

    Fuel=[]
    u='20'

    for z in zaiFuel:

        nu = nuU5
        v  = vU5

        if z == '942390':
            nu = nuPu
            v  = vPu

        elem=nucDef(z, 0, z, 'fuel', volFuel, nu=nu, v=v)
        Fuel.append(elem)

    print('Done!\n')

    return Fuel
def buildClad():

    Clad=[]
    u='2'

    for z in zaiElse:

        elem=nucDef(z, 1, z, 'clad', volClad)
        Clad.append(elem)

    print('Done!\n')

    return Clad

FUEL=buildFuel()
CLAD=buildClad()


### absorbers ###

#B = nucDef('Boron-10 (belt)', 2, '50100', 'boro', volB )
He = nucDef('Helium', 1, '20040', 'helium', volHe )

### heavy water ###

OxC = nucDef('Oxygen-coolant',  2, '80160', 'water', volCool )
H2C = nucDef('Hydrogen-coolant',  2, '10010', 'water', volCool )

#####################
### XS DICTIONARY ###
#####################


#Pm.xs['101']=[[0,0]]*len(Pm.xs['101'])
#print(Pm.xs['101'])

def buildSig():

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

sig = buildSig()

####################
### ATOMS VECTOR ###
####################

At = [a.at for a in nuc]

where = np.zeros((len(reg),len(nuc)))

for a in nuc:

    where[a.reg][a.id]=1

xs = sig.copy()

for key in xs.keys():

    xs[key] = np.zeros((ene,steps,len(At)))


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

    fig.savefig('xs/'+zai+'_'+MT+'_xs.png')




