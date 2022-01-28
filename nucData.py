import config
import serpent
import numpy as np
import serpentTools as ST
from periodictable import elements
import matplotlib.pyplot as plt
from datetime import datetime
import copy


startNuc = datetime.now()


model     = config.model                                    # INPUT MODEL
energy    = config.energy                                      # INPUT ENERGY GROUPS
dayStop   = config.dayStop
PASSI     = config.PASSI                                       # INPUT STEP NUMBER
fpSwitch  = config.fpSwitch                                        # SWITCH TO FULL NUCLIDE CHART
hetSwitch = config.hetSwitch                                        # SWITCH TO HETEROGENEOUS CORRECTION FOR FUEL AND NICHEL

### INITS ###

input = serpent.inp[model]
file = model+'/'+str(energy)+'_groups'
res = ST.read(file + '/'+input+'_res.m')
dep = ST.read(file + '/'+input+'_dep.m')
keff = ST.read(file+'/'+input+'_res.m').resdata['absKeff'][:,0]
time =  dep.days
giorni = dep.days[-1]*dayStop
T = dep.days[-1]*24*3600
MXT=ST.MicroXSTuple

### TIME STEPS ###

def buildHet(PASSI):

    pass1 = 5
    pass2 = round(giorni*0.1)
    pass3 = giorni - pass2
    pass4 = giorni - pass1

    het1 = round(PASSI*0.1)
    het2 = round(PASSI*0.15)
    het3 = PASSI -2*(het1+het2)

    tempo_het  = np.linspace(0,pass1,het1).tolist() + np.linspace(pass1,pass2,het2).tolist()[1:] + np.linspace(pass2,pass3,het3).tolist()[1:] + np.linspace(pass3, pass4, het2).tolist()[1:] + np.linspace(pass4, giorni, het1).tolist()[1:]

    return tempo_het

tempo_het = buildHet(PASSI)
tempo_homo = np.linspace(0, giorni, PASSI)
tempo = tempo_het

if config.hetSteps == False:

    tempo = tempo_homo

steps = len(tempo)
#nodo = min([int(len(dep.days)/2),int(len(tempo)/2)])
nodo = 5


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

reaz = ['fission', '(n,gamma)', '(n,2n)', '(n,3n)', '(n,p)', '(n,a)','nubar', 'removal', 'sca']
MT = ['18', '102', '16', '17', '103', '107', '452', 'removal', 'sca']

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

def getName(zai, idReg):

    z=zai[:-4]
    a=zai[-4:-1]

    if int(a) in elements[int(z)].isotopes:

        Name=elements[int(z)][int(a)].name

    elif zai == 60000:

        Name = 'carbon'

    else:

        Name = 'sarcazzo'

    Name += '-'+str(a)
    Name=Name[0].upper()+Name[1:]+ ' (' + REG[idReg] + ')'

    return Name

def xsPlot(sigma, MT, zai, t):

    x = grid
    id= ZAI.index(zai)
    name = nuc[id].name
    fig, axs = plt.subplots()

    y = [0] + [sigma[MT][e][t][id]*1E+24/(x[e]**0) for e in range(ene)]
    axs.step(x, y, where = 'pre')
    tit = name+ ' MT '+MT+' cross section\n'
    axs.set(xlabel='Energy [MeV]', ylabel='cross section [barn/MeV]', title=tit)
    axs.set_xlim(1E-12,1E+7)
    #axs.set_ylim(min(y)/10,max(y)*10)
    axs.set_yscale('log')
    axs.set_xscale('log')
    print(y)
    fig.savefig(model+'/xs/'+zai+'_'+MT+'_xs.png')

def rescaleTime(xs, xsTime, newTime):

    newXs = []

    for t in newTime:

        j = 0

        while xsTime[j]+1 < t and j < (len(xsTime)-2):

            j+=1

        width = time[j+1] - time[j]

        newXs.append(xs[j]*(xsTime[j+1]-t)/width+xs[j+1]*(t-xsTime[j])/width)
        #newXs.append(xs[nodo])

    return newXs

def readMDEP():

    print('\nScanning nuclear data')

    mdep = []

    for i in range(len(time)):

        mdep.append(ST.read(file+'/'+input+'_mdx'+str(i)+'.m'))

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

    def __init__(self, zlist, idReg, z, mat, vol, xs, id, CB, CF, CN, **kwargs):

        if z[-1] in ['0','1']:

            self.name = getName(z, idReg)

        else:

            self.name = z

        if z[-1] == '1':

            self.name += ' [metastable]'

        self.id  = id
        self.reg = idReg
        self.zai = z

        self.vol  = 1E-4
        self.at   = 0
        self.comp   = 0
        self.mat  = []

        d = []

        for i in range(len(mat)):

            if z not in ['621481', '50100']+[str(a) for a in serpent.sama]:
                if dep.materials[mat[i]].getValues('days','mdens', zai=int(z))[0][nodo] > 0 :
                    self.mat.append(mat[i])
                    self.at   += dep.materials[mat[i]].getValues('days','mdens', zai=int(z))[0][0] / getMM(z) * vol[i] * 6.022E+23
                    self.vol  += vol[i]
                d.append(dep.materials[mat[i]].getValues('days','mdens', zai=int(z))[0][0])

            else:
                self.mat = serpent.materials[model][fuelId][0]
                self.vol = serpent.volumes[model][fuelId][0]
                d.append(0)

        if self.mat == [] and idReg == fuelId:
            self.mat = serpent.materials[model][fuelId][0]

        #self.vol = vol[np.argmax(np.array(d))]

        if z in xs[UNI[idReg]].keys() and z != '50100':

            for key in xs[UNI[idReg]][z].keys():

                xs[UNI[idReg]][z][key] = (np.array(xs[UNI[idReg]][z][key]) / self.vol * VOL[idReg]).tolist()
            self.xs = {**xs[UNI[idReg]][z], **xsDef(**kwargs)}

        else:

            self.xs = xsDef(**kwargs)



        if z == '50100':

            self.mat = ['boro']

            XS = xs.copy()

            if idReg == 1:


                if model[:3] == 'LEU':

                    for key in xs[UNI[idReg]][z].keys():
                        xs[UNI[idReg]][z][key] = (np.array(XS[UNI[2]][z][key]) * 1 * np.array(CB)).tolist()

                if model[:3] == 'HEU':

                    for key in xs[UNI[idReg]][z].keys():
                        xs[UNI[idReg]][z][key] = (np.array(XS[UNI[1]][z][key]) * 1 * np.array(CB)).tolist()

                self.xs = {**xs[UNI[idReg]][z], **xsDef(**kwargs)}

                self.at = 0


            if idReg == 2:

                if model[:3] == 'LEU':

                    self.vol = serpent.volumes[model][2][serpent.materials[model][2].index('boro')]

                    for key in xs[UNI[idReg]][z].keys():
                       xs[UNI[idReg]][z][key] = (np.array(XS[UNI[2]][z][key])  * 1 * np.array(CB)).tolist()

                if model[:3] == 'HEU':

                    self.vol = serpent.volumes[model][1][serpent.materials[model][1].index('boro')]

                    for key in xs[UNI[2]][z].keys():
                        xs[UNI[idReg]][z][key] = (np.array(XS[UNI[1]][z][key]) * 1 * np.array(CB)).tolist()

                self.xs = {**xs[UNI[idReg]][z], **xsDef(**kwargs)}

                self.at = dep.materials['boro'].getValues('days', 'mdens', zai=int(z))[0][0] / getMM(z) * self.vol * 6.022E+23

            else:

                for key in xs[UNI[idReg]][z].keys():
                    xs[UNI[idReg]][z][key] = (np.array(xs[UNI[idReg]][z][key]) / self.vol * VOL[idReg]).tolist()
                self.xs = {**xs[UNI[idReg]][z], **xsDef(**kwargs)}


        if 'Fuel' in mat and int(z) > 300000 and z not in ['621481']+[str(a) for a in serpent.sama] and hetSwitch == True :

            self.mat = ['Fuel']
            self.vol = vol[mat.index('Fuel')]
            self.at  = dep.materials['Fuel'].getValues('days', 'mdens', zai=int(z))[0][0] / getMM(z) * self.vol * 6.022E+23

            for key in xs[UNI[idReg]][z].keys():
               xs[UNI[idReg]][z][key] = (np.array(xs[UNI[idReg]][z][key]) / self.vol * VOL[idReg] * np.array(CF)).tolist()

            self.xs = {**xs[UNI[idReg]][z], **xsDef(**kwargs)}

        if 'Ni' in mat and z in ['280580', '280600'] and hetSwitch == True:

            self.mat = ['Ni']
            self.vol = vol[mat.index('Ni')]
            self.at  = dep.materials['Ni'].getValues('days', 'mdens', zai=int(z))[0][0] / getMM(z) * self.vol * 6.022E+23

            for key in xs[UNI[idReg]][z].keys():
               xs[UNI[idReg]][z][key] = (np.array(xs[UNI[idReg]][z][key]) / self.vol * VOL[idReg]  * np.array(CN)).tolist()

            self.xs = {**xs[UNI[idReg]][z], **xsDef(**kwargs)}

        #print(self.vol)

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

    if model[1:3] == 'EU':

        R_out=(abs(det.detectors['bordo'].tallies[::-1])+abs(det.detectors['bordoup'].tallies[::-1])+abs(det.detectors['bordodown'].tallies[::-1]))

    RR=np.array([[R[0,1], -R[0,0], zero],[-R[0,1], (R[0,0]+R[1,1]), -R[1,0]],[zero, -R[1,1], R[1,0]+R_out]])

    F=[]

    reg=serpent.detectors[model]

    for i in range(len(reg)):

        #F.append((det.detectors[reg[i]].tallies[::-1]+np.ones(ene)*1E-5)/VOL[i])
        F.append((det.detectors[reg[i]].tallies[::-1])/VOL[i])

    S=(RR/F).transpose()

    SS=np.array([S[j].transpose() for j in range(ene)])

    F=np.array(F).transpose().ravel()

    CB=[]
    CN=[]
    CF=[]

    if model[:3] == 'LEU':

        CB=(np.array(det.detectors['boro'].tallies[::-1]/vol[2][mat[2].index('boro')])/np.array(det.detectors[DET[2]].tallies[::-1]/VOL[2])*1).tolist()
        CN=(np.array(det.detectors['Ni'].tallies[::-1]/vol[0][mat[0].index('Ni')])/np.array(det.detectors[DET[0]].tallies[::-1]/VOL[0])*1).tolist()
        CF=(np.array(det.detectors['meat'].tallies[::-1]/vol[1][mat[1].index('Fuel')])/np.array(det.detectors[DET[1]].tallies[::-1]/VOL[1])*1).tolist()

    if model[:3] == 'HEU':
        CB = (np.array(det.detectors['boro'].tallies[::-1] / vol[1][mat[1].index('boro')]) / np.array(det.detectors[DET[1]].tallies[::-1] / VOL[1]) * 0.65).tolist()
        CN = (np.array(det.detectors['Ni'].tallies[::-1] /vol[0][mat[0].index('Ni')]) / np.array(det.detectors[DET[0]].tallies[::-1] / VOL[0]) * 1).tolist()
        CF = (np.array(det.detectors['meat'].tallies[::-1] / vol[1][mat[1].index('Fuel')]) / np.array(det.detectors[DET[1]].tallies[::-1] / VOL[1]) * 1).tolist()

    return SS, F, CB, CN, CF

def cycleAlbe(time):

    Albe=[]
    Flux=[]
    corrBoro=[]
    corrNi=[]
    corrFuel=[]

    for t in range(len(time)):

        det = ST.read(file+'/'+input+'_det'+str(int(t))+'.m')
        Albedo, Flusso, CB, CN, CF = buildAlbe(det)
        Albe.append(Albedo)
        Flux.append(Flusso)
        corrBoro.append(CB)
        corrNi.append(CN)
        corrFuel.append(CF)

    return Albe, Flux, np.array(corrBoro), np.array(corrNi), np.array(corrFuel)

def buildXS(xsTime, newTime, ZAI, mdep):

    xs={}

    for u in UNI:

        xs[u]={}

        for z in ZAI:

            xs[u][z]={}

            for r in REACTIONS:

                if MXT(zai=int(z),mt=int(r),metastable=0) in mdep[nodo].xsVal[u].keys() or MXT(zai=int(z),mt=int(r),metastable=1) in mdep[1].xsVal[u].keys() :

                    xs[u][z][r]=[]

                    for i in range(len(mdep)):

                        xs[u][z][r].append(mdep[i].getXS(universe=u, isotope=int(z), reaction=int(r))[0]*1E-24)

                    xs[u][z][r]=rescaleTime(xs[u][z][r], xsTime, newTime)

                else:

                    xs[u][z][r]=np.zeros((len(newTime),ene)).tolist()
                    #xs[u][z][r]=np.zeros((len(xsTime),2)).tolist()

        #xs[u]['10020']['101'] = ( np.array(xs[u]['80160']['101'])/2 + np.array(xs[u]['10020']['101'])).tolist()


    return xs

def buildNuc(zai, xs, CB, CF, CN):

    nuc = []
    MAT = []
    NOM = []

    nuU5   = np.array([[2.65] * int((ene / 2)) + [2.43] * int((ene / 2))]*steps).tolist()
    nuPu   = np.array([[3.16] * int((ene / 2)) + [2.86] * int((ene / 2))]*steps).tolist()
    nuZero = np.array([[0] * ene]*steps).tolist()

    vU5 = 200.7E+6 * 1.6E-19  # J/fiss
    vPu = 207.0E+6 * 1.6E-19

    j = 0

    for i in range(len(UNI)):

        if i == fuelId:

            for z in zai[i]:

                nu = nuU5
                v  = vU5

                if z == '942390':
                    nu = nuPu
                    v  = vPu

                elem = Nuclide(zai, i, z, mat[i], vol[i], xs, j, CB, CF, CN, nu=nu, v=v)
                nuc.append(elem)
                MAT.append(elem.mat)
                NOM.append(elem.name)

                j+=1
        else:

            for z in zai[i]:
                elem = Nuclide(zai, i, z, mat[i], vol[i], xs, j, CB, CF, CN, nu=nuZero)
                nuc.append(elem)
                MAT.append(elem.mat)
                NOM.append(elem.name)
                j+=1


    return nuc, MAT, NOM

def buildSig(nuc):

    sig = {}

    for key in REACTIONS+['sca', 'nu']:

        sig[key]=[]

        for nuclide in nuc:

            sig[key].append(nuclide.xs[key])

        sig[key] = np.array(sig[key]).transpose()

    for key in ['gam', 'lam', 'v']:

        sig[key]=[]

        for nuclide in nuc:

            sig[key].append(nuclide.xs[key])
            #sig[key].append([nuclide.xs[key]]*steps)

        sig[key] = np.array(sig[key], dtype=object).transpose()

    sig['removal'] = np.array(sig['101'], dtype=object) + np.array(sig['18'], dtype=object)
    sig['452']     = copy.deepcopy(sig['nu'])

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

    Albe, Flux, corrBoro, corrNi, corrFuel = cycleAlbe(time)
    mdep       = readMDEP()
    zai, ZAI   = getZai(mdep)

    xs   = buildXS(time, tempo, ZAI, mdep)
    Albe = rescaleTime(Albe, time, tempo)
    corrBoro = np.array(rescaleTime(corrBoro, time, tempo))
    corrFuel = np.array(rescaleTime(corrFuel, time, tempo))
    corrNi   = np.array(rescaleTime(corrNi, time, tempo))

    ### NUCLIDES ###

    nuc, MAT, NOM = buildNuc(zai, xs, corrBoro, corrFuel, corrNi)
    sig = buildSig(nuc)

    ### VECTORS ###

    At = [a.at for a in nuc]
    comp = [a.comp for a in nuc]
    where = buildWhere(nuc)
    XS = versorXS(At,sig)

    print('Done!\n')

    return scatt, Albe, Flux, nuc, sig, zai, ZAI, MAT, NOM, where, At, comp, XS

scatt, Albe, Flux, nuc, sig, zai, ZAI, MAT, NOM, where, At, comp, xs = main()

endNuc = datetime.now()


#xsPlot(sig,'18','922350',0)
#xsPlot(sig,'102','922380',0)

#for a in nuc:

#    print([str(a.vol)] + [a.zai])

nuc[ZAI.index('922350')].name = 'Uranium-235'
nuc[ZAI.index('922380')].name = 'Uranium-238'
nuc[ZAI.index('942390')].name = 'Plutonium-239'
