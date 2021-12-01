import numpy as np
import scipy
import sympy as sp
from sympy import Symbol
import nucData
import onix
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from periodictable import elements
import copy
import serpent

### NUCLEAR DATA ###

sig   = nucData.sig
steps = nucData.steps
chi   = nucData.chi
where = nucData.where
ene = nucData.ene
reg = 3
lenPsi = ene*reg

### Results object ###

class Results:

    def __init__(self):
        self.comp = []
        self.flux = []
        self.keff = []
        self.M = []
        self.A = []
        self.source = []
        self.pow = []
        self.phi = []
        self.ind =[]
        self.pert={}
        self.homo=[]

### SIBYL ###

def buildFiss(N, sigma, t, **kwargs):

    reg = 3

    F = np.zeros((ene, ene, reg, reg))

    for i in range(ene):
        for j in range(ene):

            F[i][j] = np.diag(N.dot(sigma['18'][i][t] * sig['nu'][i])) * chi[j]

    mtx = np.vstack(np.dstack(np.array(F)))

    if 'ribalta' in kwargs.keys() :
        if kwargs['ribalta'] == True:
            mtx = ribalta(np.array(F))

    return mtx

def buildPozz(N, sigma, t, **kwargs):

    reg = 3

    A = np.zeros((ene, ene, reg, reg))

    for i in range(ene):

        A[i][i] = np.diag(N.dot(sigma['removal'][i][t]))

    mtx = np.vstack(np.dstack(np.array(A)))

    if 'ribalta' in kwargs.keys() :
        if kwargs['ribalta'] == True:
            mtx = ribalta(np.array(A))

    return mtx

def buildTrasp(Albe, t, **kwargs):

    reg = 3

    mtx = np.zeros((ene, ene, reg, reg))

    for i in range(ene):

        mtx [i][i] = Albe[t][i]

    Trasp = np.vstack(np.dstack(np.array(mtx)))

    if 'ribalta' in kwargs.keys() :
        if kwargs['ribalta'] == True:
            Trasp = ribalta(np.array(mtx))

    return Trasp

def buildSca(scatt, **kwargs):

    ene = len(scatt[0][0])
    reg = 3

    lista = []

    for i in range(ene):
        for j in range(ene):

            lista.append(np.diag([-scatt[0][i][j], -scatt[1][i][j], -scatt[2][i][j]]))

            if i == j:

                pozzo1 = sum(scatt[0][:][j])
                pozzo2 = sum(scatt[1][:][j])
                pozzo3 = sum(scatt[2][:][j])

                lista[-1] += np.diag([pozzo1, pozzo2, pozzo3])

    A = np.array(lista).reshape(ene,ene,reg,reg)

    sca = np.vstack(np.dstack(A))

    if 'ribalta' in kwargs.keys() :
        if kwargs['ribalta'] == True:
            sca = ribalta(A)

    return sca

def boltzL(N, sigma, t, **kwargs):

    S = buildSca(nucData.scatt, **kwargs)

    T = buildTrasp(nucData.Albe, t, **kwargs)

    A = buildPozz(N, sigma, t, **kwargs)

    if N.ravel()[-1] == 0:
        S=S*0
        T=T*0

    return (A+S+T)

def boltzF(N, sigma, t, **kwargs):

    F = buildFiss(N, sigma, t, **kwargs)

    return F

def Boltz(N, sigma, t, lam, **kwargs):
    A = boltzL(N, sigma, t, **kwargs) - lam * boltzF(N, sigma, t, **kwargs)

    return -A

def ribalta(LL):

    forma1 = (ene, ene, reg, reg)
    forma2 = (reg, reg, ene, ene)
    L = LL.copy()
    P = np.zeros(forma2)

    for a in range(forma1[0]):
        for b in range(forma1[1]):
            for c in range(forma1[2]):
                for d in range(forma1[3]):

                    P[c,d,a,b]=L[a,b,c,d]

    mtx = np.vstack(np.dstack(np.array(P)))

    return mtx

### SALAMECHE ###

def reshapePsi(Psi):

    return np.array(Psi.reshape(ene,reg).transpose())

reaz = ['fission', '(n,gamma)', '(n,2n)', '(n,3n)', '(n,p)', '(n,a)', 'removal', 'sca']
MT = ['18', '102', '16', '17', '103', '107', 'removal', 'sca']

def rr(sigma, psi, Phi, t):
    rr = copy.deepcopy(nucData.xs)

    rr['fission'] = 0
    rr['(n,gamma)'] = 0
    rr['(n,2n)'] = 0
    rr['(n,3n)'] = 0
    rr['(n,p)'] = 0
    rr['(n,a)'] = 0
    rr['removal'] = 0
    rr['sca'] = 0


    for i in range(3):

        for j in range(ene):

            for jj in range(len(reaz)-1):

                rr[reaz[jj]] += np.inner(np.array(sigma[MT[jj]][j][t]) * where[i], psi[i][j]) * Phi

    return rr

def buildPL(zai):

    print('\nBuilding BU Matrix')

    PL = onix.Passlist(zai)

    for i in range(len(PL.nucl_list)):

        if PL.nucl_list[i] in onix.data.default_xs_lib.keys():

            PL.passport_list[i].load_xs()

            for j in PL.passport_list[i].current_xs.keys():
                PL.passport_list[i].current_xs[j] = [0.0, 0.0]

        elif PL.nucl_list[i] not in onix.data.default_xs_lib.keys():

            PL.passport_list[i].current_xs = {}


        if PL.nucl_list[i] in onix.data.default_decay_lib_a.keys():

            PL.passport_list[i].load_decay()

        if PL.passport_list[i].get_FAM() != 'ACT' and PL.nucl_list[i] in onix.data.default_fy_lib.keys():

            PL.passport_list[i].load_fy()


    #print(PL.passport_list[5].zamid)
    #print(PL.passport_list[5].current_xs)
    #print(rr['fission'][5])
    #print(PL.passport_list[nucData.zai[nucData.fuelId].index('621490')].all_parent)
    #print(PL.passport_list[nucData.zai[nucData.fuelId].index('621490')].creation_dic)

    print('Done!\n')

    return PL

zai = nucData.zai[nucData.fuelId]

p = 0

if nucData.fuelId != 0:

    p = len(nucData.zai[nucData.fuelId-1])

    if nucData.fpSwitch == False:

        zai = zai[:-4]

pl=buildPL(zai)

def updatePL(PL, RR):

    for j in range(len(PL.nucl_list)):

        i = j+p

        for key in reaz[:-1]:

            PL.passport_list[j].current_xs[key]=[RR[key][i], 0.0]

    return PL

def onixR(PL):

    R=onix.salameche.get_xs_mat(PL)

    zeros = np.diag(np.zeros(len(nucData.nuc)))

    zeros[p:p+R.shape[0], p:p+R.shape[1]] = R

    return zeros

def onixD(PL):

    D = onix.salameche.get_decay_mat(PL)

    zeros = np.diag(np.zeros(len(nucData.nuc)))

    zeros[p:p+D.shape[0], p:p+D.shape[1]] = D

    return zeros

def Bateman(rr):

    PL=updatePL(pl, rr)

    M = np.array(onixR(PL)) + np.array(onixD(PL))

    burn=['280580', '280600', '50100']

    for i in range(len(nucData.ZAI)):

        if nucData.ZAI[i] in burn:

            M[i,i]+=-rr['removal'][i]

    return tuple(map(tuple, M))


### PTERODAx ###

def tramezzino(Ns,N,R):

    return float(np.inner(np.array(Ns), np.array(R).dot(np.array(N))))

def I(Ns, N, R, dt):

    a=[(N[1][i]-N[0][i])/dt for i in range(len(N[0]))]
    b=[(Ns[1][i]-Ns[0][i])/dt for i in range(len(Ns[0]))]

    NN=(np.array(N[1])+np.array(N[0]))/2
    NS=(np.array(Ns[1])+np.array(Ns[0]))/2

    #I = tramezzino(N[0],Ns[0],R)*t + tramezzino(N[0],b,R)*t**2/2 + tramezzino(a,Ns[0],R)*t**2/2 + tramezzino(a,b,R)*t**3/3
    I = tramezzino(NS,NN,R)*dt

    return I

def kSens(Gh, Psi, N, k, pertId, t):

    dN = np.zeros(len(N))
    dN[pertId] = 1

    num1 = tramezzino(Gh, Psi, np.matrix(Boltz(dN*where, sig, t, 1/k)))
    num = tramezzino(Gh, Psi, np.matrix(Boltz(N*where, sig, t, 1/k)))

    den1 = tramezzino(Gh, Psi, boltzF(dN*where, sig, t))/k
    den = tramezzino(Gh, Psi, boltzF(N*where, sig, t))/k

    sens = (num1*den-num*den1)/(den**2)*k
    #sens = -num1/den*k

    return sens

def Qs(N, Ns, Ps, Phi, t, dt):
    sf = np.array(sig['18'])

    Q = []

    for i in range(lenPsi):

        Psi = np.zeros(lenPsi)

        Psi[i] = 1

        if i%3 == nucData.fuelId:

            fis = Ps * np.inner(sf[int(i/3)][t]*sig['v'], np.array(N[0]))

        else:

            fis = 0

        r=rr(sig, reshapePsi(Psi), Phi, t)
        R = Bateman(r)-onixD(updatePL(pl,r))
        M = Bateman(r)
        Q.append(I(Ns,N,R,dt) - fis*Phi)

    return Q

def beta(Psi, lam, N, t):
    dN = np.zeros(len(N))

    BETA = []

    for i in range(len(N)):

        dN[i]  = 1
        dN[-1] = 0


        DN = dN * where

        BETA.append(Boltz(DN, sig, t, lam).dot(np.array(Psi)))

        dN = np.zeros(len(N))

    return np.array(BETA)

def beta2(Gh, lam, N, t):
    dN = np.zeros(len(N))

    BETA = []

    for i in range(len(N)):

        dN[i]  = 1
        dN[-1] = 0

        DN = dN * where

        BETA.append(Boltz(DN, sig, t, lam).transpose().dot(np.array(Gh)))

        dN = np.zeros(len(N))

    return np.array(BETA)

def pi(Psi, Phi):
    sf = np.array(sig['18'])

    return (sf).dot(np.array(Psi)) * Phi

def dR(Psi, Gh, N, k, t):

    dR = []

    for i in range(len(Psi)):
        dPsi = np.zeros(len(Psi))
        dPsi[i] = 1

        num1 = tramezzino(Gh, dPsi, np.matrix(Boltz(N * where, sig, t, 1/k )))
        num = tramezzino(Gh, Psi, np.matrix(Boltz(N * where, sig, t, 1/k )))

        den1 = tramezzino(Gh, dPsi, boltzF(N * where, sig, t))/k
        den = tramezzino(Gh, Psi, boltzF(N * where, sig, t))/k

        #dR.append(num1/den)
        dR.append((num1*den-num*den1)/(den**2))

    return np.array(dR)*k

def dR2(Psi, Gh, N, k, t):

    dR = []

    for i in range(len(Psi)):
        dGh = np.zeros(len(Psi))
        dGh[i] = 1

        num1 = tramezzino(dGh, Psi, np.matrix(Boltz(N * where, sig, t, 1/k )))
        num = tramezzino(Gh, Psi, np.matrix(Boltz(N * where, sig, t, 1/k )))

        den1 = tramezzino(dGh, Psi, boltzF(N * where, sig, t))/k
        den = tramezzino(Gh, Psi, boltzF(N * where, sig, t))/k

        #dR.append(num1/den)
        dR.append((num1*den-num*den1)/(den**2))

    return np.array(dR)*k

### ND ###

def betaSig(Psi, lam, pert, N, G, id, t):

    M = N.copy()
    BETA = []
    M[-1] = 0

    for e in range(ene):

        XS = copy.deepcopy(nucData.xs)
        XS[pert][e][t][id]=1
        XS['removal'][e][t][id]=1

        B = Boltz(M * where, XS, t, lam)

        BETA.append(tramezzino(G,Psi,B))

    return np.array(BETA)

def beta2Sig(Gh, lam, pert, N, G2, id, t):

    M = N.copy()
    BETA = []
    M[-1] = 0

    for e in range(ene):

        XS = copy.deepcopy(nucData.xs)
        XS[pert][e][t][id]=1
        XS['removal'][e][t][id]=1

        B = Boltz(M * where, XS, t, lam).transpose()

        BETA.append(tramezzino(G2,Gh,B))

    return np.array(BETA)

def bateSig(Psi, Phi, pert, N, Ns, id, t, dt):

    BATE = []
    PSI = reshapePsi(Psi)

    for e in range(ene):

        XS = copy.deepcopy(nucData.xs)
        XS[pert][e][t][id]=1
        XS['removal'][e][t][id]=1

        RR = rr(XS,PSI,Phi,t)
        PL = updatePL(pl,RR)

        R = onixR(PL)

        BATE.append(I(Ns,N,R, dt))

    return np.array(BATE)

def PiSig(Psi, Phi, pert, N, id, t):

    PI = []
    PSI = reshapePsi(Psi)

    for e in range(ene):

        XS = copy.deepcopy(nucData.xs)
        XS[pert][e][t][id]=1
        XS['removal'][e][t][id]=1


        RR = rr(XS,PSI,Phi,t)

        PI.append(RR['fission'][id]*sig['v'][id]*N[id])

    return np.array(PI)

def kSensSig(Gh, Psi, N, k, pertId, pertMT, t):

    M = N.copy()
    sens = []
    M[-1] = 0

    for e in range(ene):

        XS = copy.deepcopy(nucData.xs)
        XS[pertMT][e][t][pertId]=1
        XS['removal'][e][t][pertId]=1

        num1 = tramezzino(Gh, Psi, np.matrix(Boltz(M*where, XS, t, 1/k)))
        num = tramezzino(Gh, Psi, np.matrix(Boltz(N*where, sig, t, 1/k)))

        den1 = tramezzino(Gh, Psi, boltzF(M*where, XS, t))/k
        den = tramezzino(Gh, Psi, boltzF(N*where, sig, t))/k

        sens.append((num1*den-num*den1)/(den**2)*k)
        #sens = -num1/den*k

    return np.array(sens)

### RICERCA AUTOVALORE ###

def findK(N, t):
    lam = 0.5

    start = np.linalg.det(Boltz(N, sig, t, lam))

    maj = start * 1E-2

    for x in [1E-3, 1E-6, 1E-8, 1E-10, 1E-12, 1E-14]:

        while abs(np.linalg.det(Boltz(N, sig, t, lam))) > maj:
            lam = lam + x

        lam = lam - x

        maj = maj * 1E-2

    lam = lam + x

    return 1 / lam

def solveK(N, inter, t):
    x = Symbol('x')

    lam = sp.solveset(sp.Matrix(boltzL(N, t) - x * boltzF(N, t)).det(), domain=inter)

    return float(1 / list(lam)[0])

def crushK(N, sigma, t, a,b):

    def f(x):

        F = scipy.linalg.det(boltzL(N, sigma, t) - x * boltzF(N, sigma, t))

        return F

    lam = scipy.optimize.brentq(f, a, b)

    def plotDet():

        x = np.linspace(0, 2, 1000)
        y = [f(a) for a in x]
        fig, ax1 = plt.subplots()
        ax1.plot(x, y)
        fig.savefig('det/'+ str(t) + '.png')

    #plotDet()

    return 1 / lam

def kInf(At, rr):

    fis = np.inner(rr['18'], At) * 2.43

    cap = np.inner(rr['102'], At)

    return fis / cap

def moveCR(N,t, k):

    p = nucData.ZAI.index('280580')
    q = nucData.ZAI.index('280600')

    a = 0
    b = 3

    MM = N.copy()

    def g(x):

        M = N.copy()

        M[p] = x*N[p]
        M[q] = x*N[q]

        return M

    def f(x):

        F = scipy.linalg.det(Boltz(g(x)*where, sig, t, 1/k))

        return F

    #print(f(a))
    #print(f(b))
    #print(Boltz(g(a)*where, sig, t, 1/k))
    #print(Boltz(g(b)*where, sig, t, 1/k))


    y = scipy.optimize.brentq(f, a, b, rtol=1E-10)

    MM[p]=N[p]*y
    MM[q]=N[q]*y

    def buildSM():

        SM = np.diag(np.ones(len(N)))

        SM[p][p] = y
        SM[q][q] = y

        #al = nucData.AlC.id
        #D2 = nucData.D2C.id
        #H2 = nucData.H2C.id
        #Ox = nucData.OxC.id

        #SM[al][al] = y
        #SM[D2][D2] += (1/y-1)*(nucData.volNi+nucData.volAl)/nucData.volHwC
        #SM[H2][H2] += (1/y-1)*(nucData.volNi+nucData.volAl)/nucData.volHwC
        #SM[Ox][Ox] += (1/y-1)*(nucData.volNi+nucData.volAl)/nucData.volHwC

        return SM

    SM = buildSM()

    def plotDet():

        x = np.linspace(0, 2, 1000)
        y = [f(a) for a in x]
        fig, ax1 = plt.subplots()
        ax1.plot(x, y)
        fig.savefig('det/'+ str(t) + '.png')

    return SM



def getMM(zai):

    z=zai[0:2]
    a=zai[2:5]
    return elements[int(z)][int(a)].mass

def plotBU(A, name):

    C=abs(A)

    fig, axs = plt.subplots()

    im = axs.imshow(C, cmap='Reds',  norm=LogNorm(vmin=1E-4, vmax=1E+5))
    fig.colorbar(im, orientation='vertical')

    if nucData.energy == 44:

        plt.setp(axs, xticks=[0, 21, 43, 65, 87, 109, 131], xticklabels=['0','22','44','22','44','22','44'])
        plt.setp(axs, yticks=[0, 21, 43, 65, 87, 109, 131], yticklabels=['0','22','44','22','44','22','44'])

    fig.savefig(name+'.png')


    return

def plotCovx(A, name):

    fig, axs = plt.subplots()

    im = axs.imshow(A, cmap='RdYlGn')
    fig.colorbar(im, orientation='vertical')

    #plt.setp(axs, xticks=[0, 21, 43, 65, 87, 109, 131], xticklabels=['0','22','44','22','44','22','44'])
    #plt.setp(axs, yticks=[0, 21, 43, 65, 87, 109, 131], yticklabels=['0','22','44','22','44','22','44'])

    fig.savefig(name+'.png')


    return