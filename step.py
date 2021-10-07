from datetime import datetime
import copy
import numpy as np
from scipy.integrate import odeint
import scipy
import fun
import nucData
import matplotlib.pyplot as plt
import sys
import os
import math
import serpentTools as ST
import onix
from matplotlib.lines import Line2D

start = datetime.now()

plt.rcParams.update({'figure.figsize': (15, 10)})
plt.rcParams.update({'font.size': 22})
plt.rcParams.update({'lines.linewidth': 4})
plt.rcParams.update({'axes.formatter.limits' : (-3,3)})
plt.rcParams.update({'figure.max_open_warning': 60})

##################
### INITIALIZE ###
##################

ZAI = nucData.ZAI
zeroAt = nucData.At
sig = nucData.sig
ene = nucData.ene
reg = 3

where = nucData.where
nodo = nucData.nodo

delta = '\u0394'
square= '\u00b2'

### depletion init ###

P = 1.55E+2   # J/fiss
T = nucData.T
n = nucData.steps
s = 1000

#####################
### PERTURBATION  ###
#####################88

PERT     = ['922350', '922380']                                     # INPUT PERTURBATION NUCLIDE
MT       =  '102'                                                   # INPUT PERTURBATION XS
RESP_NUC =  '942390'                                                # OUTPUT RESPONSE NUCLIDE
RESPONSE =  'keff'                                                  # OUTPUT NUCLIDE OR KEFF
ND       =  False                                                   # SWITCH ND PERTURBATION


pert = 1.01
PERTid = nucData.ZAI.index(PERT[0])
reac = fun.reaz[fun.MT.index(MT)]
respId = nucData.ZAI.index(RESP_NUC)

##############
### DEFINE ###
##############

def zeroStep(At, sig):

    res = fun.Results()

    N = np.array(At) * where

    # flux shape algebra

    k = fun.crushK(N, 0, 0.5, 2)
    #k=fun.findK(N, 0)
    print(k)

    A  = fun.Boltz(N, sig, 0, 1 / k)
    AA = fun.Boltz(N, sig, 0, 1 / k, ribalta=True)


    #fun.plotBU(np.matrix(AA), 'Regions_mtx')
    #fun.plotBU(np.matrix(A), 'Transport_mtx')

    A[1] = np.ones(ene*reg)

    B = np.zeros(len(A[0])).tolist()

    B[1] = 1

    Psi = np.linalg.inv(A).dot(B)

    psi = fun.reshapePsi(Psi)

    #Phi = P / (np.inner(np.inner(At * where, np.array([sig['18'][j][0] for j in [0, 1]])).reshape(6, ), psi.reshape(6, )))
    #Phi = P / (np.inner(np.inner(At * where, np.array([sig['18'][j][0]*sig['v'][0] for j in [0, 1]])).reshape(6, ), psi.reshape(6, )))
    Phi = P / (np.inner(np.inner(At * where, np.array([sig['18'][j][0] * sig['v'] for j in range(ene)])).ravel(), psi.ravel()))

    ### burnup ###

    rr = fun.rr(nucData.sig, psi, Phi, 0)

    C = fun.Bateman(rr)

    # results

    res.flux.append(Psi)

    res.phi.append(Phi)

    res.keff.append(k)

    res.M.append(C)

    return res

def directStep(At, sig):

    res = zeroStep(At,sig)

    res.comp.append(At)


    for i in range(n-1):

        # depletion ODE

        dt = (nucData.tempo[i+1]-nucData.tempo[i])*24*3600

        old_stdout = sys.stdout  # backup current stdout
        sys.stdout = open(os.devnull, "w")
        #At1 = odeint(fun.ODE2, At, np.linspace(0, T / n, s), res.M[i]).tolist()
        At1 = onix.salameche.CRAM16(np.matrix(res.M[i])*(dt), np.array(At))
        sys.stdout = old_stdout  # reset old stdout


        res.comp.append(At1.tolist())
        At = At1

        N = np.array(At) * where

        # updated flux shape algebra

        k = fun.crushK(N, i, 0.3, 3)
        #k=fun.findK(N, i)

        A = fun.Boltz(N, sig, i, 1 / k)

        A[1] = np.ones(ene*reg)

        B = np.zeros(len(A[0])).tolist()

        B[1] = 1

        Psi = np.linalg.inv(A).dot(B)

        # power normalization

        psi = fun.reshapePsi(Psi)


        Phi = P / (np.inner(np.inner(At * where, np.array([sig['18'][j][i]* sig['v'] for j in range(ene)])).ravel(), psi.ravel()))

        # updated transmuation matrix

        rr = fun.rr(fun.sig, psi, Phi, i)

        C = fun.Bateman(rr)

        # results

        res.flux.append(Psi)

        res.phi.append(Phi)

        res.keff.append(k)

        res.M.append(C)

        #print(str(i + 1) + '\t%.5f' % k)

        print('\t'+str(math.ceil(i/ n *100))+"% complete", end='\r')
        sys.stdout.flush()

    print("100% complete                ", end='\r')


    return res

def pertBlock(res, **kwargs):

    j = 0

    res.pert['atoms'] = []
    res.pert['keff'] = []

    ND = kwargs['ND']

    if ND == False:

        for z in PERT:
            zeroAt = nucData.At.copy()

            pertId = nucData.ZAI.index(z)

            print('\n\nPerturbed calculation ' + str(j + 1) + ' (' + nucData.ZAI[pertId] + ')\n')

            Bdir = zeroAt[pertId]
            Bper = zeroAt[pertId] * pert

            zeroAt[pertId] = Bper

            pertRes = directStep(zeroAt, sig)

            Npert = pertRes.comp[-1]
            No = res.comp[-1]
            res.pert['atoms'].append(Npert[respId] - No[respId])
            res.pert['keff'].append((pertRes.keff[-1] - res.keff[-1]))

            j = j + 1
            # printRes(pertRes, keff=True, comp = False, phi = False, flux = False)

    elif ND == True:

        for e in range(ene):

            zeroAt = nucData.At.copy()
            pertId = PERTid

            print('\n\nPerturbed calculation ' + str(j + 1) + '/'+str(ene)+' (' + nucData.ZAI[pertId] + ')\n')

            sigma = copy.deepcopy(sig)

            for t in range(nucData.steps):

                sigma[MT][e][t][pertId] = sigma[MT][e][t][pertId] * pert
                sigma['removal'][e][t][pertId] += sigma[MT][e][t][pertId] * (pert-1)

            pertRes = directStep(zeroAt, sigma)

            Npert = pertRes.comp[-1]
            No = res.comp[-1]
            res.pert['atoms'].append(Npert[respId] - No[respId])
            res.pert['keff'].append((pertRes.keff[-1] - res.keff[-1]))

            j = j + 1

def adjoStep(res, **kwargs):

    adjoRes = fun.Results()

    adjoRes.source.append(np.zeros(len(res.flux[0].tolist())).tolist())
    adjoRes.flux.append(np.zeros(len(res.flux[0].tolist())).tolist())

    Ns = np.zeros(len(res.comp[0]))
    Ns1 = Ns.copy()

    sk=np.zeros(len(res.comp[0])).tolist()
    ind_1=np.zeros(len(res.comp[0])).tolist()
    ind_2=np.zeros(len(res.comp[0])).tolist()
    dir_1=np.zeros(len(res.comp[0])).tolist()

    S = np.zeros(ene).tolist()

    resp = kwargs['resp']

    if resp == 'keff':

        ind_3 = np.zeros(len(res.comp[0])).tolist()
        SS=Ns.copy()

        for i in range(n-1):

            v = -2 - i

            Psi = res.flux[v]

            Phi = float(res.phi[v])

            C = res.M[v+1]

            k = res.keff[v]

            N = res.comp[v+1]
            No = res.comp[v]

            psi = fun.reshapePsi(Psi)

            rr = fun.rr(fun.sig, psi, Phi, v)

            # homogeneous adjoint

            A = fun.Boltz(No * where, sig, v, 1 / k).transpose()

            A[1] = Psi

            B = np.zeros(len(A[0])).tolist()

            B[1] = 1

            Gh = np.linalg.inv(A).dot(B)

            adjoRes.homo.append(Gh)

            # keff sens

            ssk = [(fun.kSens(Gh, Psi, No, k, j, v)) for j in range(len(No))]
            #ssk[nucData.ZAI.index('922350')]=-ssk[nucData.ZAI.index('922350')]

            dR = Psi.copy() * 0
            dR2 = Psi.copy() * 0
            ai = 0

            if  i == 0:
                Ns  = np.array(ssk)
                Nss = Ns.copy()

                dR = fun.dR(Psi, Gh, N, k, v)
                ai = (Psi).dot(dR)

                dR2 = fun.dR2(Psi, Gh, N, k, v)
                dR2[1] = 0

                SS = Ns.copy()
                sk = [sk[j] + ssk[j] for j in range(len(N))]

            adjoRes.comp.append(Nss)

            adjoRes.ind.append([ind_2, dir_1, ind_1, ind_3, sk, Nss])

            dt = (nucData.tempo[v+1] - nucData.tempo[v]) * 24 * 3600

            SS = Ns1.copy()

            old_stdout = sys.stdout  # backup current stdout
            sys.stdout = open(os.devnull, "w")
            Ns1 = onix.salameche.CRAM16(np.matrix(C).transpose() * (dt), np.array(Ns))
            sys.stdout = old_stdout  # reset old stdout

            # adjoint power normalization

            PL = fun.updatePL(fun.pl, rr )
            R = fun.Bateman(rr)#-fun.onixD(PL)
            Ps = (fun.I([Ns1, SS], [No, N],  R, dt) + ai) / P

            adjoRes.pow.append(Ps)

            # adjoint source

            B = -dR -np.array(fun.Qs([No, N], [Ns1, SS], Ps, Phi, v, dt))

            adjoRes.source.append(np.array(B))

            # adjoflux-shape algebra

            A = fun.Boltz(No * where, sig, v, 1 / k).transpose()

            #A[1]=fun.boltzF(No*where, sig, v).dot(Psi)
            A[1] = Gh

            B[1] = 0

            Gp = np.linalg.inv(A).dot(np.array(B))

            b = (fun.boltzF(No * where, sig, v).dot(Psi).dot(Gp)) / (fun.boltzF(No * where, sig, v).dot(Psi).dot(Gh))

            G = Gp - b * Gh

            adjoRes.flux.append(G)

            A = fun.Boltz(No * where, sig, v, 1 / k)
            A[1] = Psi
            G2 = np.linalg.inv(A).dot(np.array(-dR2))

            # step condition

            lam = 1 / k

            Beta  = fun.beta(Psi, lam, No, v)
            Beta2 = fun.beta(Gh, lam, No, v)

            Pi = rr['fission'].copy()*sig['v']

            ind_1 = [ind_1[j] + np.inner(G, Beta[j]) for j in range(Beta.shape[0])]
            ind_2 = [ind_2[j] - Ps * Pi[j] for j in range(Beta.shape[0])]
            ind_3 = [ind_3[j] + np.inner(G2, Beta2[j]) for j in range(Beta.shape[0])]
            dir_1 = [dir_1[j] + Ns1[j] - Ns[j] for j in range(Beta.shape[0])]

            # Ns = [Ns1[j] - (np.inner(G, Beta[j]) + (Ps * Pi[j])) + (fun.kSens(Gh, Psi,  N, k , j, v)) for j in range(Beta.shape[0])]

            Ns = [Ns1[j] + (np.inner(G, Beta[j]) - Ps * Pi[j]  + np.inner(G2, Beta2[j])) for j in range(Beta.shape[0])]
            Nss = [Ns[j] for j in range(Beta.shape[0])]

            # adjoRes.ind.append([ind_2,dir_1,ind_1, sk, Nss])

            print('\t' + str(math.ceil(i / n * 100)) + "% complete", end='\r')
            sys.stdout.flush()

        adjoRes.ind.append([ind_2, dir_1, ind_1, ind_3, sk, Ns])

    elif resp == 'nuclide' and 'xs' not in kwargs.keys():

        Ns[respId] = 1
        Nss = Ns.copy()
        Ns1 = Ns.copy()

        SS=Ns.copy()

        for i in range(n-1):

            v = -2-i

            Psi = res.flux[v]

            Phi = float(res.phi[v])

            C = res.M[v]

            k = res.keff[v]

            N = res.comp[v+1]
            No = res.comp[v]

            psi = fun.reshapePsi(Psi)

            rr = fun.rr(fun.sig, psi, Phi,v)

            # homogeneous adjoint

            A = fun.Boltz(No * where, sig, v, 1/k).transpose()

            A[1] = Psi

            B = np.zeros(len(A[0])).tolist()

            B[1] = 1

            Gh = np.linalg.inv(A).dot(B)

            adjoRes.homo.append(Gh)

            # keff sens

            adjoRes.comp.append(Ns)

            adjoRes.ind.append([ind_2,dir_1,ind_1, Ns])

            dt = (nucData.tempo[v+1] - nucData.tempo[v]) * 24 * 3600

            SS=Ns1.copy()

            old_stdout = sys.stdout  # backup current stdout
            sys.stdout = open(os.devnull, "w")
            #Ns1 = odeint(fun.ODE2b, Ns, time, C).tolist()[::-1]
            Ns1 = onix.salameche.CRAM16(np.matrix(C).transpose()*(dt), np.array(Ns))
            sys.stdout = old_stdout  # reset old stdout

            # adjoint power normalization

            PL=fun.updatePL(fun.pl, rr)
            R=fun.onixR(PL)
            Ps=fun.I([Ns1,SS],[No,N],R, dt)/P

            adjoRes.pow.append(Ps)

            # adjoint source

            B = -np.array(fun.Qs([No,N],[Ns1,SS], Ps, Phi, v, dt))

            adjoRes.source.append(np.array(B))

            # adjoflux-shape algebra

            A = fun.Boltz(No * where, sig, v, 1/k).transpose()

            A[1]=fun.boltzF(No*where, sig, v).dot(Psi)
            #A[1]=Gh

            B[1]=0

            Gp = np.linalg.inv(A).dot(np.array(B))

            b = (fun.boltzF(No*where, sig, v).dot(Psi).dot(Gp))/(fun.boltzF(No*where, sig, v).dot(Psi).dot(Gh))

            G = Gp - b*Gh*0

            adjoRes.flux.append(G)

            # step condition

            lam = 1 / k

            Beta = fun.beta(Psi, lam, No, v)

            Pi = rr['fission'].copy()*sig['v']

            ind_1= [ind_1[j] + np.inner(G, Beta[j]) for j in range(Beta.shape[0])]
            ind_2= [ind_2[j] - Ps * Pi[j] for j in range(Beta.shape[0])]
            dir_1= [dir_1[j] + Ns1[j]-Ns[j] for j in range(Beta.shape[0])]

            Ns = [Ns1[j]  + (np.inner(G, Beta[j]) - (Ps * Pi[j]))  for j in range(Beta.shape[0])]

            print('\t'+str(math.ceil(i / n * 100)) + "% complete", end='\r')
            sys.stdout.flush()

        adjoRes.ind.append([ind_2, dir_1, ind_1, Ns])

    elif resp == 'nuclide' and 'xs' in kwargs.keys() :

        ind_1 = np.zeros(ene).tolist()
        ind_2 = np.zeros(ene).tolist()
        dir_1 = np.zeros(ene).tolist()

        RESP = res.comp[-1][respId]

        xs_pert = kwargs['xs']

        Ns[respId] = 1
        Nss = Ns

        SS=Ns.copy()

        for i in range(n-1):

            v = -2-i

            Psi = res.flux[v]

            Phi = float(res.phi[v])

            C = res.M[v]

            k = res.keff[v]

            N = res.comp[v+1]
            No = res.comp[v]

            psi = fun.reshapePsi(Psi)

            rr = fun.rr(fun.sig, psi, Phi,v)

            # homogeneous adjoint

            A = fun.Boltz(No * where, sig, v, 1/k).transpose()

            A[1] = Psi

            B = np.zeros(len(A[0])).tolist()

            B[1] = 1

            Gh = np.linalg.inv(A).dot(B)

            adjoRes.homo.append(Gh)

            # keff sens

            adjoRes.comp.append(Ns)

            adjoRes.ind.append([dir_1, ind_1, ind_2, S])

            dt = (nucData.tempo[v+1] - nucData.tempo[v]) * 24 * 3600

            SS=Ns1.copy()

            old_stdout = sys.stdout  # backup current stdout
            sys.stdout = open(os.devnull, "w")
            #Ns1 = odeint(fun.ODE2b, Ns, time, C).tolist()[::-1]
            Ns1 = onix.salameche.CRAM16(np.matrix(C).transpose()*(dt), np.array(Ns))
            sys.stdout = old_stdout  # reset old stdout

            # adjoint power normalization

            PL=fun.updatePL(fun.pl, rr)
            R=fun.Bateman(rr)
            Ps=fun.I([Ns1,SS],[No,N],R, dt)/P

            adjoRes.pow.append(Ps)

            # adjoint source

            B = -np.array(fun.Qs([No,N],[Ns1,SS], Ps, Phi, v, dt))

            adjoRes.source.append(np.array(B))

            # adjoflux-shape algebra

            A = fun.Boltz(No * where, sig, v, 1/k).transpose()

            A[1]=fun.boltzF(No*where, sig, v).dot(Psi)
            #A[1]=Gh

            B[1]=0

            Gp = np.linalg.inv(A).dot(np.array(B))

            b = (fun.boltzF(No*where, sig, v).dot(Psi).dot(Gp))/(fun.boltzF(No*where, sig, v).dot(Psi).dot(Gh))

            G = Gp - b*Gh*0

            adjoRes.flux.append(G)

            # step condition

            lam = 1 / k

            Beta = fun.beta(Psi, lam, No, v)

            Pi = rr['fission'].copy()*sig['v']

            Ns = [Ns1[j]  + (np.inner(G, Beta[j]) - (Ps * Pi[j]))  for j in range(Beta.shape[0])]

            # SENSITIVITY

            Beta_sig = fun.betaSig(Psi, lam, xs_pert, No, G, PERTid, v)*Phi
            Bate_sig = fun.bateSig(Psi, Phi, xs_pert, [No,N],[Ns1,SS], PERTid, v, dt)
            Pi_sig   = -fun.PiSig(Psi, Phi, xs_pert, No, PERTid, v)*Ps

            # *sig['18'][e][v][PERTid]/RESP

            S = [S[e] + (Bate_sig[e] + Beta_sig[e] + Pi_sig[e])*(sig[xs_pert][e][v][PERTid]/RESP)**0 for e in range(ene)]

            ind_1= [ind_1[e] + Beta_sig[e]*(sig[xs_pert][e][v][PERTid]/RESP)**0 for e in range(ene)]
            ind_2= [ind_2[e] + Pi_sig[e]*(sig[xs_pert][e][v][PERTid]/RESP)**0 for e in range(ene)]
            dir_1= [dir_1[e] + Bate_sig[e]*(sig[xs_pert][e][v][PERTid]/RESP)**0 for e in range(ene)]

            print('\t'+str(math.ceil(i / n * 100)) + "% complete", end='\r')
            sys.stdout.flush()

        adjoRes.ind.append([dir_1, ind_1, ind_2, S])

    print("100% complete                ", end='\r'+'\n\n' )

    adjoRes.comp = adjoRes.comp[::-1]
    adjoRes.pow = adjoRes.pow[::-1]
    adjoRes.flux = adjoRes.flux[::-1]
    adjoRes.source = adjoRes.source[::-1]
    adjoRes.ind = adjoRes.ind[::-1]

    return adjoRes


def printRes(res, **kwargs):

    print('\n\nResults')

    for key,value in kwargs.items():

        if value == True:

            print('\n' + key)

            try:

                getattr(res, key)[0][0][1]

            except:

                #print('%.5f' % getattr(res, key)[0])
                #print('%.5f' % getattr(res, key)[-1])
                print(getattr(res, key)[0])
                print(getattr(res, key)[-1])


            else:

                print(getattr(res, key)[0][0])
                print(getattr(res, key)[-1][-1])

        def pertBlock():

            j = 0

            res.pert['atoms'] = []
            res.pert['keff'] = []

            for z in PERT:
                zeroAt = nucData.At.copy()

                pertId = nucData.ZAI.index(z)

                print('\n\nPerturbed calculation ' + str(j + 1) + ' (' + nucData.ZAI[pertId] + ')\n')

                Bdir = zeroAt[pertId]
                Bper = zeroAt[pertId] * pert

                zeroAt[pertId] = Bper

                pertRes = directStep(zeroAt)

                Npert = pertRes.comp[-1]
                No = res.comp[-1]
                res.pert['atoms'].append(Npert[respId] - No[respId])
                res.pert['keff'].append((pertRes.keff[-1] - res.keff[-1]))

                j = j + 1
                # printRes(pertRes, keff=True, comp = False, phi = False, flux = False)

def massPlot(res):

    dep = ST.read(nucData.file+'/REP_dep.m')

    isoFuel=nucData.zaiOld
    isoElse=[]

    prezz = []

    mat = nucData.MAT

    #x = np.linspace(0, T, n+1) / 3600 / 24

    x = nucData.tempo
    x2 = dep.days

    for i in range(len(isoFuel)):

        k=nucData.ZAI.index(isoFuel[i])
        name = nucData.nuc[k].name

        y1 = np.array([a[k] for a in res.comp]) * nucData.getMM(isoFuel[i]) / 6.022E+23

        fig, ax1 = plt.subplots()

        ax1.set(xlabel='BU (days)', ylabel='Nuclide mass [g]', title=name)
        ax1.grid()

        ax1.plot(x, y1, 'b', label = 'SIBYL')

        if ZAI[k] not in ['922390', '932390', '942400']:

            y2 = dep.materials[mat[nucData.nuc[k].id]].getValues('days','mdens', zai=int(ZAI[k]))[0]*nucData.nuc[k].vol
            ax1.plot(x2, y2, 'r', label = 'SERPENT')

        ax1.legend(loc='best')

        #jump=abs(max(y1)-min(y1))
        #ax1.set_ylim(min(y1)-0.1*jump, min(y1)+1.1*jump)

        fig.savefig('masse/'+ isoFuel[i] + '.png')

    for i in range(len(isoElse)):

        k=nucData.ZAI.index(isoElse[i])

        y1 = np.array([a[k] for a in res.comp]) * nucData.getMM(isoElse[i]) / 6.022E+23

        y2 = dep.materials[mat[nucData.nuc[k].id]].getValues('days','mdens', zai=int(ZAI[k]))[0]*nucData.nuc[k].vol

        fig, ax1 = plt.subplots()

        if nucData.nuc[k].zai in prezz:

            y2 = np.ones(len(x2))*y2[0]

            ax1.set_ylim(y2[0]*0.5, y2[0]*1.1)

        if nucData.nuc[k].xmat:

            y2 += dep.materials[nucData.nuc[k].xmat].getValues('days','mdens', zai=int(ZAI[k]))[0]*nucData.nuc[k].vol


        ax1.set(xlabel='BU (days)', ylabel='Nuclide mass [g]', title=nucData.nuc[k].name)
        ax1.grid()

        ax1.plot(x, y1, 'b', label = 'SIBYL')
        ax1.plot(x2, y2, 'r', label = 'SERPENT')
        ax1.legend(loc='center right')

        #jump=abs(max(y1)-min(y1))
        #ax1.set_ylim(min(y1)-0.1*jump, min(y1)+1.1*jump)

        fig.savefig('masse/'+ nucData.nuc[k].name + '.png')

def paraPlot(para, name, **kwargs):

    #x = np.linspace(0, T/24/3600, len(para))
    y1 = np.array([a for a in para ])

    fig, ax1 = plt.subplots()

    if 'serp' in kwargs.keys():
        if kwargs['serp'] == True:

            ax2 = ax1.twinx()
            ax2.set(xlabel='BU (days)', ylabel='SERPENT '+name)

            x2 = kwargs['paraserp'][0]
            y2 = np.array([a for a in kwargs['paraserp'][1] ])

            jump = abs(max(y2) - min(y2))
            ax1.set_ylim(max(y1) - jump * 1.1, max(y1) * 1.005)
            ax2.set_ylim(max(y2) - jump * 1.1, max(y2) * 1.005)

            ax2.plot(x2, y2, 'r', label = 'SERPENT')
            ax2.legend(loc='upper center')


    x = nucData.tempo[:len(y1)]

    ax1.set(xlabel='BU (days)', ylabel='SIBYL '+name, title = name + ' BU evolution')
    ax1.grid()

    ax1.plot(x, y1, 'b', label = 'SIBYL')
    ax1.legend(loc='upper right')

    fig.savefig(name + '.png')

def adjoPlot(res, adjoRes, **kwargs):

    x = nucData.tempo
    g = []
    X = []
    c = 1
    pcm = 1

    lab = ['power', 'at. density', 'flux shape']
    col = ['brown', 'green', 'orange']

    resp = kwargs['resp']

    if resp=='nuclide':

        c = nucData.getMM(ZAI[respId]) / 6.022E+23
        resp = nucData.nuc[respId].name
        sibyl = res.pert['atoms']
        title = resp + ' mass change [g]'


    if resp == 'keff':
        sibyl = res.pert['keff']
        lab.extend(['flux adjoint','k-sens'])
        col.extend(['gold','m'])
        pcm = 1E+5
        title = 'reactivity [pcm]'

    erro = []

    j = 0


    lab.extend(['TOTAL'])
    col.extend(['b--'])

    for z in nucData.zaiOld:

        k = nucData.ZAI.index(z)

        nn = zeroAt[k] * (pert-1)

        if nn == 0:

            nn = res.comp[-1][k]*0.01

        fig, ax1 = plt.subplots()
        ax1.grid()

        name = nucData.nuc[k].name

        ax1.set(xlabel='BU (days)', ylabel='Gross Sensitivity', title=name+' contributions to EOL '+title+'  (rel. error)')

        i = 0
        b = 0

        for a in np.array(adjoRes.ind).transpose()[k]:
            ax1.plot(x, a*nn*c, col[i], label=lab[i])
            i += 1
            b += a

        # b=np.array(adjoRes.ind).transpose()[k][0]-sum(np.array(adjoRes.ind).transpose()[k][1:4])
        # ax1.plot(x, b, 'green', label='direct')

        if z in PERT:
            kk = PERT.index(z)
            y = np.ones(len(x)) * sibyl[kk] * c
            ax1.plot(x, y, 'r--', markersize=15, label='pert')
            ax1.vlines(0, np.array(adjoRes.ind).transpose()[k][-1][0]*nn*c, y[0], linewidth=10, label='error')

            g.append([np.array(adjoRes.ind)[0][j][k] * c * nn for j in range(len(lab))] + [sibyl[kk] *c ])
            X.append(name)

            e =(np.array(adjoRes.ind).transpose()[k][-1][0]*nn*c - sibyl[kk]*c) / sibyl[kk]/c * 100
            erro.append(e)
            string = '%.1f' % abs(e) + '%'
            ax1.annotate(string, (0, y[0]*(1 - 0.05*np.sign(e))))

            j = j + 1

        ax1.legend(loc='best')

        fig.savefig('adjomasse/' + z + '.png')

    col[-1]='b'
    lab.append('SIBYL')
    col.append('r')

    fig3, ax3 = plt.subplots()
    xx = np.arange(len(PERT)) * 5

    ax3.grid()

    # ax3.set_yscale('symlog', linthresh=0.01)
    i = 0
    w = 0.5

    for a in np.array(g).transpose():
        ax3.bar(xx - w * len(lab) / 2 + w * i, a*pcm, width=w, color=col[i], label=lab[i])
        i = i + 1

    for j in range(len(erro)):

        string = '%.1f' % (erro[j]) + '%'

        if pcm == 1E+5:
            string = '%.1f' % (erro[j]) + '% ('+str(abs(int(erro[j]*g[j][-1]*pcm/100)))+' pcm)'

        if abs(erro[j]) > 100:
            string = '      NR'

        if g[j][-2] > 0:
            ax3.annotate(string, (xx[j] -2*w, max(g[j]) * 1.1 * pcm))

        else:
            ax3.annotate(string, (xx[j] -2*w, min(g[j]) * 1.6 * pcm))

    ax3.set_yscale('symlog', linthresh=pcm*min([abs(a[-1]+a[0]) for a in np.array(g).transpose()+1E-6]))
    ax3.set(xlabel='\n1% Perturbed BOL isotope', ylabel='EOL '+title,
            title='PTERODAx contributions to EOL '+title+'  (rel. error)')
    ax3.legend(loc='best')
    ax3.set_xticks(xx)
    ax3.set_xticklabels(X)
    ax3.tick_params(axis='both', which='major', labelsize=18)
    ax3.margins(0.1)

    fig3.savefig('adjohisto_'+resp)

def fluxPlot(flux, name, UM, **kwargs):

    if 'phi' in kwargs.keys():
        phi = kwargs['phi']

    else:
        phi=np.ones(len(flux))
        #flux.append(flux[-1])

    t = nucData.tempo[:len(flux)]
    #t.append(50)

    fig, axs = plt.subplots(1, 3, sharey=True)
    # Remove vertical space between axes
    fig.subplots_adjust(wspace=0)

    # Plot each graph, and manually set the y tick values
    for i in range(int(len(flux[0])/2)):

        therm=[phi[j]*flux[j][i+3] for j in range(len(flux))]
        fast=[phi[j]*flux[j][i] for j in range(len(flux))]

        axs[i].plot(t, therm, 'b', label='thermal')
        axs[i].plot(t, fast, 'r', label='fast')
        axs[i].set_xlim(0, nucData.giorni)
        axs[i].set(xlabel='BU (days)', title=nucData.reg[i].name)

    if 'serp' in kwargs.keys():

        flux = kwargs['serp']

        t = nucData.time

        for i in range(int(len(flux[0]) / 2)):
            therm = [ flux[j][i + 3] for j in range(len(flux))]
            fast = [ flux[j][i] for j in range(len(flux))]

            axs[i].plot(t, therm, 'b:', label='SERP thermal')
            axs[i].plot(t, fast, 'r:', label='SERP fast')

    axs[2].legend(loc='upper right')

    #axs[1].set(xlabel='BU (days)', title='Neutron flux evolution in central-fuel-reflector regions \n')
    axs[0].set( ylabel=name + ' [' + UM +']')


    fig.savefig('flux/' + name + '.png')

def fluxSnap(flux, name, UM, **kwargs):

    if 'phi' in kwargs.keys():
        phi = kwargs['phi']

    else:
        phi=np.ones(len(flux))
        #flux.append(flux[-1])

    x = nucData.grid

    fig, axs = plt.subplots()

    therm=[phi*flux[j] for j in range(len(flux))] + [0]

    axs.step(x, therm, 'b', where='pre', label='SIBYL')
    #axs.set_xlim(0, ene)
    axs.set(xlabel='Energy [MeV]', ylabel=name + ' [' + UM +']', title='Flux spectrum')

    if 'serp' in kwargs.keys():

        flux = kwargs['serp']

        therm = [flux[j] for j in range(len(flux))] + [0]

        axs.step(x, therm, 'r', where='pre', label='SERPENT')

    #axs.set_yscale('log')
    axs.set_xscale('log')
    axs.legend(loc='upper right')
    axs.set_xlim(1E-9, 1E+1)

    fig.savefig('flux/' + name + '_snap.png')

def bunSnap(resu, res, name, xs):

    x = nucData.grid

    fig, axs = plt.subplots()

    lin = [(0, ()), (0, ()), (0, ()), (0, (2,2))]
    lab = ['at. density', 'flux shape', 'power', 'TOTAL']
    col = ['green', 'orange', 'brown', 'blue']

    j = 0

    for flux in resu:

        y = flux + [0]
        axs.step(x, y, col[j], linestyle= lin[j], where = 'pre', label=lab[j])
        tit = 'EOL '+nucData.nuc[respId].name+' Sensitivity to '+nucData.nuc[PERTid].name+' '+xs+' cross section\n'
        axs.set(xlabel='Energy [MeV]', ylabel='sensitivity', title=tit)
        #axs.set_yscale('log')
        axs.set_xscale('log')

        j+=1

    RESP = res.comp[-1][respId]
    #y2 = [res.pert['atoms'][e]/RESP/(pert-1) for e in range(ene)] + [0]
    y2 = [res.pert['atoms'][e]/(sig[MT][e][nodo][PERTid]*(pert-1)) for e in range(ene)] + [0]
    axs.step(x, y2, 'red', linestyle=lin[-1], where='pre', label='SIBYL')
    axs.legend(loc='best')
    axs.set_xlim(1E-9, 1E+1)

    fig.savefig('flux/sensitivity_to_'+name+'_snap.png')

def printTime():
    end = datetime.now()
    diff = end - start
    chrono = divmod(diff.total_seconds(), 60)

    print('\nExecution time: ' + str(int(chrono[0])) + ':' + '%02d' % (float(chrono[1]),) + ':%02d\t\t(min:sec:dec)\n' % (
            float(chrono[1]) * 100 % 100,))

### MAIN ###

def main(**kwargs):

    print('\nNominal calculation\n')

    res=directStep(zeroAt, sig)

    printRes(res, keff=True, comp = False, phi = False, flux = False)

    paraPlot(res.keff, 'keff', serp=True, paraserp=[nucData.time, nucData.keff])

    massPlot(res)

    #fluxPlot(res.flux,  'Neutron Flux', 'n/cm'+square+'s', phi=res.phi, serp=nucData.Flux)
    psi = fun.reshapePsi(res.flux[0])
    serpsi = fun.reshapePsi(nucData.Flux[0])
    fluxSnap(psi[0],  'Fuel Flux', 'n/cm'+square+'s', phi=res.phi[0], serp=serpsi[0])
    fluxSnap(psi[2],  'Reflector Flux', 'n/cm'+square+'s', phi=res.phi[0], serp=serpsi[2])

    if kwargs['ptero'] != None:

        ### perturbed solution ###

        pertBlock(res, ND=kwargs['ND'])

        ### adjoint solution ###

        resp = kwargs['ptero']

        if kwargs['ND'] == True :

            print('\n\nAdjoint calculation\n')

            adjoRes=adjoStep(res, resp=resp, xs=MT)

            bunSnap(adjoRes.ind[0], res, PERT[0], reac)

        else:

            print('\n\nAdjoint calculation\n')

            adjoRes=adjoStep(res, resp=resp)

        adjoPrint = True
        if adjoPrint == True:

            #print('\nAdjoint res\n')

            #print(adjoRes.comp[0][0])
            #print(adjoRes.comp[-1][-1])

            #paraPlot(nucData.sig['cap'][0][0], 'xs_U5')
            #paraPlot(nucData.keff, 'SERPENT keff')

            flu = adjoRes

            fluxPlot(flu.flux, 'Adjoint Flux', '')
            fluxPlot(flu.source, 'Adjoint Source', '')
            fluxPlot(flu.homo, 'Homogeneous adjoint Flux', '')

        adjoPlot(res, adjoRes, resp=resp)

main(ptero=RESPONSE, ND=ND)

### chrono ###

printTime()

sys.exit(0)
raise SystemExit

