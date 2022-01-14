import json
import serpentTools as ST
import matplotlib.pyplot as plt
import numpy as np
from periodictable import elements

PERT = [('922350', '18'), ('922350', '102'), ('922380', '102')]

zailist  = ['922350', '922380', '942390', '80160', '10010']
pertlist = ['tot', 'ela', 'sab', 'inl', '102', '18', 'nxn']

sens = ST.read('covx/GPT/INP_sens0.m')
SSK = sens.sensitivities['keff']
with open('covx/SENS.json') as fp:
    res_sens = json.load(fp)
printKeys = False

def tramezzino(Ns,N,R):

    return float(np.inner(np.array(Ns), np.array(R).dot(np.array(N))))
def getCovx(PERT):

    mtx = PERT[0][2:5] + '_' + PERT[1]

    with open('covx/mtx/' + mtx + '.json') as fp:
        covx = json.load(fp)

    return covx
def plotCovx(A, name):

    fig, axs = plt.subplots()

    im = axs.imshow(A, cmap='Reds')
    fig.colorbar(im, orientation='vertical')
    axs.set(xlabel='Energy groups', title='Covariance matrix for '+name[0]+' MT='+name[1])

    fig.savefig('covx/'+name[0]+' MT='+name[1]+'.png')
def fluxSnap(flux, name, UM, **kwargs):

    x = sens.energies[::-1]

    fig, axs = plt.subplots()

    therm=  [0] + [flux[j] for j in range(len(flux))]

    axs.step(x, therm, 'b', where='pre', label='PTERODAx')
    #axs.set_xlim(0, ene)
    axs.set(xlabel='Energy [MeV]', ylabel='sensitivity' + ' [' + UM +']', title='keff sensitivity to '+name[0]+' MT='+name[1])

    if 'serp' in kwargs.keys():

        flux = kwargs['serp'][::-1]

        therm =  [0] + [flux[j] for j in range(len(flux))]

        axs.step(x, therm, 'r', where='pre', label='SERPENT')

    #axs.set_yscale('log')
    axs.set_xscale('log')
    axs.legend(loc='upper right')
    axs.set_xlim(1E-9, 1E+2)

    fig.savefig('covx/sens_comparison_'+name[0]+' MT='+name[1]+'.png')
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

    return Name
def main(res_sens, PERT):

    print('\nUAM-benchmark keff uncertainties to ND:\n')

    for i in range(len(PERT)):

        zai = PERT[i][0]
        MT = PERT[i][1]

        covx = getCovx(PERT[i])

        plotCovx(covx, PERT[i])

        bun = res_sens['keff'][zai][MT]

        unc = tramezzino(bun[::-1], bun, covx) ** 0.5

        print(getName(zai)+' MT='+MT+':\t'+str(int(unc*1E+5))+' pcm')

        idZai = zailist.index(zai)
        idPert = pertlist.index(MT)

        ssk = SSK[0, idZai, idPert, :, 0]
        fluxSnap(bun, PERT[i], '%/%', serp=ssk)

    if printKeys == True:

        print(sens.sensitivities['keff'].shape)
        print(sens.materials)
        print(sens.zais)
        print(sens.perts)
        print(SSK[0,0,5,:,0])


main(res_sens, PERT)
