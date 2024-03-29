import json
import serpentTools as ST
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from periodictable import elements

PERT = [('922350', '18'), ('922350', '452'), ('922350', '102'), ('922380', '18'), ('922380', '452'), ('922380', '102')]
#PERT = [('922350', '18'), ('922350', '452'), ('922350', '102'), ('922380', '18'), ('922380', '452'), ('922380', '102')]

zailist  = ['10010',  '80160','922350', '922360', '922380', ]
pertlist = ['tot', 'ela', 'sab', 'inl', '102', '18', '452', 'nxn']

sens = ST.read('adrien/inp_sens0.m')
SSK = sens.sensitivities['keff']
with open('COVX/SENS_het_bigrun.json') as fp:
    res_sens = json.load(fp)
printKeys = False

resp = 'keff'



def tramezzino(Ns,N,R):

    return float(np.inner(np.array(Ns), np.array(R).dot(np.array(N))))
def getCovx(PERT):

    z  = PERT[0][:-1]
    MT = PERT[1]

    with open('COVX/mtx.json') as fp:
        cvx = json.load(fp)

    mtx = cvx[z][MT]

    return mtx
def fluxSnap(flux, name, UM, **kwargs):

    x = sens.energies[::-1]

    fig, axs = plt.subplots()

    therm=  [0] + [flux[j] for j in range(len(flux))]

    axs.step(x, therm, 'b', where='pre', label='PTERODAx DPT')
    #axs.set_xlim(0, ene)
    axs.set(xlabel='Energy [MeV]', ylabel='sensitivity' + ' [' + UM +']', title='keff sensitivity to '+name[0]+' MT='+name[1])

    if 'serp' in kwargs.keys():

        flux = kwargs['serp'][::-1]

        therm =  [0] + [flux[j] for j in range(len(flux))]

        axs.step(x, therm, 'r', where='pre', label='SERPENT GPT')

    #axs.set_yscale('log')
    axs.set_xscale('log')
    axs.legend(loc='upper right')
    axs.set_xlim(1E-9, 1E+2)

    fig.savefig('COVX/PLOT/sens_comparison_'+name[0]+' MT='+name[1]+'.png')
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

    tot = 0
    pcm = 1
    udm = ' %'

    if resp == 'keff':
        pcm = 1E+5
        udm = ' pcm'

    for i in range(len(PERT)):

        zai = PERT[i][0]
        MT = PERT[i][1]

        covx = getCovx(PERT[i])

        idZai = zailist.index(zai)
        idPert = pertlist.index(MT)

        bun = SSK[0, idZai, idPert, :, 0][::-1]

        unc = abs(tramezzino(bun[::-1], bun, covx)) ** 0.5

        tot += unc**2

        print(SSK[0, idZai, idPert, :, 0][::-1])
        print(getName(zai)+' MT='+MT+':\t'+str(unc*pcm)+udm)




    print('\ntotal uncertainty: \t'+str(int(tot**0.5*pcm))+udm+'\n')

    if printKeys == True:

        print(sens.sensitivities['keff'].shape)
        print(sens.materials)
        print(sens.zais)
        print(sens.perts)
        print(SSK[0,0,5,:,0])

main(res_sens, PERT)



