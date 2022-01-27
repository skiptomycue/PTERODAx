import os
import configparser as CP
import json
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import nucData

plt.rcParams.update({'font.size': 22})
plt.rcParams.update({'lines.linewidth': 4})
plt.rcParams['lines.markersize'] = 10
plt.rcParams.update({'figure.figsize': (15, 10)})
plt.rcParams.update({'figure.max_open_warning': 60})
plt.rcParams.update({'axes.formatter.limits' : (-3,3)})

def resDict(key1, key2, name):

    res_sens={}

    for a in key1:
        res_sens[a]= {}

        for b in key2:
            res_sens[a][b] = {}


    with open(name+'.json', 'w') as fp:
        json.dump(res_sens, fp)

def spesaSENS(RESP, PERT, MT):

    resDict(RESP, PERT, 'COVX/SENS')

    start = datetime.now()

    i   = 0
    tot = len(RESP)*len(PERT)*len(MT)

    for r in RESP:

        for p in PERT:

            for m in MT:

                print('\nCalculation '+str(i+1)+'/'+str(tot))

                resp = 'nuclide'
                s    = '%s' % r

                if r == 'keff':

                    resp = 'keff'
                    s    = '942390'

                config.set('sibyl', 'energy', '44')
                config.set('sibyl', 'PASSI', '125')
                config.set('sibyl', 'hetSteps', 'True')

                config.set('pterodax', 'direct', 'True')
                config.set('pterodax', 'PERT_NUC', p)
                config.set('pterodax', 'RESP_NUC', s)
                config.set('pterodax', 'RESPONSE', resp)
                config.set('pterodax', 'MT', m)

                with open(r"configfile.ini", 'w') as configfile:
                    config.write(configfile)

                os.system('python step.py > /dev/null')

                end = datetime.now()
                diff = end - start
                chrono = divmod(diff.total_seconds(), 60)

                print('Done!')
                print('Running time: ' + str(int(chrono[0])) + ':' + '%02d' % (
                float(chrono[1]),) + ':%02d\t\t(min:sec:dec)\n' % (
                          float(chrono[1]) * 100 % 100,))

                i += 1

    os.system('python covx.py > /dev/null')

def hetStudy(maglie, resp, pert):

    start = datetime.now()

    resDict(maglie, [], model+'/HET')

    config.set('pterodax', 'PERT_NUC', pert)
    config.set('pterodax', 'RESP_NUC', '942390')
    config.set('pterodax', 'RESPONSE', resp)
    config.set('pterodax', 'ND', 'False')
    config.set('pterodax', 'direct', 'True')

    i   = 0
    tot = len(maglie)*2

    for h in ['True', 'False']:

        config.set('sibyl', 'hetSteps', h)

        for m in maglie:

            print('\nCalculation ' + str(i + 1) + '/' + str(tot))

            config.set('sibyl', 'PASSI', str(m))

            with open(r"configfile.ini", 'w') as configfile:
                config.write(configfile)

            os.system('python step.py > /dev/null')

            end = datetime.now()
            diff = end - start
            chrono = divmod(diff.total_seconds(), 60)

            print('Done!')
            print('Running time: ' + str(int(chrono[0])) + ':' + '%02d' % (
                float(chrono[1]),) + ':%02d\t\t(min:sec:dec)\n' % (
                      float(chrono[1]) * 100 % 100,))

            i += 1

def plotHET(maglie, resp, pert):

    with open(model+'/HET.json') as fp:
        res_sens = json.load(fp)

    hetero = np.array([res_sens[str(m)]['hetero'] for m in maglie])
    homo = np.array([res_sens[str(m)]['homo'] for m in maglie])

    x = maglie

    x2 = np.array([nucData.giorni/m for m in maglie])

    fig, ax1 = plt.subplots()

    #ax1.set(xlabel='steps number', ylabel='error [%]', title = 'DPT convergence for '+resp+', perturbation: '+pert+' MT='+MT)
    ax1.set(xlabel='average step length [days] ', ylabel='DPT error [%]', title = 'DPT convergence for '+resp+', perturbation: '+pert)

    ax1.grid()

    ax1.plot(x2, hetero, 'b', marker='o', label = 'hetero')
    ax1.plot(x2, homo, 'r', marker='o', label = 'homo')

    ax1.legend(loc='upper left')

    fig.savefig(model + '/hetStudy.png')

def evoPtero(giorni, resp, pert):

    start = datetime.now()

    resDict(giorni, [], model+'/EVO')

    config.set('pterodax', 'PERT_NUC', pert)
    config.set('pterodax', 'RESP_NUC', '942390')
    config.set('pterodax', 'RESPONSE', resp)
    config.set('pterodax', 'ND', 'False')
    config.set('pterodax', 'direct', 'True')

    i   = 0
    tot = len(giorni)
    EOL = nucData.time[-1]


    for m in giorni:

        print('\nCalculation ' + str(i + 1) + '/' + str(tot))

        steps = round(EOL*m/25)

        config.set('sibyl', 'PASSI', str(steps))
        config.set('sibyl', 'daystop', str(m))

        with open(r"configfile.ini", 'w') as configfile:
            config.write(configfile)

        #os.system('python step.py > /dev/null')
        os.system('python step.py')

        end = datetime.now()
        diff = end - start
        chrono = divmod(diff.total_seconds(), 60)

        print('Done!')
        print('Running time: ' + str(int(chrono[0])) + ':' + '%02d' % (
            float(chrono[1]),) + ':%02d\t\t(min:sec:dec)\n' % (
                  float(chrono[1]) * 100 % 100,))

        i += 1

def plotEVO(giorni, resp, pert, EOL):

    with open(model+'/EVO.json') as fp:
        res_sens = json.load(fp)

    ptero = np.array([res_sens[str(m)]['ptero'] for m in giorni])
    sibyl = np.array([res_sens[str(m)]['sibyl'] for m in giorni])

    x = np.array(giorni) * EOL

    fig, ax1 = plt.subplots()

    #ax1.set(xlabel='steps number', ylabel='error [%]', title = 'DPT convergence for '+resp+', perturbation: '+pert+' MT='+MT)
    ax1.set(xlabel='Calculation length [days] ', ylabel='Keff sensitivity', title = resp+' sensitivity to '+pert+' concentration')

    ax1.grid()

    ax1.axhline(y=0, color='k')
    ax1.plot(x, ptero, 'b', marker='o', label = 'PTERODAx')
    ax1.plot(x, sibyl, 'r', marker='o', label = 'SIBYL')
    ax1.ticklabel_format(useOffset=False, style='plain')

    ax1.legend(loc='upper left')

    fig.savefig(model + '/evoPtero.png')



model = 'UO2/NEW'

config = CP.ConfigParser()

config.add_section('sibyl')

config.set('sibyl', 'model', model)
config.set('sibyl', 'energy', '2')
config.set('sibyl', 'hetSteps', 'True')
config.set('sibyl', 'fpSwitch', 'False')
config.set('sibyl', 'hetSwitch', 'False')

config.add_section('pterodax')

config.set('pterodax', 'ND', 'True')
config.set('pterodax', 'MT', '18')
config.set('pterodax', 'pert', '1.01')
config.set('pterodax', 'resetK', 'False')
config.set('pterodax', 'direct', 'False')
config.set('pterodax', 'sens_formula', 'False')


def main(run):

    if run == 'spesa':

        config.set('pterodax', 'run', run)
        RESP = ['keff', '942390', '922350']
        PERT = ['922350', '922380', '942390']
        MT = ['18', '102', '452']

        spesaSENS(RESP, PERT, MT)

    elif run == 'hetero':

        config.set('pterodax', 'run', run)
        maglie = [25, 50, 75, 100, 125, 150, 200, 250]
        hetStudy(maglie, 'keff', '922380')
        plotHET(maglie, 'keff', '922380')

    elif run == 'evoPtero':

        EOL = nucData.time[-1]
        config.set('pterodax', 'run', run)
        giorni = [0.2, 0.4, 0.6, 0.8, 1.0]
        evoPtero(giorni, 'keff', '942390')
        plotEVO(giorni, 'keff', '942390', EOL)

main('evoPtero')