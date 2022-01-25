import os
import configparser as CP
import json
from datetime import datetime


def sensDict(RESP, PERT):

    res_sens={}

    for a in RESP:
        res_sens[a]= {}

        for b in PERT:
            res_sens[a][b] = {}


    with open('COVX/SENS.json', 'w') as fp:
        json.dump(res_sens, fp)

config = CP.ConfigParser()

config.add_section('sibyl')

config.set('sibyl', 'model', 'UO2/NEW')
config.set('sibyl', 'energy', '44')
config.set('sibyl', 'PASSI', '100')
config.set('sibyl', 'fpSwitch', 'False')
config.set('sibyl', 'hetSwitch', 'False')

config.add_section('pterodax')

config.set('pterodax', 'ND', 'True')
config.set('pterodax', 'pert', '1.01')
config.set('pterodax', 'resetK', 'False')
config.set('pterodax', 'sens_formula', 'False')

RESP = ['keff','942390']
PERT = ['922350', '922380', '942390']
MT   = ['18', '102', '452']
sensDict(RESP, PERT)

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