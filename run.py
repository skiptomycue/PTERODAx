import os
import configparser as CP
import json

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

config.set('sibyl', 'model', 'UO2')
config.set('sibyl', 'energy', '44')
config.set('sibyl', 'PASSI', '10')
config.set('sibyl', 'fpSwitch', 'False')
config.set('sibyl', 'hetSwitch', 'False')

config.add_section('pterodax')

config.set('pterodax', 'ND', 'True')
config.set('pterodax', 'pert', '1.01')
config.set('pterodax', 'resetK', 'False')
config.set('pterodax', 'sens_formula', 'False')

RESP = ['keff']
PERT = ['922350', '922380']
MT   = ['18', '102', '452']
sensDict(RESP, PERT)

i   = 0
tot = len(RESP)*len(PERT)*len(MT)

for r in RESP:

    for p in PERT:

        for m in MT:

            print('\nCalculation '+str(i)+'/'+str(tot))

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

            print('Done!')

            i += 1