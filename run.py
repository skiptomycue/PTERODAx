import os
import configparser as CP
import sys

config = CP.ConfigParser()

config.add_section('sibyl')

config.set('sibyl', 'model', 'UO2')
config.set('sibyl', 'energy', '2')
config.set('sibyl', 'PASSI', '3')
config.set('sibyl', 'fpSwitch', 'False')
config.set('sibyl', 'hetSwitch', 'False')

config.add_section('pterodax')

config.set('pterodax', 'ND', 'True')
config.set('pterodax', 'pert', '1.01')
config.set('pterodax', 'resetK', 'False')
config.set('pterodax', 'sens_formula', 'False')

RESP = ['keff', '942390']
PERT = ['922350', '922380']
MT   = ['18', '102', '452']

i   = 0
tot = len(RESP)*len(PERT)*len(MT)

for r in RESP:

    for p in PERT:

        for m in MT:

            print('\nCalculation '+str(i)+'/'+str(tot))

            resp = 'nuclide'

            if r == 'keff':

                resp = 'keff'
                r    = '942390'

            config.set('pterodax', 'PERT_NUC', p)
            config.set('pterodax', 'RESP_NUC', r)
            config.set('pterodax', 'RESPONSE', resp)
            config.set('pterodax', 'MT', m)

            with open(r"configfile.ini", 'w') as configfile:
                config.write(configfile)

            os.system('python step.py > /dev/null')

            print('Done!')

            i += 1