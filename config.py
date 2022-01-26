import configparser as CP
import json

config = CP.ConfigParser()
config.read("configfile.ini")
sibyl  = config['sibyl']
ptero  = config['pterodax']

def bull(stri):

    if stri == 'False' or stri == 'None' or stri == '0' or stri == 0:

        return False

    else:

        return True

### MODEL INPUT ###

model        =      sibyl['model']                                     # INPUT MODEL
energy       =      sibyl['energy']                                    # INPUT ENERGY GROUPS
PASSI        =  int(sibyl['passi'])                                    # INPUT STEP NUMBER
hetSteps     = bull(sibyl['hetSteps'])
fpSwitch     = bull(sibyl['fpswitch'])                                 # SWITCH TO FULL NUCLIDE CHART
hetSwitch    = bull(sibyl['hetswitch'])                                # SWITCH TO HETEROGENEOUS CORRECTION FOR FUEL AND NICHEL

### SENSITIVITY INPUT ###

PERT_NUC     =       [ptero['pert_nuc']]                               # INPUT PERTURBATION NUCLIDE
RESP_NUC     =        ptero['resp_nuc']                                # OUTPUT RESPONSE NUCLIDE
RESPONSE     =        ptero['response']                                # OUTPUT NUCLIDE, KEFF OR NONE
ND           =   bull(ptero['nd'])                                     # SWITCH ND PERTURBATION
MT           =        ptero['mt']                                      # INPUT PERTURBATION XS
pert         =  float(ptero['pert'])                                   # INPUT PERTURBATION %
resetK       =   bull(ptero['resetk'])                                 # SWITCH K-RESET
direct       =   bull(ptero['direct'])                                 # SWITCH K-RESET
sens_formula =   bull(ptero['sens_formula'])

