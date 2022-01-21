### MODEL INPUT ###

model        = 'UO2'                                                  # INPUT MODEL
energy       =  44                                                    # INPUT ENERGY GROUPS
PASSI        =  5                                                    # INPUT STEP NUMBER
fpSwitch     =  0                                                     # SWITCH TO FULL NUCLIDE CHART
hetSwitch    =  0                                                     # SWITCH TO HETEROGENEOUS CORRECTION FOR FUEL AND NICHEL

### CALCULATION INPUT

PERT         = ['922350', '922380']#, '280580', '50100']#, '10020']                                    # INPUT PERTURBATION NUCLIDE
RESP_NUC     =  '942390'                                              # OUTPUT RESPONSE NUCLIDE
RESPONSE     =  'keff'                                                # OUTPUT NUCLIDE, KEFF OR NONE
ND           =  True                                                  # SWITCH ND PERTURBATION
MT           =  '452'                                                 # INPUT PERTURBATION XS
pert         =  1.01                                                  # INPUT PERTURBATION %
resetK       =  0                                                     # SWITCH K-RESET
sens_formula =  False

