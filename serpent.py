import json

volumes  ={}
regions  ={}
unis     ={}
detectors={}
materials={}
zais     ={}
power    ={}
inp      ={}

def sensDict(RESP, PERT):

    res_sens={}

    for a in RESP:
        res_sens[a]= {}

        for b in PERT:
            res_sens[a][b] = {}


    with open('COVX/SENS.json', 'w') as fp:
        json.dump(res_sens, fp)
def buildDict(model, VOL,REG,UNI,DET,ZAI, MAT, POW, INP):

    volumes[model]   = VOL
    regions[model]   = REG
    unis[model]      = UNI
    detectors[model] = DET
    zais[model]      = ZAI
    materials[model] = MAT
    power[model]     = POW
    inp[model]       = INP

REACTIONS = ['101', '102', '103', '107', '16', '17', '18' ]
sama =[]
#sama1 = [601470, 601471, 611470, 611471, 611480, 611481, 611491, 621480, 621491, 621501, 621511, 631491, 631501, 631511, 641521]
#sama2 = ['621480', '621500', '621510', '631490', '641520', '621481', '621501', '621511', '631491', '641521', '631511', '611490', '631490', '641530', '621500', '631500', '611491', '631491', '641531', '621491', '621501', '631501']
#sama = list(set(sama1 + [int(a) for a in sama2]))

model1 = 'UO2'

volZr = 2.16817E-01
volHe = 2.84600E-02
volClad = volZr + volHe
volFuel = 6.92643E-01
volCool = 1.14346E+00

VOL=[[volFuel], [volClad], [volCool]]
REG=['Fuel', 'Cladding', 'Coolant']
UNI=['20','2','50']
DET=['fuel', 'clad', 'cool']
MAT=[['fuel'], ['clad'], ['water']]

ZAI1 = ['531350','541350', '601490', '611490', '621490', '922340', '922350', '922380', '922390', '932390', '942390', '942400']
ZAI2 = ['400900', '80160']
ZAI3 = ['10010', '80160']
ZAI  = [ZAI1,ZAI2,ZAI3]

POW = 1.55E+2
INP = 'REP'

buildDict(model1, VOL, REG, UNI, DET, ZAI, MAT, POW, INP)


model2 = 'LEU'

#volMeat = 1.1821E+4

volClad = 1.33388E+04
volMeat = 1.214E+4
volB    = 1.066E+2
volNi   = 6.541E+3
volAl   = 3.94478E+03
volHwC  = 1.311E+5
volHwF  = 3.625E+4
volHwR  = 1.10417E+06
volAG3  = 3.94032E+03

volCent = volNi + volHwC + volAl
volFuel = volMeat + volHwF + volClad + volAG3
volRefl = volHwR + volB

#VOL=[[volCent], [volMeat], [volRefl]]
VOL=[[volNi , volHwC ], [volMeat , volHwF , volClad], [volHwR, volB]]
REG=[ 'Central', 'Fuel', 'Reflector']
UNI=['10','2','11']
DET=['cent', 'fuel', 'refl']
MAT=[['Ni','hwat_ring'], ['Fuel', 'hwat_in', 'Clad'], ['hwat_ref','boro']]
#MAT=[['hwat_ring'], ['Fuel'], ['hwat_ref']]

ZAI1 = ['280580', '280600', '10010', '10020', '80160']
ZAI2 = ['531350','541350', '601490', '611490', '621490', '922340', '922350', '922380', '922390', '932390', '942390', '130270', '10010', '10020', '80160']
ZAI3 = ['50100', '10010', '10020', '80160']
ZAI  = [ZAI1,ZAI2,ZAI3]

POW = 58.3E+6
INP = 'RHF'

buildDict(model2, VOL, REG, UNI, DET, ZAI, MAT, POW, INP)


model3 = 'UO2/NEW'

volZr = 2.16817E-01
volHe = 2.84600E-02
volClad = volZr + volHe
volFuel = 6.92643E-01
volCool = 1.14346E+00

VOL=[[volFuel], [volClad], [volCool]]
REG=['Fuel', 'Cladding', 'Coolant']
UNI=['1000','1001','1002']
DET=['REG1', 'REG2', 'REG3']
MAT=[['Fuel'], ['Clad'], ['Water']]

ZAI1 = ['531350','541350', '601490', '611490', '621490', '922340', '922350', '922380', '922390', '932390', '942390', '942400']
ZAI2 = ['400900', '80160']
ZAI3 = ['10010', '80160']
ZAI  = [ZAI1,ZAI2,ZAI3]

POW = 1.55E+2
INP = 'INP'

buildDict(model3, VOL, REG, UNI, DET, ZAI, MAT, POW, INP)


model4 = 'LEU/NEW'

volClad = 1.29039E+04
volMeat = 1.14367E+04
volB    = 1.76479E+02
volNi   = 6.51249E+03
volAl   = 3.94478E+03
volHwC  = 7.50280E+04
volHwF  = 3.46045E+04
volHwR  = 1.10711E+06

volCent = volNi + volHwC + volAl
volFuel = volMeat + volHwF + volClad
volRefl = volHwR + volB

#VOL=[[volCent], [volMeat], [volRefl]]
VOL=[[volNi , volHwC ], [volMeat , volHwF , volClad], [volHwR, volB]]
REG=[ 'Central', 'Fuel', 'Reflector']
UNI=['1000','1001','1002']
DET=['REG1', 'REG2', 'REG3']
MAT=[['Ni','hwat_ring'], ['Fuel', 'hwat_in', 'Clad'], ['hwat_ref','boro']]
#MAT=[['hwat_ring'], ['Fuel'], ['hwat_ref']]

ZAI1 = ['280580', '280600', '10010', '10020', '80160']
ZAI2 = ['531350','541350', '601490', '611490', '621490', '922340', '922350', '922380', '922390', '932390', '942390', '942400', '130270', '10010', '10020', '80160']
ZAI3 = ['50100', '10010', '10020', '80160']
ZAI  = [ZAI1,ZAI2,ZAI3]

POW = 58.3E+6
INP = 'INP'

buildDict(model4, VOL, REG, UNI, DET, ZAI, MAT, POW, INP)



model5 = 'LEU/OLD'

#volMeat = 1.1821E+4

volClad = 1.33388E+04
volMeat = 1.214E+4
volB    = 1.066E+2
volNi   = 6.541E+3
volAl   = 3.94478E+03
volHwC  = 1.311E+5
volHwF  = 3.625E+4
volHwR  = 1.10417E+06

volCent = volNi + volHwC + volAl
volFuel = volMeat + volHwF + volClad #+ volAG3
volRefl = volHwR + volB

#VOL=[[volCent], [volMeat], [volRefl]]
VOL=[[volNi , volHwC ], [volMeat , volHwF , volClad], [volHwR, volB]]
REG=[ 'Central', 'Fuel', 'Reflector']
UNI=['10','2','11']
DET=['cent', 'fuel', 'refl']
MAT=[['Ni','hwat_ring'], ['Fuel', 'hwat_in', 'Clad'], ['hwat_ref','boro']]
#MAT=[['hwat_ring'], ['Fuel'], ['hwat_ref']]

ZAI1 = ['280580', '280600', '10010', '10020', '80160']
ZAI2 = ['531350','541350', '601490', '611490', '621490', '922340', '922350', '922380', '922390', '932390', '942390', '942400', '130270', '10010', '10020', '80160']
ZAI3 = ['50100', '10010', '10020', '80160']
ZAI  = [ZAI1,ZAI2,ZAI3]

POW = 58.3E+6
INP = 'RHF'

buildDict(model5, VOL, REG, UNI, DET, ZAI, MAT, POW, INP)


model6 = 'HEU'

volClad = 1.61035E+04
volMeat = 7.68092E+03
volB    = 5.74833E+02
volNi   = 6.51249E+03
volAl   = 3.92740E+03
volHwC  = 8.37647E+04
volHwF  = 3.45781E+04
volHwR  = 1.10347E+06

volCent = volNi + volHwC + volAl + volB
volFuel = volMeat + volHwF + volClad
volRefl = volHwR

#VOL=[[volCent], [volMeat], [volRefl]]
VOL=[[volNi , volHwC ], [volMeat , volHwF , volClad, volB], [volHwR]]
REG=[ 'Central', 'Fuel', 'Reflector']
UNI=['1000','1001','1002']
DET=['REG1', 'REG2', 'REG3']
MAT=[['Ni','hwat_ring'], ['Fuel', 'hwat_in', 'Clad', 'boro'], ['hwat_ref']]
#MAT=[['hwat_ring'], ['Fuel'], ['hwat_ref']]

ZAI1 = ['280580', '280600', '10010', '10020', '80160']
ZAI2 = ['50100', '531350','541350', '601490', '611490', '621490', '922340', '922350', '922380', '130270', '10010', '10020', '80160']
ZAI3 = ['10010', '10020', '80160']
ZAI  = [ZAI1,ZAI2,ZAI3]

POW = 58.3E+6
INP = 'INP'

buildDict(model6, VOL, REG, UNI, DET, ZAI, MAT, POW, INP)