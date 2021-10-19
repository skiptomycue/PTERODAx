volumes  ={}
regions  ={}
unis     ={}
detectors={}
materials={}
zais     ={}
power    ={}
inp      ={}

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
sama = [611491, 621491, 621501, 621511, 631491, 631501, 631511, 641521]


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

ZAI1 = ['531350','541350', '601490', '611490', '621490', '922340', '922350', '922380', '922390', '932390', '942390']
ZAI2 = ['400900', '80160']
ZAI3 = ['10010', '80160']
ZAI  = [ZAI1,ZAI2,ZAI3]

POW = 1.55E+2
INP = 'REP'

buildDict(model1, VOL, REG, UNI, DET, ZAI, MAT, POW, INP)


model2 = 'LEU'

volClad = 1.33388E+04
volMeat = 1.1821E+4
volB    = 1.066E+2
volNi   = 6.541E+3
volAl   = 3.94478E+03
volHwC  = 1.311E+5
volHwF  = 3.625E+4
volHwR  = 1.10417E+06

volCent = volNi + volHwC + volAl
volFuel = volMeat + volHwF + volClad
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

