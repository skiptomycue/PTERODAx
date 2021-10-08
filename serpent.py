volumes  ={}
regions  ={}
unis     ={}
detectors={}
materials={}
zais     ={}

def buildDict(model, VOL,REG,UNI,DET,ZAI, MAT):

    volumes[model]   = VOL
    regions[model]   = REG
    unis[model]      = UNI
    detectors[model] = DET
    zais[model]      = ZAI
    materials[model] = MAT

REACTIONS = ['101', '102', '103', '107', '16', '17', '18' ]
sama =[]
#sama = [621481, 621501, 621511, 631491, 641521, 631511, 611491, 631491, 621491, 621501, 631501]


model1 = 'UO2-pin'

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

buildDict(model1, VOL, REG, UNI, DET, ZAI, MAT)