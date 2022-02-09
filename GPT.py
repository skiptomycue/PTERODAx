import serpentTools as ST
from periodictable import elements
import plotly.graph_objects as go
import numpy as np
import os

day = '0'

sens  = ST.readDataFile('k-reset/LEU/'+day+'/branch/INP_sens0.m')
sensH = ST.readDataFile('k-reset/HEU/'+day+'/branch/INP_sens0.m')

SSK =  sens.energyIntegratedSens['keff']
SSH = sensH.energyIntegratedSens['keff']

def getName(zai):

    z=zai[:-4]
    a=zai[-4:-1]

    if int(a) in elements[int(z)].isotopes:

        Name=elements[int(z)][int(a)].name

    elif zai == 60000:

        Name = 'carbon'

    else:

        Name = 'sarcazzo'

    Name += '-'+str(int(a))

    Name = Name[0].capitalize() + Name[1:]

    return Name

#ssk = SSK[0, idZai, idPert, :, 0]

print(sens.energyIntegratedSens['keff'].shape)
print(sens.materials.keys())
print(sens.zais.keys())
print(sens.perts.keys())

print('\n\nEnergy integrated sensitivities: \n')

m = 0

MATS = ['Fuel       ', 'Reflector  ', 'Poison     ', 'Coolant    ', 'Control Rod']

print('Material\tNuclide    \tPerturbation\t\tLEU      \tHEU     \tDifference\n')

header = ['Material', 'Nuclide', 'Perturbation', 'LEU', 'HEU', 'Difference', 'Uncertainty']

values = []

for M in MATS:

    z = 0

    for Z in sens.zais.keys():

        p = 0

        for P in sens.perts.keys():

            ssk  = int(SSK[m, z, p, 0] >= 0)*"+" + "{:.2e}".format(SSK[m, z, p, 0])
            ssH  = int(SSH[m, z, p, 0] >= 0)*"+" + "{:.2e}".format(SSH[m, z, p, 0])

            if abs(SSK[m, z, p, 0]) > 1E-4 or abs(SSH[m, z, p, 0]) > 1E-4 :

                name = getName(str(Z))

                d    = round((SSK[m, z, p, 0] - SSH[m, z, p, 0]) / SSK[m, z, p, 0] * 100)
                diff = int(d > 0)*"+" + str(d)

                if abs(d) >= 100:

                    diff = '$>$100'

                values.append([M, name, P, str(ssk), str(ssH), str(diff)+' \%', str(int(SSK[m, z, p, 1]*100))])
                print(M + '\t' + name + '\t' + P + '\t\t' + str(ssk) + '\t' + str(ssH) + '\t' + str(diff)  +'%')

            p += 1

            if p == 1:

                break

        #print()
        z +=1

    print()
    m +=1


vals = np.array(values).transpose().tolist()



def latexTable(caption, header, vals):

    cols = len(header)

    print('\\begin{table}[H]\n\\centering\n\\caption{'+ caption +'}')
    print('\\begin{tabular}{|'+ 'c|'*cols +'} \hline')

    head = ['\\textbf{'+h+'}\t&' for h in header]
    joined = ' '.join(x for x in head )
    print(joined[:-1] + " \\\ \hline")

    for v in vals:

        head = [h + '  &' for h in v]
        joined = ' '.join(x for x in head)
        print(joined[:-1] + "\%  \\\ ")

    print('\\hline\n\\end{tabular}\n\\label{table}\n\\end{table}')

latexTable('Energy integrated sensitivities at day '+day, header, values)

