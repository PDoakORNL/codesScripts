import phonopyPostProcess  
import matplotlib.pyplot as plt
import itertools
import os, sys
import re
import numpy as np

dos=phonopyPostProcess.densityOfStates()

folderPath = {"WithH": ['/Users/kw1/Documents/ResearchProjects/IonicConductivity/Data/Perovskites/Analysis/BZO_Phonons/BZOH_ISIF2',
                          '/Users/kw1/Documents/ResearchProjects/IonicConductivity/Data/Perovskites/Analysis/BZO_Phonons/BZYOH_3X3X3_NOCluster_ISIF2_Oindex102',
                          '/Users/kw1/Documents/ResearchProjects/IonicConductivity/Data/Perovskites/Analysis/BZO_Phonons/BZYOH_3X3X3_Cluster_ISIF2_Oindex66']}

inFileName = 'partial_dos.dat'
outPath = '/Users/kw1/Documents/ResearchProjects/IonicConductivity/Data/Perovskites/Analysis/BZO_Phonons/'
pltFileName = 'partial-dos-H.eps'
atomIndexList =[135,]
energyRange = [0.0, 40.0]
colors = itertools.cycle(["k", "r", "g"])

for folder in folderPath["WithH"]:
    strVal = '-'.join(re.split('/|_',folder)[11:13])
    file = os.path.join(folder, inFileName)
    dosList = dos.partial_dos(file, atomIndexList, energyRange)
    plt.plot(dosList[:,0],dosList[:,1],color=next(colors),label=strVal)
    
plt.title('P-DoS WithH')
plt.legend(loc = 'lower right', fontsize='x-small')
pltFile = os.path.join(outPath, pltFileName)
plt.savefig(pltFile,format='eps',dpi=1000)
plt.close()
