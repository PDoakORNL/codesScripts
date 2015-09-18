#  Authors: Jonathan Anchell, Ram Balachandran.
#  Center for Nanophase Materials Sciences, Oak Ridge National Laboratory, Oak Ridge, TN.

import sys
import math
from operator import itemgetter
import numpy as np
import pymatgen as mg
import copy
from collections import namedtuple
from enum import Enum

class densityOfStates:
    #create a constructor with the structure information (pymatgen structure).
    # Everything needs to be reimplemented by importing phonopy into python
    def total_dos(self, tdosFileName, energyRange):
        temp_list = []
        with open(tdosFileName, 'r') as f:
            next(f)
            for line in f:
                val = [float(i) for i in line.split()]
                temp_list.append(val)

        total_dos_list = []
        for cnt, val in enumerate(temp_list):
            if (energyRange[0]<=val[0]<=energyRange[1]):
                total_dos_list.append(val)
        return np.asarray(total_dos_list)

    # pass which types of atom we are interested. The index in turn can be obtained from the structure created in constructor
    def partial_dos(self, pdosFileName, atomIndexList, energyRange):
        temp_list = []
        energy_list = []
        with open(pdosFileName, 'r') as f:
            next(f)
            for line in f:
                lineVal = line.split()
                energyVal = float(lineVal[0])
                energy_list.append(energyVal)
                pdosVal = [float(i) for i in lineVal[1:]]
                temp_list.append(pdosVal)

        partial_dos_list=[]
        for cnt, val in enumerate(energy_list):
            pdosVal = 0
            if (energyRange[0]<= val <= energyRange[1]):
                for index in atomIndexList:
                    try:
                        pdosVal += temp_list[cnt][index]
                    except IndexError:
                        print 'Index does not exist'
                partial_dos_list.append([val,pdosVal])
        return np.asarray(partial_dos_list)
        
