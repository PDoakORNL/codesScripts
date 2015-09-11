#  Sample Runs of analyzeStructure3D module.
import pymatgen as mg
from analyzeStructure3D import analyzePerovskite3D
from analyzeStructure3D import makePerovskite3D


#  Make poscar structures.
poscar_01 = mg.Structure.from_file('POSCAR_01')
poscar_03 = mg.Structure.from_file('POSCAR_03')


#  Run analyzePerovskite3D class.
#  analyzePerovskite3D.find_angle(self, poscar, ref_atom, angles, radius, verbose=False)
analyzePerovskite3D().find_angle(poscar_03, ['O', 102], [['O', 'Al', 'O'], ['O', 'Zr', 'O']], 5, verbose=False)

#  analyzePerovskite3D.find_displacement(self, poscar_start, poscar_relax, ref_atom, neighbor_types, radius, verbose=False)
analyzePerovskite3D().find_displacement(poscar_01, poscar_03, ['O', 102], ['Zr', 'Al'], 5, verbose=False)

#  analyzePerovskite3D().find_distance(self, poscar, ref_atom, neighbor_types, radius=False, number_of_neighbors=False, verbose=False)
analyzePerovskite3D().find_distance(poscar_03, ['O', 102], ['Zr', 'Al'], radius=5, number_of_neighbors=False, verbose=False)


#  Run makePerovskite3D class.
#  makePerovskite3D().create_pure_structure(self, oxidation_states, basic_atoms, noCells, lattice_system = 'cubic', lattice_lengths = None, lattice_angles = None, supercell=False):
structure, keyVal = makePerovskite3D().create_pure_structure([2, 4, -2], ['Ba', 'Zr', 'O'], [2, 2, 2], lattice_system='cubic', lattice_lengths=None, lattice_angles=None, supercell=False)

#  makePerovskite3D().create_substitutional_dopant_structure(self, origStruc, existingAtomType, atomIndices, dopantAtomType)
makePerovskite3D().create_substitutional_dopant_structure(structure, keyVal, [3, 4], 'Y')

#  makePerovskite3D().def create_vacancy_defect_structure(self, origStruc, keyVal, vacancyIndices):
makePerovskite3D().create_vacancy_defect_structure(structure, keyVal, [2, 3])

#  makePerovskite3D().create_interstitial_dopant_structure(self, origStruc, keyVal, interstitialAtoms, interstitialIndices, distances, directions):
makePerovskite3D().create_interstitial_dopant_structure(structure, keyVal, ['H'], [2, 3], [.5, .2], [[0, 1, 0], [1, 0, 0]])

#  makePerovskite3D.move_interstitial_dopant(self, origStruc, keyVal, atoms, indicesToMove, distancesFromIndices, directionsFromIndices)
makePerovskite3D().move_interstitial_dopant(structure, keyVal, ['H'], [3, 4], [.6, .7], [[0, 0, 1], [1, 2, 3]])

