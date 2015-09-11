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


class analyzePerovskite3D:
    def get_distance(self, cart_start, cart_end):
        """Given two sets of coordinates, calculates distance."

        """

        dist_x = cart_end[0] - cart_start[0]
        dist_y = cart_end[1] - cart_start[1]
        dist_z = cart_end[2] - cart_start[2]

        distance = (dist_x ** 2 + dist_y ** 2 + dist_z ** 2) ** .5

        return distance

    def get_vector(self, atom1_coords, atom2_coords):
        """Returns normalized vector that points from atom1 to atom2.

        """

        coords = ['x', 'y', 'z']

        #  Create unnormalized vector.
        unnormalized_vector = {}
        count = 0
        for coord in coords:
            unnormalized_vector[coord] = atom1_coords[count] - atom2_coords[count]
            count += 1

        # Create normalized vector.
        vector_length = (unnormalized_vector['x'] ** 2 +
                         unnormalized_vector['y'] ** 2 + unnormalized_vector['z'] ** 2) ** .5

        normalized_vector = {}
        for coord in coords:
            normalized_vector[coord] = unnormalized_vector[coord] / vector_length

        return normalized_vector

    def get_dot_product(self, vector1, vector2):
        """Returns dot product of two vectors.

        """

        dot_product = vector1['x'] * vector2['x'] + vector1['y'] * vector2['y'] + vector1['z'] * vector2['z']

        return dot_product

    def get_angle(self, dot_product):
        """Returns angle given by dot product.

        """

        angle_degrees = math.degrees(math.acos(dot_product))

        return angle_degrees

    def find_angle(self, poscar, ref_atom, angles, radius, verbose=False):
        """Find all angles of type angles within specified radius.

        """

        #  Get left, middle and right atom set.
        left_atom_set = []
        middle_atom_set = []
        right_atom_set = []
        for angle_list in angles:
            left_atom = angle_list[0]
            middle_atom = angle_list[1]
            right_atom = angle_list[2]

            left_atom_set.append(left_atom)
            middle_atom_set.append(middle_atom)
            right_atom_set.append(right_atom)

        # Create global list of sites.
        sites = poscar.sites

        #  Get reference atom index.
        ref_index = int(ref_atom[1])

        #  Get reference atom species from index.
        species_list = poscar.species
        ref_species_from_index = str(species_list[ref_index])

        #  Get reference atom species from input.
        ref_species_from_input = ref_atom[0]

        #  Verify correctness of reference atom input.
        if ref_species_from_index != ref_species_from_input:
            print "Reference index does not match reference species."
            sys.exit(1)

        # Verify 'ref_atom' species is the first atom in 'angle' input.
        for entry in left_atom_set:
            if ref_species_from_index != entry:
                print "Reference atom does not match first entry in angle list."
                sys.exit(1)

        # Verify radius is a float and not a string.
        radius = float(radius)

        # Get reference site and reference coordinates.
        ref_site = sites[ref_index]
        ref_coords = ref_site.coords

        #  Get neighbors to reference site within specified radius.
        ref_neighbors = poscar.get_neighbors(ref_site, radius, include_index='true')
        ref_neighbors.sort(key=itemgetter(1))

        #  Instantiate angle list.
        angle_results = []

        #  Loop through neighbors to reference site.
        for ref_neighbor1 in ref_neighbors:
            ref_neighbor1_index = ref_neighbor1[2]
            ref_neighbor1_site = ref_neighbor1[0]
            ref_neighbor1_species = str(ref_neighbor1_site.specie)
            ref_neighbor1_coords = ref_neighbor1_site.coords

            #  Verify neighbor1_species matches input specification.
            if ref_neighbor1_species in right_atom_set:
                distance_list = []

                #  Start a second loop through neighbors to reference site.
                for ref_neighbor2 in ref_neighbors:
                    ref_neighbor2_index = ref_neighbor2[2]
                    ref_neighbor2_site = ref_neighbor2[0]
                    ref_neighbor2_species = str(ref_neighbor2_site.specie)
                    ref_neighbor2_coords = ref_neighbor2_site.coords

                    #  Verify neighbor2_species matches input specification.
                    if (ref_neighbor1 != ref_neighbor2) and (ref_neighbor2_species in middle_atom_set):
                        #  Find distance between ref_neighbor1 and ref_neighbor2.
                        refNeighbor1_refNeighbor2_dist = self.get_distance(ref_neighbor1_coords, ref_neighbor2_coords)

                        #  Find distance between ref atom and ref_neighbor2.
                        refAtom_refNeighbor2_dist = self.get_distance(ref_coords, ref_neighbor2_coords)

                        #  Sum distances.
                        total_distance = refNeighbor1_refNeighbor2_dist + refAtom_refNeighbor2_dist

                        #  Create distance list which holds: total distance, ref_neighbor2_index, ref_neighbor1_index,
                        #  ref_neighbor2_coords, ref_neighbor1_coords.
                        temp_list = [total_distance, ref_neighbor2_index, ref_neighbor1_index,
                                     ref_neighbor2_coords, ref_neighbor1_coords]
                        distance_list.append(temp_list)

                # Sort distance_list by total_distance.
                sorted_distance_list = sorted(distance_list)

                #  Find normalized vector between ref_neighbor1 and ref_neighbor2.
                refNeighbor1_refNeighbor2_vector = self.get_vector(sorted_distance_list[0][4],
                                                                   sorted_distance_list[0][3])

                #  Find normalized vector between ref atom and ref_neighbor2.
                refAtom_refNeighbor2_vector = self.get_vector(ref_coords, sorted_distance_list[0][3])

                #  Use vectors to find dot product and subsequent angle.
                dot_product = self.get_dot_product(refNeighbor1_refNeighbor2_vector, refAtom_refNeighbor2_vector)
                angle = self.get_angle(dot_product)

                #  Append to angle_results.
                angle_results.append(angle)

                #  Print results if specified.
                if verbose:
                    print ("%s - %s - %s: %s" % (ref_index, sorted_distance_list[0][1],
                                                 sorted_distance_list[0][2], angle))

        # Find number of angles.
        number_of_angles = len(angle_results)

        #  Find average angle.
        average_angle = np.mean(angle_results)

        #  Find angle standard deviation.
        stdev_angle = np.std(angle_results)

        #  Print results.
        if verbose:
            print("\nNumber of angles: %s" % number_of_angles)
            print("Average: %s" % average_angle)
            print("Standard Deviation: %s" % stdev_angle)

    def find_displacement(self, poscar_start, poscar_relax, ref_atom, neighbor_types, radius, verbose=False):
        """Find displacement of atoms in neighbor_types after relaxation.

        """

        #  Create lists of sites.
        sites_start = poscar_start.sites
        sites_relax = poscar_relax.sites

        #  Get reference atom index.
        ref_index = int(ref_atom[1])

        #  Get reference atom species from index.
        species_list = poscar_start.species
        ref_species_from_index = str(species_list[ref_index])

        #  Get reference atom species from input.
        ref_species_from_input = ref_atom[0]

        #  Verify correctness of reference atom input.
        if ref_species_from_index != ref_species_from_input:
            print "Reference index does not match reference species."
            sys.exit(1)

        # Verify neighbor_types is a list.  This is only necessary if neighbor_types is a single atom.
        if isinstance(neighbor_types, list) is False:
            print "neighbor_types must be a list."
            sys.exit(1)

        # Verify radius is a float and not a string.
        radius = float(radius)

        # Get reference sites and reference coordinates.
        ref_site_start = sites_start[ref_index]
        ref_site_relax = sites_relax[ref_index]

        ref_coords_start = ref_site_start.coords
        ref_coords_relax = ref_site_relax.coords

        #  Get neighbors to reference sites within specified radius.
        ref_neighbors_start = poscar_start.get_neighbors(ref_site_start, radius, include_index='true')
        ref_neighbors_start.sort(key=itemgetter(1))
        ref_neighbors_relax = poscar_relax.get_neighbors(ref_site_relax, radius, include_index='true')

        #  Instantiate results and displacement lists.
        results_list = []
        displacement_list = []

        #  Loop through neighbors to reference atom.
        for count in range(0, len(ref_neighbors_start)):
            species_start = str(ref_neighbors_start[count][0].specie)

            #  Verify species in neighbor_types.
            if species_start in neighbor_types:
                index = ref_neighbors_start[count][2]
                coordinates_start = ref_neighbors_start[count][0].coords

                for neighbor in ref_neighbors_relax:

                    #  Get coordinates corresponding to index in poscar_relax.
                    if index == neighbor[2]:
                        coordinates_relax = neighbor[0].coords

                # Get start and relaxed distances.
                distance_start = self.get_distance(ref_coords_start, coordinates_start)
                distance_relax = self.get_distance(ref_coords_relax, coordinates_relax)

                #  Get displacement.
                displacement = distance_relax - distance_start

                #  Append to displacement_list.
                displacement_list.append(displacement)

                #  Create results matrix.
                temp_list = [index, distance_start, distance_relax, displacement]
                results_list.append(temp_list)

        # Find average displacement.
        average_displacement = np.mean(displacement_list)

        #  Print results.
        if verbose:
            for entry in results_list:
                print ("Start distance %s - %s: %s" % (ref_index, entry[0], entry[1]))
                print ("Relax distance %s - %s: %s" % (ref_index, entry[0], entry[2]))
                print ("Displacement: %s\n" % entry[3])
            print ("Average displacement: %s" % average_displacement)

    def find_distance(self, poscar, ref_atom, neighbor_types, radius=False, number_of_neighbors=False, verbose=False):
        """If radius=number is specified, distances are found between ref_atom and all
         neighbor_type atoms within one radius of ref_atom.  If number_of_neighbors=number is specified,
         distances are found between reference atom and closest neighbor_type atoms up until
         number_of_neighbors is reached.

        """

        #  Create list of sites.
        sites = poscar.sites

        #  Get reference atom index.
        ref_index = int(ref_atom[1])

        #  Get reference atom species from index.
        species_list = poscar.species
        ref_species_from_index = str(species_list[ref_index])

        #  Get reference atom species from input.
        ref_species_from_input = ref_atom[0]

        #  Verify correctness of reference atom input.
        if ref_species_from_index != ref_species_from_input:
            print "Reference index does not match reference species."
            sys.exit(1)

        # Verify neighbor_types is a list.  This is only necessary if neighbor_types is a single atom.
        if isinstance(neighbor_types, list) is False:
            print "neighbor_types must be a list."
            sys.exit(1)

        # Verify radius and number of neighbors are a float and integer and not strings.
        if radius:
            radius = float(radius)
        if number_of_neighbors:
            number_of_neighbors = int(number_of_neighbors)

        # Get reference site and reference coordinates.
        ref_site = sites[ref_index]

        ref_coords = ref_site.coords

        #  Instantiate results and distance lists.
        results_list_radius = []
        results_list_neighbors = []
        distance_list_radius = []
        distance_list_neighbors = []

        #  If radius is specified find neighbors in neighbor_types.
        if radius:
            #  Get neighbors to reference site within specified radius.
            ref_neighbors_radius = poscar.get_neighbors(ref_site, radius, include_index='true')
            ref_neighbors_radius.sort(key=itemgetter(1))

            for neighbor in ref_neighbors_radius:
                neighbor_species = str(neighbor[0].specie)

                #  Check if neighbor_species is in neighbor_types.
                if neighbor_species in neighbor_types:
                    #  Get neighbor index and coordinates
                    neighbor_index = neighbor[2]
                    neighbor_coords = neighbor[0].coords

                    #  Get distance between neighbor and ref_atom.
                    distance = self.get_distance(ref_coords, neighbor_coords)

                    #  Append distance to distance_list_radius.
                    distance_list_radius.append(distance)

                    #  Create results matrix.
                    temp_list = [neighbor_index, distance]
                    results_list_radius.append(temp_list)

            # Find average distance.
            average_distance_radius = np.mean(distance_list_radius)

            #  Print results.
            if verbose:
                print "RADIUS RESULTS."
                for result in results_list_radius:
                    print("distance %s - %s: %s" % (ref_index, result[0], result[1]))
                print("\nNeighbors in radius: %s" % len(distance_list_radius))
                print ("Average distance: %s" % average_distance_radius)

        # If number of neighbors is specified, find distance /
        # from reference atom to specified number of neighbors in neighbor_types.
        if number_of_neighbors:
            used_index_list = []
            temp_radius = 1

            while len(distance_list_neighbors) < number_of_neighbors:
                #  Get neighbors to reference site within specified radius.
                ref_neighbors = poscar.get_neighbors(ref_site, temp_radius, include_index='true')
                #  Sort ref_neighbors by distance.
                ref_neighbors.sort(key=itemgetter(1))
                for neighbor in ref_neighbors:
                    neighbor_index = neighbor[2]
                    neighbor_species = str(neighbor[0].specie)
                    if (neighbor_species in neighbor_types) and (neighbor_index not in used_index_list) and len(
                            distance_list_neighbors) < number_of_neighbors:
                        #  Get neighbor coordinates.
                        neighbor_coords = neighbor[0].coords

                        #  Get distance between neighbor and ref_atom.
                        distance = self.get_distance(ref_coords, neighbor_coords)

                        #  Append distance to distance_list_neighbors.
                        distance_list_neighbors.append(distance)

                        #  Create results matrix.
                        temp_list = [neighbor_index, distance]
                        results_list_neighbors.append(temp_list)

                        #  Append to used index_list.
                        used_index_list.append(neighbor_index)

                # Increase radius.
                temp_radius += 1

            # Find average distance.
            average_distance_neighbors = np.mean(distance_list_neighbors)

            #  Print results.
            if verbose:
                if radius:
                    print ""
                print "NUMBER OF NEIGHBORS RESULTS."
                for result in results_list_neighbors:
                    print("distance %s - %s: %s" % (ref_index, result[0], result[1]))
                print ("\nAverage distance: %s" % average_distance_neighbors)


class atomSites(Enum):
    siteA = 0
    siteB = 1
    siteC = 2


class makePerovskite3D:
    """
    1. This class can create cubic perovskite structures ABO3 (of symmetric supercell sizes nxnxn) employing pymatgen.

    2. This class can create substitutional dopants (tested in siteB for 1 dopant).  The number, position of dopant
        hardcoded -- 1, max(siteB).

    3. This class can create vacancies (tested for siteC for 1 vacancy). The number, position of vacancy is
        hardcoded -- 1, max(siteC).

    4. This class can create simple interstitial atoms (tested for C-B bonds for one atom (pure), 2 atoms (doped)).
        The number, position of interstitials is hardcoded --
        a. Pure or Only Vacant defects --> [1, max(siteC)],
        b. Only Doped or Doped+Vacant defects --> [2, max(siteC-1), max(siteC-2)] where siteC-1 are sites whose
            the two nearest neighbors (for cubic perovskites) are only pure siteB, while siteC-2 are sites which have at
            least one of the two nearest neighbors is a dopant atom.

    5. The structures are saved as a dictionary in crystalStrucData. Each entry in crystalStrucData has a unique key
        associated with it whose format is given by structureKey.

    """

    def __init__(self):

        self.structureKey = namedtuple("structureKey",
                                       "latticeSystem basicAtoms dopantAtoms noDopants "
                                       "vacancyAtoms noVacancies noCells InterstitialAtoms noInterstitials")

    def return_site_index(self, noAtomSites, atomComposition, keyValue, siteName):
        """ Can return the atomic site index that needs to be doped (substituional)
        or needs to be removed (Vacancy). Currently can only return the maximum of site id. Needs to be generalized

        """

        noCells = keyValue[5]
        noAtomSites = noAtomSites
        noAtoms = [0] * noAtomSites
        atomsUnitCell = [0] * noAtomSites
        totalNoCells = 1
        siteIndex = atomSites[siteName]

        for i in noCells:
            totalNoCells *= i

        for key, val in atomComposition.iteritems():
            site = atomSites[key]
            atomsUnitCell[site.value] = val
            noAtoms[site.value] = val * totalNoCells

        atomSiteIndex = []
        indexStartVal = 0

        for i in range(0, noAtomSites):
            atomSiteIndex.append(range(indexStartVal, indexStartVal + noAtoms[i]))
            indexStartVal += noAtoms[i]

        return max(atomSiteIndex[siteIndex.value])

    def return_lattice_constant(self, noAtomSites, latticeSystem, keyValue):
        #  HARDCODING ALERT-- 1. keyValue 2.oxidationStates
        oxidationStates = [2, 4, -2]
        atoms = []
        ionicRadii = []

        for i in range(0, len(keyValue[1])):
            atoms.append(mg.Element(keyValue[1][i]))
            ionicRadii.append(atoms[i].ionic_radii[oxidationStates[i]])

        if latticeSystem == "cubic":
            latticeVal1 = 2 * (ionicRadii[0] + ionicRadii[1]) / math.pow(3, 1 / 2)
            latticeVal2 = 2 * (ionicRadii[0] + ionicRadii[2]) / math.pow(2, 1 / 2)
            latticeConst = max(latticeVal1, latticeVal2)
        else:
            print "ERROR: Only Cubic Perovskite Currently Implemented"
            latticeConst = float('NaN')

        return latticeConst

    #  Code adopted from https://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
    def rotation_matrix(self, axis, theta):
        """
        Return the rotation matrix associated with counterclockwise rotation about
        the given axis by theta radians.

        """

        axis = np.asarray(axis)
        theta = np.asarray(theta)
        axis /= math.sqrt(np.dot(axis, axis))
        a = math.cos(theta / 2)
        b, c, d = -axis * math.sin(theta / 2)
        aa, bb, cc, dd = a * a, b * b, c * c, d * d
        bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
        return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                         [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                         [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

    def create_pure_structure(self, oxidation_states, basic_atoms, noCells, lattice_system='cubic',
                              lattice_lengths=None, lattice_angles=None, supercell=False):
        """Creates pure structure given input parameters.

        """

        oxidationA = oxidation_states[0]
        oxidationB = oxidation_states[1]
        oxidationC = oxidation_states[2]

        dopantAtoms = []
        noDopants = []
        noVacancies = []
        noInterstitials = []
        vacancyAtoms = []
        interstitialAtoms = []

        siteA = basic_atoms[0]
        siteB = basic_atoms[1]
        siteC = basic_atoms[2]

        atomA = mg.Element(siteA)
        atomB = mg.Element(siteB)
        atomC = mg.Element(siteC)

        if not lattice_lengths:
            ionic_radiiA = atomA.ionic_radii[oxidationA]
            ionic_radiiB = atomB.ionic_radii[oxidationB]
            ionic_radiiC = atomC.ionic_radii[oxidationC]
            latticeVal1 = 2 * (ionic_radiiA + ionic_radiiB) / math.pow(3, 1 / 2)
            latticeVal2 = 2 * (ionic_radiiA + ionic_radiiC) / math.pow(2, 1 / 2)
            latticeConst = max(latticeVal1, latticeVal2)
            latticeLength = [latticeConst, latticeConst, latticeConst]
            latticeAngles = [90, 90, 90]  # WARNING: HARD-CODE ALERT
        else:
            latticeLength = lattice_lengths
            latticeAngles = lattice_angles

        lattice = mg.Lattice.from_lengths_and_angles(latticeLength, latticeAngles)
        structure = mg.Structure(lattice, [siteA, siteB, siteC, siteC, siteC],
                                 [[0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0.5, 0],
                                  [0.5, 0, 0.5], [0, 0.5, 0.5]])
        if supercell:
            structure.make_supercell(noCells)

        keyVal = self.structureKey(lattice_system, basic_atoms,
                                   dopantAtoms, noDopants, vacancyAtoms, noVacancies,
                                   noCells, interstitialAtoms, noInterstitials)

        return structure, keyVal

    def create_substitutional_dopant_structure(self, origStruc, keyVal, atomIndices, dopantAtomType):
        """Substitutes dopant at specified indices in structure.

        """

        #  Create structure.
        structure = origStruc

        #  Make substitution.
        for index in atomIndices:
            structure.replace(index, dopantAtomType)

        # Get initial key values.
        lattice_system = keyVal[0]
        vacancyAtoms = keyVal[4]
        noVacancies = keyVal[5]
        noCells = keyVal[6]
        interstitialAtoms = keyVal[7]
        noInterstitials = keyVal[8]

        #  Create dopantAtoms list.
        dopantAtoms = keyVal[2]
        noDopants = []

        if dopantAtomType not in dopantAtoms:
            dopantAtoms.append(dopantAtomType)

        # Find number of Dopants.
        for dopant in dopantAtoms:
            dopant_count = 0
            for site in structure:
                species = site.specie
                if str(dopant) == str(species):
                    dopant_count += 1
            noDopants.append(dopant_count)

        # Create basic atom list.
        basic_atoms = keyVal[1]

        for site in structure:
            species = str(site.specie)
            if species not in basic_atoms and species not in dopantAtoms:
                basic_atoms.append(species)

        # Create key.
        keyVal = self.structureKey(lattice_system, basic_atoms,
                                   dopantAtoms, noDopants, vacancyAtoms, noVacancies,
                                   noCells, interstitialAtoms, noInterstitials)

        return structure, keyVal

    def create_vacancy_defect_structure(self, origStruc, keyVal, vacancyIndices):
        """Creates vacancies at specified indices.

        """

        # Verify vacancyIndices is a list.
        if isinstance(vacancyIndices, list) is False:
            print "vacancyIndices must be a list."
            sys.exit(1)

        # Create structure.
        structure = origStruc

        #  Get initial key values.
        lattice_system = keyVal[0]
        basic_atoms = keyVal[1]
        dopantAtoms = keyVal[2]
        noDopants = keyVal[3]
        vacancyAtoms = keyVal[4]
        noVacancies = keyVal[5]
        noCells = keyVal[6]
        interstitialAtoms = keyVal[7]
        noInterstitials = keyVal[8]

        #  Create vacancyAtoms list.
        for index in vacancyIndices:
            site = structure[index]
            species = str(site.specie)
            if species not in vacancyAtoms:
                vacancyAtoms.append(species)

        # Create noVacancies list.
        preexisting_vacancy_index = 0
        for atom in vacancyAtoms:
            vacancy_count = 0
            index = 0
            for site in structure:
                species = str(site.specie)
                if str(atom) == str(species) and index in vacancyIndices:
                    vacancy_count += 1
                index += 1
            # Append vacancy_count to noVacancies.  If previous vacancies exist for atom,
            # add vacancy_count to previous vacancy_count.
            try:
                if noVacancies[preexisting_vacancy_index] > 0:
                    new_vacancy_count = noVacancies[preexisting_vacancy_index] + vacancy_count
                    noVacancies[preexisting_vacancy_index] = new_vacancy_count
            except:
                noVacancies.append(vacancy_count)
            preexisting_vacancy_index += 1

        # Remove vacancy sites.
        structure.remove_sites(vacancyIndices)

        #  Create key.
        keyVal = self.structureKey(lattice_system, basic_atoms,
                                   dopantAtoms, noDopants, vacancyAtoms, noVacancies,
                                   noCells, interstitialAtoms, noInterstitials)

        return structure, keyVal

    def create_interstitial_dopant_structure(self, origStruc, keyVal, atomsToInsert, neighboringIndices,
                                             distancesFromIndices, directionsFromIndices):
        """Creates interstitial dopants near neighboring indices.  Interstitial dopants are created at
        a specified distance and direction from the index.

        """

        # Verify atoms is a list.
        if isinstance(atomsToInsert, list) is False:
            print "atomsToInsert must be a list."
            sys.exit(1)

        # Verify interstitialIndices, distanceFromIndices and directionsFromIndices have the same length.
        if (len(neighboringIndices) != len(distancesFromIndices)) or (
                    len(neighboringIndices) != len(directionsFromIndices)) or (
                    len(distancesFromIndices) != len(directionsFromIndices)):
            print "neighboringIndices, distanceFromIndices and directionsFromIndices must be the same length."
            sys.exit(1)

        # Create structure.
        structure = origStruc

        #  Get initial key values.
        lattice_system = keyVal[0]
        basic_atoms = keyVal[1]
        dopantAtoms = keyVal[2]
        noDopants = keyVal[3]
        vacancyAtoms = keyVal[4]
        noVacancies = keyVal[5]
        noCells = keyVal[6]
        interstitialAtoms = keyVal[7]
        noInterstitials = keyVal[8]

        # Get interstitial atom species.  Currently Hardcoded to assume their is one only one interstitial atom species.
        interstitialAtom = atomsToInsert[0]

        #  Insert interstitials.
        for key, index in enumerate(neighboringIndices):
            coordinatesAtIndex = structure[index].coords
            distance = distancesFromIndices[key]
            directionVector = directionsFromIndices[key]

            #  Check that direction is normalized.
            vectorLength = (directionVector[0] ** 2 +
                            directionVector[1] ** 2 + directionVector[2] ** 2) ** .5

            #  If not normalized, normalize vector.
            if vectorLength != 1:
                directionVector[0] /= vectorLength
                directionVector[1] /= vectorLength
                directionVector[2] /= vectorLength

            # Get coordinates to place interstitial atom.
            interstitialCoordinates = []
            for key2, component in enumerate(directionVector):
                distanceTimesComponent = component * distance
                interstitialComponent = coordinatesAtIndex[key2] + distanceTimesComponent
                interstitialCoordinates.append(interstitialComponent)

            # Insert interstitial atom into structure at interstitialCoordinates.
            structure.append(interstitialAtom, interstitialCoordinates, True, True)

        # Create key.
        #  Add interstitialAtom to interstitialAtoms if not in list.
        if interstitialAtom not in interstitialAtoms:
            interstitialAtoms.append(interstitialAtom)

        # Find number of interstitials.
        for atom in interstitialAtoms:
            atom_count = 0
            for site in structure:
                species = site.specie
                if str(atom) == str(species):
                    atom_count += 1
            noInterstitials.append(atom_count)

        keyVal = self.structureKey(lattice_system, basic_atoms,
                                   dopantAtoms, noDopants, vacancyAtoms, noVacancies,
                                   noCells, interstitialAtoms, noInterstitials)

        return structure, keyVal

    def move_interstitial_dopant(self, origStruc, keyVal, atomsToMove, indicesToMove,
                                             distancesFromIndices, directionsFromIndices):
        """Move an interstitial atom to another location.

        """

        # Verify atomsToMove is a list.
        if isinstance(atomsToMove, list) is False:
            print "atomsToMove must be a list."
            sys.exit(1)

        if (len(indicesToMove) != len(distancesFromIndices)) or (
                    len(indicesToMove) != len(directionsFromIndices)) or (
                    len(distancesFromIndices) != len(directionsFromIndices)):
            print "indicesToMove, distanceFromIndices and directionsFromIndices must be the same length."
            sys.exit(1)

        # Create structure.
        structure = origStruc

        # Get interstitial atom species.  Currently Hardcoded to assume their is one only one interstitial atom species.
        interstitialAtom = atomsToMove[0]

        #  Verify that indicesToMove match atoms.
        for index in indicesToMove:
            species = structure[index].specie
            if str(interstitialAtom) != str(species):
                print "indicesToMove does not match atomsToMove"
                sys.exit(1)

        #  Get initial key values.
        lattice_system = keyVal[0]
        basic_atoms = keyVal[1]
        dopantAtoms = keyVal[2]
        noDopants = keyVal[3]
        vacancyAtoms = keyVal[4]
        noVacancies = keyVal[5]
        noCells = keyVal[6]
        interstitialAtoms = keyVal[7]
        noInterstitials = keyVal[8]

        #  Verify atomsToMove are interstitial.
        if interstitialAtom not in interstitialAtoms:
            print ("%s is not an interstitial atom." % str(interstitialAtom))
            sys.exit(1)

        #  Verify that you are not trying to move more interstitials than currently exist.
        if len(indicesToMove) > noInterstitials:
            print "IndicesToMove is greater than the number of interstitials."
            sys.exit(1)

        #  Insert interstitials.
        for key, index in enumerate(indicesToMove):
            coordinatesAtIndex = structure[index].coords
            distance = distancesFromIndices[key]
            directionVector = directionsFromIndices[key]

            #  Check that direction is normalized.
            vectorLength = (directionVector[0] ** 2 +
                            directionVector[1] ** 2 + directionVector[2] ** 2) ** .5

            #  If not normalized, normalize vector.
            if vectorLength != 1:
                directionVector[0] /= vectorLength
                directionVector[1] /= vectorLength
                directionVector[2] /= vectorLength

            # Get coordinates to move interstitial atom.
            interstitialCoordinates = []
            for key2, component in enumerate(directionVector):
                distanceTimesComponent = component * distance
                interstitialComponent = coordinatesAtIndex[key2] + distanceTimesComponent
                interstitialCoordinates.append(interstitialComponent)

            # Insert interstitial atom into structure at interstitialCoordinates.
            structure.append(interstitialAtom, interstitialCoordinates, True, True)

        structure.remove_sites(indicesToMove)

        #  Create unchanged key.
        keyVal = self.structureKey(lattice_system, basic_atoms,
                                   dopantAtoms, noDopants, vacancyAtoms, noVacancies,
                                   noCells, interstitialAtoms, noInterstitials)

        return structure, keyVal

