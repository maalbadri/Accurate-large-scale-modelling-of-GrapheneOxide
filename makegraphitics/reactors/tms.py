import numpy as np
from base import Reactor


class CarboxylTMSSubstitution(Reactor):
    def __init__(self, conversion):
        assert 0 <= conversion <= 1, "conversion to azide must be between 0 and 1"
        self.conversion = conversion

    def react(self, sim):
        print "WARNING: Partial charges for TMS groups are approximate,"
        print "         we strongly suggest calculating charges for each system"

        sim.atom_labels = np.array(sim.atom_labels)
        n_carboxyl = np.sum(sim.atom_labels == 8)
        assert n_carboxyl > 0, "no carboxyls exist, graphene must be oxidised to substitute"

        n_tms_to_add = int(round(n_carboxyl * self.conversion))

        sim.bond_graph = sim.generate_bond_graph(sim.bonds)

        for i in range(n_tms_to_add):
            site = self.find_site(sim)
            sim = self.substitute(sim, site)

        print "Added ", n_tms_to_add, "tms"
        sim.generate_connections()
        sim.vdw_defs.update(
                {16: 412, 17: 407, 18: 408, 19: 431,
                 20: 410, 21: 769, 22: 866, 23: 871,
                 24: 870, 25: 413}
        )
        sim.validate()
        return sim

    def find_site(self, sim):
        carboxyls_loc = np.where(sim.atom_labels == 8)
        carboxyl_loc = np.random.choice(carboxyls_loc[0])
        return carboxyl_loc

    def find_neighbour_of_type(self, sim, centre, type_):
        neighbours = sim.bond_graph[centre]
        found = False
        for n in neighbours:
            if sim.atom_labels[n] == type_:
                found = True
                break
        assert found
        return n

    def substitute(self, sim, carboxyl_loc):
        o_acid_loc = self.find_neighbour_of_type(sim, carboxyl_loc, 10)
        o_ketone_loc = self.find_neighbour_of_type(sim, carboxyl_loc, 9)
        h_acid_loc = self.find_neighbour_of_type(sim, o_acid_loc, 5)
        c_benz_loc = self.find_neighbour_of_type(sim, carboxyl_loc, 11)

        oh_vector = sim.coords[h_acid_loc] - sim.coords[o_acid_loc]
        oh_vector /= np.linalg.norm(oh_vector)

        """
        8: 209,  # Cc, Carboxylic carbon
        9: 210,  # Oc, Ketone oxygen
        10: 211,  # Oa, alcohol
        11: 108,  # Cb, Benzyl


        C-C=O -O-C-C(triple)C-Si-CH3
        16 = C-OOR 412
        17 = O=C   407
        18 = CO-O-R 408
        19 = OO-C- 431
        20 = OOC-H 410
        21 = Alkyne 769
        22 = Si 866
        23 = CH3-Si 871
        24 = SiC-H 870 # bond type 46 (OPLS had it listed as 45?)
        25 = C-COOR 413
        """

        tms = np.array([[[0, 1, 0], [20]],   # H 0
                        [[0, -1, 0], [20]],    # H 1
                        [[1.4, 0, 0], [21]],   # C 2
                        [[2.8, 0, 0], [21]],   # C 3
                        [[4.2, 0, 0], [22]],   # Si 4
                        [[4.2, 1, 0], [23]],   # C 5
                        [[5.2, 1, 0], [24]],   # H 6
                        [[4.2, 2, 0], [24]],   # H 7
                        [[4.2, 1, 1], [24]],  # H 8
                        [[4.2, -1, 0], [23]],   # C 9
                        [[5.2, -1, 0], [24]],   # H 10
                        [[4.2, -2, 0], [24]],   # H 11
                        [[4.2, -1, 1], [24]],  # H 12
                        [[4.2, 0, -1], [23]],  # C 13
                        [[5.2, 0, -1], [24]],  # H 14
                        [[4.2, 0, -2], [24]],  # H 15
                        [[3.2, 0, -1], [24]]])  # H 16
        tms_template_coords = np.array([np.array(coord) for coord in tms[:, 0]])
        tms_template_vector = np.array([1, 0, 0])

        costheta = np.dot(oh_vector, tms_template_vector)
        sintheta = np.linalg.norm(np.cross(oh_vector, tms_template_vector))

        rot_mat = np.array([[costheta, -sintheta, 0],
                            [sintheta, costheta,  0],
                            [0,        0,         1]])
        tms_rot_coords = np.empty(tms_template_coords.shape)

        for i in range(len(tms_template_coords)):
            tms_rot_coords[i] = np.dot(rot_mat, tms_template_coords[i])

        # correct for 360 rotation
        if oh_vector[1] < 0:
            tms_rot_coords[:, 1] *= -1

        tms_coords = tms_rot_coords + sim.coords[h_acid_loc]

        sim.coords = np.vstack((sim.coords, tms_coords))

        offset = np.max(sim.bonds) + 1
        new_bonds1 = np.array([[h_acid_loc + 1, 0 + offset],
                               [h_acid_loc + 1, 1 + offset],
                               [h_acid_loc + 1, 2 + offset]])
        new_bonds2 = np.array([[2, 3], [3, 4],
                               [4, 5], [4, 9], [4, 13],
                               [5, 6], [5, 7], [5, 8],
                               [9, 10], [9, 11], [9, 12],
                               [13, 14], [13, 15], [13, 16]])
        new_bonds = np.vstack((new_bonds1, new_bonds2 + offset))
        sim.bonds = np.vstack((sim.bonds, new_bonds))

        def reassign_atom(loc, label, charge):
            sim.atom_labels[loc] = label
            sim.atom_charges[loc] = charge

        reassign_atom(carboxyl_loc, 16, 0.625)
        reassign_atom(o_ketone_loc, 17, -0.430)
        reassign_atom(o_acid_loc, 18, -0.330)
        reassign_atom(h_acid_loc, 19, 0.190)
        reassign_atom(c_benz_loc, 25, 0.135)

        # These charges have been fudged to balance to zero
        # OPLS has not implememnted alyne si bonds
        sim.atom_labels = np.append(sim.atom_labels,
                                 np.array([np.array(xi) for xi in tms[:, 1]])
                )
        sim.atom_charges = np.append(
            sim.atom_charges, [0.030, 0.030, -0.210, -0.210, 1.400]
        )
        sim.atom_charges = np.append(
            sim.atom_charges, [-0.2, -0.07, -0.07, -0.07]*3
        )
        this_molecule = sim.molecule_labels[carboxyl_loc]
        sim.molecule_labels = np.append(sim.molecule_labels,
                                        [this_molecule] * len(tms))

        return sim
