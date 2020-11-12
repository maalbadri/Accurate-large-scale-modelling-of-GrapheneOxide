import numpy as np
from base import Reactor


class EpoxyAzideSubstitution(Reactor):
    def __init__(self, conversion):
        assert 0 <= conversion <= 1, "conversion to azide must be between 0 and 1"
        self.conversion = conversion

    def react(self, sim):
        print "WARNING: Partial charges for azide groups are not approximate,"
        print "         we strongly suggest calculating charges for each system"

        sim.atom_labels = np.array(sim.atom_labels)
        n_epoxy = np.sum(sim.atom_labels == 6)
        assert n_epoxy > 0, "no epoxies exist, graphene must be oxidised to substitute"

        n_azide_to_add = int(round(n_epoxy * self.conversion))

        sim.bond_graph = sim.generate_bond_graph(sim.bonds)

        for i in range(n_azide_to_add):
            site = self.find_site(sim)
            sim = self.substitute(sim, site)

        print "Added ", n_azide_to_add, "azides"
        sim.generate_connections()
        sim.vdw_defs.update(
            {12: 699, 13: 254, 14: 204, 15: 204}  # Tertiary alkyl nitrile  # NC  # NZ
        )
        sim.validate()
        return sim

    def find_site(self, sim):
        epoxies_loc = np.where(sim.atom_labels == 6)
        epoxy_loc = np.random.choice(epoxies_loc[0])
        return epoxy_loc

    def substitute(self, sim, epoxy_loc):
        carbons = list(sim.bond_graph[epoxy_loc])  # carbons connected to epoxy
        np.random.shuffle(carbons)  # random which carbon has nucleophillic attack
        C_a_loc, C_b_loc = carbons

        # assumes sheet is in xy plane
        above = 0
        if sim.coords[epoxy_loc][2] - sim.coords[C_a_loc][2] > 0:
            above = 1
        else:
            above = -1

        CH_bond = 1
        CO_bond = 1.4
        # move epoxy oxygen to alcohol position
        cc_vector = sim.coords[C_a_loc] - sim.coords[C_b_loc]
        sim.coords[epoxy_loc] += 0.5 * cc_vector  # move O above carbon a
        sim.coords[epoxy_loc, 2] = CO_bond * above

        H_coord = np.array([sim.coords[epoxy_loc] + np.array([0, 0, above * CH_bond])])
        N_coords = np.array(
            [
                sim.coords[C_b_loc] + [0, 0, -1.45 * above],
                sim.coords[C_b_loc] + [0, 0, -2.69 * above],
                sim.coords[C_b_loc] + [0, 0, -3.81 * above],
            ]
        )
        sim.coords = np.vstack((sim.coords, H_coord, N_coords))

        H_id = len(sim.atom_labels) + 1
        N1_id = H_id + 1
        N2_id = H_id + 2
        N3_id = H_id + 3

        epoxy_id = epoxy_loc + 1
        C_a_id = C_a_loc + 1
        C_b_id = C_b_loc + 1

        # remove OCb bond
        removed = False
        for i, bond in enumerate(sim.bonds):
            if list(bond) == [epoxy_id, C_b_id] or list(bond) == [C_b_id, epoxy_id]:
                sim.bonds = np.delete(sim.bonds, i, 0)
                removed = True
                break
        if not removed:
            raise RuntimeError

        new_bonds = np.array(
            [[epoxy_id, H_id], [C_b_id, N1_id], [N1_id, N2_id], [N2_id, N3_id]]
        )
        sim.bonds = np.vstack((sim.bonds, new_bonds))

        # reasign types, charges
        def reassign_atom(loc, label, charge):
            sim.atom_labels[loc] = label
            sim.atom_charges[loc] = charge

        """
        C-OH
        3 = Ct, tertiary C-OH
        4 = Oa, C-OH
        5 = Ha, C-OH

        C-N3
        12 = C
        13 = N1
        14 = N2
        15 = N3

        azide charges doi:10.1007/s00894-009-0488-z
        these do not sum to zero,
        for simplicity we subtrcat the difference from the azide carbon
        """

        reassign_atom(epoxy_loc, 4, -0.585)
        reassign_atom(C_a_loc, 3, 0.15)
        reassign_atom(C_b_loc, 12, 0.2708)

        # H, N1, N2, N3
        sim.atom_labels = np.append(sim.atom_labels, [5, 13, 14, 15])
        sim.atom_charges = np.append(
            sim.atom_charges, [0.435, -0.4926, 0.6096, -0.3878]
        )
        this_molecule = sim.molecule_labels[epoxy_loc]
        sim.molecule_labels = np.append(sim.molecule_labels, [this_molecule] * 4)

        return sim
