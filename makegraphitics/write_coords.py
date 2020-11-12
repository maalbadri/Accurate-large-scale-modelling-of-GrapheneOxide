import os
import numpy as np


class Writer(object):
    def __init__(self, sim, system_name="comment line"):
        # Takes a numpy 3xN array of atom coordinates and outputs
        # them in different formats for viewing/modelling
        self.coords = sim.coords
        self.atom_labels = sim.atom_labels
        self.molecule = sim.molecule_labels
        self.charges = sim.atom_charges

        self.bonds = sim.bonds
        self.bond_labels = sim.bond_labels
        self.angles = sim.angles
        self.angle_labels = sim.angle_labels
        self.dihedrals = sim.dihedrals
        self.dihedral_labels = sim.dihedral_labels
        self.impropers = sim.impropers
        self.improper_labels = sim.improper_labels

        self.box_dimensions = sim.box_dimensions
        self.system_name = system_name

        if hasattr(sim, "atom_ids"):
            self.atom_ids = sim.atom_ids
        if hasattr(sim, "masses"):
            self.masses = sim.masses
        if hasattr(sim, "bond_coeffs"):
            self.bond_coeffs = sim.bond_coeffs
        if hasattr(sim, "angle_coeffs"):
            self.angle_coeffs = sim.angle_coeffs
        if hasattr(sim, "dihedral_coeffs"):
            self.dihedral_coeffs = sim.dihedral_coeffs
        if hasattr(sim, "improper_coeffs"):
            self.improper_coeffs = sim.improper_coeffs
        if hasattr(sim, "pair_coeffs"):
            self.pair_coeffs = sim.pair_coeffs

        type_lists = {
            "atom": "pair_coeffs",
            "bond": "bond_coeffs",
            "angle": "angle_coeffs",
            "dihedral": "dihedral_coeffs",
            "improper": "improper_coeffs",
        }
        for type_list in type_lists:
            maxlabel = max(getattr(self, type_list + "_labels"))
            coeff = type_lists[type_list]
            if hasattr(sim, coeff):
                maxcoeff = max(getattr(self, coeff))
            else:
                maxcoeff = 0
            setattr(self, "n" + type_list + "_types", max(maxcoeff, maxlabel))

    def write_xyz(self, filename="out.xyz", option="w"):
        with open(filename, option) as outfile:
            outfile.write(str(len(self.coords)) + "\n" + self.system_name + "\n")
            for i in range(len(self.coords)):
                xyz = (
                    str(self.coords[i][0])
                    + " "
                    + str(self.coords[i][1])
                    + " "
                    + str(self.coords[i][2])
                )
                if hasattr(self, "masses"):
                    try:
                        mass = self.masses[self.atom_labels[i]]
                    except KeyError:
                        mass = 100
                    if abs(mass - 12.0) < 0.5:
                        atom_label = "C "
                    elif abs(mass - 1.0) < 0.5:
                        atom_label = "H "
                    elif abs(mass - 14.0) < 0.5:
                        atom_label = "N "
                    elif abs(mass - 16.0) < 0.5:
                        atom_label = "O "
                    else:
                        atom_label = str(self.atom_labels[i]) + " "
                else:
                    atom_label = str(self.atom_labels[i]) + " "

                outfile.write(atom_label + xyz + "\n")
            print "Coords written to " + str(filename)

    def write_reaxff(self, filename="data.lammps"):
        # atom_type charge
        masses = np.unique(self.masses.values())
        nreax_types = len(masses)
        reax_types = {12.011: 1, 1.008: 2, 15.999: 3, 14.007: 4, 28.086: 5}
        with open(filename, "w") as outfile:
            outfile.write(
                "# "
                + self.system_name
                + "\n"
                + str(len(self.coords))
                + " atoms \n"
                + "\n"
                + str(nreax_types)
                + " atom types \n"
                + "\n"
                + str(self.box_dimensions[0, 0])
                + "\t"
                + str(self.box_dimensions[0, 1])
                + "\t xlo xhi \n"
                + str(self.box_dimensions[1, 0])
                + "\t"
                + str(self.box_dimensions[1, 1])
                + "\t ylo yhi \n"
                + str(self.box_dimensions[2, 0])
                + "\t"
                + str(self.box_dimensions[2, 1])
                + "\t zlo zhi \n"
                + "\n"
            )
            if hasattr(self, "masses"):
                outfile.write("\n Masses \n \n")
                for mass in reax_types:
                    outfile.write(str(reax_types[mass]) + "\t" + str(mass) + "\n")

            outfile.write("\n Atoms \n \n")

            for i in range(len(self.coords)):
                atom_type = self.atom_labels[i]
                reax_type = reax_types[self.masses[atom_type]]
                outfile.write(
                    str(i + 1)
                    + "\t "
                    + str(reax_type)  # atom ID
                    + "\t "
                    + str(self.charges[i])  # atom type
                    + "\t "
                    + str(self.coords[i][0])  # atomcharg
                    + "\t "
                    + str(self.coords[i][1])  # x
                    + "\t "
                    + str(self.coords[i][2])  # y
                    + "\n "  # z
                )

            print "Coords written to " + filename

    def write_lammps(self, filename="data.lammps"):
        # atom_type full
        with open(filename, "w") as outfile:
            outfile.write(
                "# "
                + self.system_name
                + "\n"
                + str(len(self.coords))
                + " atoms \n"
                + str(len(self.bonds))
                + " bonds \n"
                + str(len(self.angles))
                + " angles \n"
                + str(len(self.dihedrals))
                + " dihedrals \n"
                + str(len(self.impropers))
                + " impropers \n"
                "\n"
                + str(self.natom_types)
                + " atom types \n"
                + str(self.nbond_types)
                + " bond types \n"
                + str(self.nangle_types)
                + " angle types \n"
                + str(self.ndihedral_types)
                + " dihedral types \n"
                + str(self.nimproper_types)
                + " improper types \n"
                + "\n"
                + str(self.box_dimensions[0, 0])
                + "\t"
                + str(self.box_dimensions[0, 1])
                + "\t xlo xhi \n"
                + str(self.box_dimensions[1, 0])
                + "\t"
                + str(self.box_dimensions[1, 1])
                + "\t ylo yhi \n"
                + str(self.box_dimensions[2, 0])
                + "\t"
                + str(self.box_dimensions[2, 1])
                + "\t zlo zhi \n"
                + "\n"
            )
            if hasattr(self, "masses"):
                outfile.write("\n Masses \n \n")
                for i in self.masses:
                    outfile.write(str(i) + "\t" + str(self.masses[i]) + "\n")

            if hasattr(self, "pair_coeffs"):
                outfile.write("\n Pair Coeffs \n \n")
                for i in self.pair_coeffs:
                    outfile.write(
                        str(i)
                        + "\t"
                        + str(self.pair_coeffs[i][1])
                        + "\t"
                        + str(self.pair_coeffs[i][2])
                        + "\n"
                    )

            if hasattr(self, "bond_coeffs"):
                outfile.write("\n Bond Coeffs \n \n")
                for i in self.bond_coeffs:
                    outfile.write(
                        str(i)
                        + "\t"
                        + str(self.bond_coeffs[i][1])
                        + "\t"
                        + str(self.bond_coeffs[i][2])
                        + "\n"
                    )
            if hasattr(self, "angle_coeffs"):
                outfile.write("\n Angle Coeffs \n \n")
                for i in self.angle_coeffs:
                    outfile.write(
                        str(i)
                        + "\t"
                        + str(self.angle_coeffs[i][1])
                        + "\t"
                        + str(self.angle_coeffs[i][2])
                        + "\n"
                    )
            if hasattr(self, "dihedral_coeffs"):
                outfile.write("\n Dihedral Coeffs \n \n")
                for i in self.dihedral_coeffs:
                    outfile.write(
                        str(i)
                        + "\t"
                        + str(self.dihedral_coeffs[i][1])
                        + "\t"
                        + str(self.dihedral_coeffs[i][2])
                        + "\t"
                        + str(self.dihedral_coeffs[i][3])
                        + "\t"
                        + str(self.dihedral_coeffs[i][4])
                        + "\n"
                    )

            if hasattr(self, "improper_coeffs"):
                outfile.write("\n Improper Coeffs \n \n")
                for i in self.improper_coeffs:
                    outfile.write(
                        str(i)
                        + "\t"
                        + str(self.improper_coeffs[i][1])
                        + "\t"
                        + str(self.improper_coeffs[i][2])
                        + "\n"
                    )

            outfile.write("\n Atoms \n \n")

            if hasattr(self, 'atom_ids'):
                ids = np.argsort(self.atom_ids)
            else:
                ids = np.arange(len(self.coords))

            for serial, i in enumerate(ids):
                outfile.write(
                    str(serial+1)
                    + "\t "
                    + str(self.molecule[i])  # atom ID
                    + "\t "
                    + str(self.atom_labels[i])  # molecule ID
                    + "\t "
                    + str(self.charges[i])  # atom type
                    + "\t "
                    + str(self.coords[i][0])  # atomcharg
                    + "\t "
                    + str(self.coords[i][1])  # x
                    + "\t "
                    + str(self.coords[i][2])  # y
                    + "\n "  # z
                )

            if len(self.bonds):
                outfile.write("\n Bonds \n \n")
                for i in range(len(self.bonds)):
                    outfile.write(
                        str(i + 1)
                        + "\t "
                        + str(self.bond_labels[i])  # bond ID
                        + "\t "
                        + str(self.bonds[i][0])
                        + "\t "
                        + str(self.bonds[i][1])  # atom 1
                        + "\n"  # atom 2
                    )

            if len(self.angles):
                outfile.write("\n Angles \n \n")
                for i in range(len(self.angles)):
                    outfile.write(
                        str(i + 1)
                        + "\t "
                        + str(self.angle_labels[i])  # angle ID
                        + "\t "
                        + str(self.angles[i][0])
                        + "\t "
                        + str(self.angles[i][1])
                        + "\t "
                        + str(self.angles[i][2])
                        + "\n"
                    )

            if len(self.dihedrals):
                outfile.write("\n Dihedrals \n \n")
                for i in range(len(self.dihedrals)):
                    outfile.write(
                        str(i + 1)
                        + "\t"
                        + str(self.dihedral_labels[i])
                        + " \t"
                        + str(self.dihedrals[i][0])
                        + " \t"
                        + str(self.dihedrals[i][1])
                        + " \t"
                        + str(self.dihedrals[i][2])
                        + " \t"
                        + str(self.dihedrals[i][3])
                        + " \n"
                    )

            if len(self.impropers):
                outfile.write("\n Impropers \n \n")
                for i in range(len(self.impropers)):
                    outfile.write(
                        str(i + 1)
                        + "\t"
                        + str(self.improper_labels[i])
                        + "\t"
                        + str(self.impropers[i][0])
                        + " \t"
                        + str(self.impropers[i][1])
                        + " \t"
                        + str(self.impropers[i][2])
                        + " \t"
                        + str(self.impropers[i][3])
                        + " \n"
                    )

            print "Coords written to " + filename
