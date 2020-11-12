import makegraphitics as mg

motif = mg.molecules.Rectangle_Graphene(50, 30)
flake = mg.Crystal(motif, [1, 1, 1])

oxidiser = mg.reactors.Oxidiser(
    ratio=2.2, video_xyz=20, new_island_freq=1e15, method="rf", 
    edge_carboxyl_ratio=0.1
)
flake = oxidiser.react(flake)

mg.Parameterise(flake, flake.vdw_defs)

name = "go_50x30"
output = mg.Writer(flake, name)
output.write_xyz(name + ".xyz")
output.write_lammps(name + ".data")
output.write_reaxff(name + "reax.data")
