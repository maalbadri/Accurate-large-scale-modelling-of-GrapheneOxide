# Supplementary Code and Data

This repository includes data and code from the `Accurate large scale modelling of GrapheneOxide:  ion trapping and chaotropic potential atthe interface` paper. It gives a breakdown of using the `makegraphitics` code <sup>1</sup> for making the graphene-oxide sheet structure and forcefield. 

# Graphene Oxide Builder

MakeGraphitics is a library to ceate atomistic graphene oxide structures for molecular dynamics.

Output:
- .xyz 
- lammps data file

Automatically parameterise by OPLS forcefield.

## Install

Clone this repository. Install using Python2.7. Run the tests to check the installation has worked.
```
git clone https://github.com/velocirobbie/make-graphitics
cd make-graphitics
python setup.py install
pytest
```
A conda environment is provided if you do not have the right packages. If you have conda set up, execute these commands to create a working python environment before the install setp.
```
conda env create --file graphene-env.yml
conda activate graphene
```

## Instructions 

1) Make a rectangular graphene sheet that extends through periodic boundaries. Parameterised with OPLS and outputs to .xyz for easy veiwing with VMD and a LAMMPS data file.
```
python2.7 GO_rect.py
```
There are several tunable parameters in `GO_rect.py` that you may be interested in. Including:
- sheet dimensions
- C/O target ratio, `ratio`
- Rate at which new nodes are added, `new_island_freq`
- output snapshots of the oxidation process every N steps with `video_xyz=N`. Viewed in VMD with `topo readvarxyz out.xyz`

## DDEC parametrisation

DFT calculations of the above `.xyz` output can be performed using the ONETEP electronic structure calculation package. The ONETEP code version 6 is available from www.onetep.org. Instructions for implementing Density Derived Electrostatic and Chemical (DDEC) electron density partitioning are available [here](https://www.onetep.org/pmwiki/uploads/Main/Documentation/ddec.pdf).

The Lennard-Jones parameters can be calculated using the Tkatchenko-Scheffler relations using the [QUBEKit](https://github.com/qubekit/QUBEKit) package. Installation and instructions are available [here](https://github.com/qubekit/QUBEKit#qubekit-commands-custom-start-and-end-points-single-molecule). 

## Supplementary Data

The DDEC forcefield file for the DDEC GO used in this paper is provided in the `example_input` folder.

### Notes on the Oxidiser 

The Oxidiser takes a graphitic structure and attempts to oxidise it by the process described in (Yang, Angewandte Chemie, 2014; Sinclair, 2019). The algorithm proceeds as follows:

1) If hydrogens exist (i.e. edge of a flake), 1/4 are changed to alcohol groups and 1/4 to carboxyl groups (Lerf and Klinowski model). These values can be changed by passing the Oxidiser object the optional arguments: ` edge_OHratio = 0.25, edge_carboxyl_ratio = 0.5`.

2) The reactivity of every possible site is calculated. This is done by using a ranodom forest approach to extend the data set of GO reactivites given by Yang et al. We do not take into account the reactivity of the edges.

3) A site is oxidised at random weighted by each site's reactivity. The chance of an oxidisation producing an alcohol or epoxy group on the surface is by default 50:50, but can be specified by passing Oxidiser the optional argument: `surface_OHratio = 0.5`

4) The time elapsed between oxidations is estimated from the reactivity of the site that has been oxidised.

5) New nodes are added proportionally to `time_elapsed * new_island_freq`. Note this can be 0. The reasons for doing this are outlined in (Sinclair 2019).

6) Steps 2-5 are repeated until the target C/O ratio is reached or no new sites are available, usually C/O ~ 1.7 . We recommend setting the target, `ratio`, to over 2 as this is what is seen experimentally.

# Citing

<sup>1</sup>: The work contained here has been published in the following papers: 

 - Modeling Nanostructure in Graphene Oxide: Inhomogeneity and the Percolation Threshold https://doi.org/10.1021/acs.jcim.9b00114

 - Accurate large scale modelling of Graphene Oxide:  ion trapping and chaotropic potential at the interface (*in preparation*) 
```
@article{albadri2020accurate,
  title={Accurate large scale modelling of GrapheneOxide:  ion trapping and chaotropic potential atthe interface},
  author={al-Badri, Mohamed Ali and Smith, Paul and Sinclair, Robert C and al-Jamal, Khuloud and Lorenz, Christian D},
  journal = {Carbon},
  volume = {},
  number = {},
  pages = {},
  year = {2020},
  doi = {},
}
@misc{make-graphitics-github,
    url = {https://github.com/maalbadri/Accurate-large-scale-modelling-of-GrapheneOxide},
    howpublished = {\url{https://github.com/maalbadri/Accurate-large-scale-modelling-of-GrapheneOxide}},
    note = {Accessed: \today},
    author = {al-Badri, Mohamed Ali and Sinclair, Robert C.},
    year = {2020}
}
@article{sinclair2019modelling,
  title={Modelling nanostructure in graphene oxide: inhomogeneity and the percolation threshold},
  author={Sinclair, Robert Callum and Coveney, Peter Vivian},
  journal = {Journal of Chemical Information and Modeling},
  volume = {59},
  number = {6},
  pages = {2741-2745},
  year = {2019},
  doi = {10.1021/acs.jcim.9b00114},
}
@misc{make-graphitics-github,
    url = {https://github.com/velocirobbie/make-graphitics},
    howpublished = {\url{https://github.com/velocirobbie/make-graphitics}},
    note = {Accessed: \today},
    author = {Sinclair, Robert C.},
    year = {2019}
}
