# NGOpt Copyright(C) 2020 Krzysztof Zberecki
# This program comes with ABSOLUTELY NO WARRANTY; for details type 'show w'.
# This is free software, and you are welcome to redistribute it under certain
# conditions; type 'show c' for details.

# path to vasp pseudo files
VASP_PSEUDO_PATH = "/path/to/your/vasp/pseudos"

# number of model evaluations during structure generation
RANDOM_ATTEMPTS = 100

# distance from substrate during structure generation on surface
SURF_DIST = 7.5
# lenght of cell to be put on surface
CELL_Z = 25.

# vasp, siesta and openmx exec paths
VASP_EXEC = "mpirun -np 8 /path/to/your/vasp_bin/vasp_std > out.tmp"
SIESTA_EXEC = "mpirun -np 8 /path/to/your/siesta_bin/siesta < siesta.fdf > out.tmp"
OPENMX_EXEC = "mpirun -np 8 /path/to/your/openmx_bin/openmx openmx.dat > out.tmp"

# minimal distance between atoms
DMIN = 1.2


