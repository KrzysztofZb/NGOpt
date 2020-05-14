# path to vasp pseudo files
VASP_PSEUDO_PATH = "/home/kz/projects/new_materials_2D_2018/pseudos/vasp_new2"

# number of model evaluations during structure generation
RANDOM_ATTEMPTS = 100

# distance from substrate during structure generation on surface
SURF_DIST = 7.5
# lenght of cell to be put on surface
CELL_Z = 25.

# vasp, siesta and openmx exec paths
VASP_EXEC = "mpirun -np 8 /opt/vasp/5.4.4.intel.mpi/vasp_std > out.tmp"
SIESTA_EXEC = "mpirun -np 8 /opt/siesta/4.1-b4.gnu.mpi/bin/siesta < siesta.fdf > out.tmp"
OPENMX_EXEC = ""

# minimal distance between atoms
DMIN = 1.2


