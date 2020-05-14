from NGOptDriver import *
from NGOptDeepLearnTools import *

from pymatgen.io.ase import AseAtomsAdaptor
from megnet.models import MEGNetModel

chem_sym = ['Mo', 'S']
stoichio = [1, 2]

driver = Driver(Npop=24, chem_sym=chem_sym, stoichio=stoichio, model_file="megnet_model_1000", calculator="VASP")

amin, amax  = get_amin_amax(chem_sym=chem_sym, stoichio=stoichio, model_file="lattice_model_1000")

print("lattice constants taken from interval: (" + str(amin) + "," + str(amax) + ")")
driver.amin = amin
driver.amax = amax
driver.dmin = 1.6

driver.exec = "mpirun -np 8 /opt/vasp/5.4.4.intel.mpi/vasp_std > out.tmp"

populations = []
Ngen = 10
population = driver.init_population_random(model=1)
populations.append(population)

for i in range(Ngen):
    driver.evaluate_population(populations[i], "POP" + str(i + 1))
    driver.population_report(populations[i], i+1)
    new_population = driver.generate_new_population_alg1(populations[i])
    populations.append(new_population)

os.system("cp populations.db populations.finished.db")

