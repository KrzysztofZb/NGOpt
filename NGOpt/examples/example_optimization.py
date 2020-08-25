import numpy as np

from NGOptDriver import *
from NGOptDeepLearnTools import *
from NGOptIndividual import Individual

from pymatgen.io.ase import AseAtomsAdaptor
from megnet.models import MEGNetModel

# MoS2
chem_sym = ['Mo', 'S']
stoichio = [1, 2]

# uncomment for CoGa2S4
#chem_sym = ['Co', 'Ga', 'S']
#stoichio = [1, 2, 4]

# uncomment for Cr2Te3
#chem_sym = ['Cr', 'Te']
#stoichio = [2, 3]

# driver - governs the GA algorithm
driver = Driver(Npop=24, chem_sym=chem_sym, stoichio=stoichio, model_file="megnet_model_1000", calculator="VASP")

# amin and amax are taken from a simple model, which estimates the lattice constants
amin, amax  = get_amin_amax(chem_sym=chem_sym, stoichio=stoichio, model_file="lattice_model_1000")

print("lattice constants taken from interval: (" + str(amin) + "," + str(amax) + ")")
driver.amin = amin
driver.amax = amax
driver.dmin = 1.6

# change path to match your vasp path
driver.exec = "mpirun -np 8 vasp > out.tmp"


populations = []
Ngen = 20
population = driver.init_population_random(model=1)
populations.append(population)

# main GA loop
for i in range(Ngen):
    driver.evaluate_population(populations[i], "POP" + str(i + 1))
    driver.population_report(populations[i], i+1)
    new_population = driver.generate_new_population_alg1(populations[i])
    populations.append(new_population)

# data from optimization run are collected in a database
os.system("cp populations.db populations.finished.db")

