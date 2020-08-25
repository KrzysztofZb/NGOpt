import numpy as np

from pymatgen.io.vasp import Poscar
from pymatgen.io.ase import AseAtomsAdaptor

from NGOptDBTools import get_c2db_prototypes
from NGOptDeepLearnTools import *

from megnet.models import MEGNetModel
from megnet.data.graph import GaussianDistance
from megnet.data.crystal import CrystalGraph

# individuals taken from c2db (https://cmr.fysik.dtu.dk/c2db/c2db.html)
individuals = get_c2db_prototypes(dbfile='c2db.db',prototype="",stability=True)
#print(individuals)

# prepare new model base on stable compounds from c2db, 1000 epochs
prepare_model_megnet(individuals, 1000, "model_1010", "")



