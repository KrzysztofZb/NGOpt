# deep learning models preparation
# at present works only with the Megnet library

from pymatgen.io.ase import AseAtomsAdaptor

from NGOptDBTools import *
from NGOptData import *
from NGOptIOTools import *

from megnet.models import MEGNetModel
from megnet.data.graph import GaussianDistance
from megnet.data.crystal import CrystalGraph

from keras.models import model_from_json

from collections import Counter
from ase.data import atomic_numbers, covalent_radii


def prepare_model_megnet(individuals, epochs, outfile, excl=[]):
    # prepares model file
    # prepares Megnet model based on list of individuals
    # uses total energy per atom
    # excl - excluding particular stoichiometry
    structures = []
    energies = []
    adapt = AseAtomsAdaptor()
    empty = 0
    if not excl:
        empty = 1

    i = 0
    for ind in individuals:
        struct_ase = ind.get_init_structure()
        chem_sym = struct_ase.get_chemical_symbols()
        e_tot = ind.e_tot
        struct_pymatgen = adapt.get_structure(struct_ase)
        flag = 1
        if empty == 0 and chem_sym == excl:
            flag = 0

        if flag == 1:
            structures.append(struct_pymatgen)
            energies.append(e_tot)
            i = i + 1

    print("read data of " + str(i) + " structures total")

    # standard vales as taken from Megnet manual
    nfeat_bond = 100
    nfeat_global = 2
    r_cutoff = 5
    gaussian_centers = np.linspace(0, r_cutoff + 1, nfeat_bond)
    gaussian_width = 0.5
    distance_converter = GaussianDistance(gaussian_centers, gaussian_width)
    graph_converter = CrystalGraph(bond_converter=distance_converter, cutoff=r_cutoff)
    model = MEGNetModel(nfeat_bond, nfeat_global, graph_converter=graph_converter)

    # model training
    model.train(structures, energies, epochs=epochs)

    model.save_model(outfile)


def prepare_dataset(symbols):
    # prepares dataset for simple model for lattice constant prediction
    occur = Counter(symbols)
    symbols_set = list(dict.fromkeys(symbols))

    nums = [0, 0, 0, 0]  # no. of atoms of each specie
    covs = [0., 0., 0., 0.]  # covalent radius
    ionics = [0., 0., 0., 0.]  # ionic radius
    electronegs = [0., 0., 0., 0.]  # electronegativity (Allen)
    ionizations = [0., 0., 0., 0.]  # ionization energy (kJ/mol)

    i = 0
    for s in symbols_set:
        n = occur[s]
        nums[i] = n
        Z = atomic_numbers[s]
        cov = covalent_radii[Z]
        covs[i] = cov
        ionic = ionic_radii[str(Z)]
        ionics[i] = ionic
        electroneg = electronegativity_Allen[str(Z)]
        electronegs[i] = electroneg
        ionization = ionization_energy[str(Z)]
        ionizations[i] = ionization
        i = i + 1

    dataset = [nums[0], nums[1], nums[2], nums[3],
               covs[0], covs[1], covs[2], covs[3],
               ionics[0], ionics[1], ionics[2], ionics[3],
               electronegs[0], electronegs[1], electronegs[2], electronegs[3],
               ionizations[0], ionizations[1], ionizations[2], ionizations[3]]

    return dataset


def get_amin_amax(chem_sym, stoichio, model_file):
    # predicts lattice constants based on stoichiometry
    json_file = open(model_file + '.json', 'r')
    loaded_model_json = json_file.read()
    json_file.close()
    loaded_model = model_from_json(loaded_model_json)

    # load weights into model
    loaded_model.load_weights(model_file + ".h5")
    symbols = get_symbols(chem_sym, stoichio)
    x = prepare_dataset(symbols)
    a = loaded_model.predict(np.array([x, ]))

    return a[0][0]-1.0, a[0][0]+0.5
