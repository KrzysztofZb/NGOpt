# generation of random structures

from random import random, uniform, randrange
import random
import string

from ase import Atoms, Atom
from ase.build import surface

from pymatgen.io.ase import AseAtomsAdaptor

from megnet.models import MEGNetModel

from pyxtal.crystal import random_crystal_2D
from pyxtal.crystal import Tol_matrix

from NGOptIOTools import *
from NGOptConfig import *

# from memory_profiler import profile


def crystal_to_atoms(crystal):
    # pyxtal crystal class object to ase Atoms class object, returns ase Atoms
    if os.path.isfile('POSCAR_TMP'):
        os.system("rm POSCAR_TMP")
    crystal.to_file(fmt="poscar", filename="POSCAR_TMP")
    iot = IOTools(calculator="VASP")
    struct = iot.read_structure(structfile="POSCAR_TMP")
    os.system("rm POSCAR_TMP")
    return struct


def random_str(lenght=10):
    # returns a random string of fixed length
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(lenght))


def dist(a, b):  # returns distance of two vectors
    d = np.sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]))
    return d


def check_dist(Nat, positions, dmin):
    # checks if atoms' positions fullfill dmin distance constrain, returns flag 1/0
    flag = 1
    for i in range(Nat):
        for j in range(Nat):
            if i != j:
                ds = dist(positions[i], positions[j])
                if ds < dmin:
                    flag = 0
    return flag


def distz(z1, z2):
    return np.abs(z1 - z2)


def random_cell(amin, amax, rnd=0):
    # returns random cell, cubic or hexagonal
    cell_types = {1: "cubic_simple", 2: "hexagonal"}
    indx = 0
    if rnd == 0:
        indx = randrange(2) + 1
    if rnd > 0:
        indx = rnd
    # cubic    
    if cell_types[indx] == "cubic_simple":
        a = uniform(amin, amax)
        b = uniform(amin, amax)
        a1 = [a / 1., 0., 0.]
        a2 = [0., b / 1, 0.]
        a3 = [0., 0., 20.]
        cell = [a1, a2, a3]
    # hexagonal
    if cell_types[indx] == "hexagonal":
        a = uniform(amin, amax)
        a1 = [a, 0., 0.]
        a2 = [-0.5 * a, a * np.sqrt(3.) / 2., 0.]
        a3 = [0., 0., 20.]
        cell = [a1, a2, a3]

    return cell


def random_structure(stoichio, amin, amax, dmin, iwrite=0):
    # returns pure random structure (ase Atoms)
    struct = Atoms(stoichio)
    cell = random_cell(amin, amax)
    struct.set_cell(cell)

    symbols = struct.get_chemical_symbols()

    Nat = len(symbols)

    flag = 0
    niter = 0
    while flag == 0:
        positions = []

        for i in range(Nat):
            pos = [uniform(0., 0.90), uniform(0., 0.90), uniform(-0.2, 0.2)]
            positions.append(pos)

        struct.set_scaled_positions(positions)
        struct2x2x1 = struct * (2, 2, 1)
        positions_ang = struct2x2x1.get_positions()
        flag = check_dist(Nat * 2 * 2, positions_ang, dmin)

        niter = niter + 1

    if flag == 1 and iwrite == 1:
        write(filename="POSCAR", images=struct, format="vasp")
        write(filename="POSCAR.in", images=struct, format="espresso-in")

    return struct


# @profile
def random_structure_model(stoichio, amin, amax, dmin, model_file, Natt=RANDOM_ATTEMPTS):
    # returns random structure (ase Atoms) with lowest e_tot according to megnet model
    adapt = AseAtomsAdaptor()
    model = MEGNetModel.from_file(model_file)
    e_tot_min = 0.
    flag = 0

    for i in range(Natt):
        struct = random_structure(stoichio, amin, amax, dmin)
        # Nat = len(struct.get_chemical_symbols())
        struct_pymatgen = adapt.get_structure(struct)
        try:
            e_tot = model.predict_structure(struct_pymatgen)
        except:
            e_tot = 0.
            print("isolated molecule exception handled")
        if e_tot < e_tot_min:
            struct_out = struct
            e_tot_min = e_tot
            flag = 1
    if flag == 0:
        print("Warning: structure not generated!")
        struct_out = Atoms(stoichio)
    if flag == 1:
        print("e_tot/atom: " + str(e_tot_min))
        write(filename='POSCAR_' + random_str(5) + '.in', images=struct_out, format="espresso-in")

    del model

    return struct_out


def random_structure_group(symbols, composition, thickness, tol_factor, model_file, dmin=2.0, Natt=RANDOM_ATTEMPTS):
    # returns pyxtal generated structure (ase Atoms) with lowest e_tot according to megnet model
    tol_m_1 = Tol_matrix(prototype="atomic", factor=tol_factor)
    adapt = AseAtomsAdaptor()
    model = MEGNetModel.from_file(model_file)
    e_tot_min = 0.

    for i in range(Natt):
        group_id = randrange(80) + 1
        my_crystal = random_crystal_2D(group_id, symbols, composition, 1.0, thickness=thickness, tm=tol_m_1)
        flag = 0
        if my_crystal.valid == True:
            struct = crystal_to_atoms(my_crystal)
            Nat = len(struct.get_chemical_symbols())
            struct_pymatgen = adapt.get_structure(struct)
            try:
                e_tot = model.predict_structure(struct_pymatgen)
            except:
                e_tot = 0.
                print("isolated molecule exception handled")
            struct2x2x1 = struct * (2, 2, 1)
            positions = struct2x2x1.get_positions()
            # positions = struct.get_positions()
            # print(struct)
            flag = check_dist(Nat * 2 * 2, positions, dmin)
            # print(flag)

            if (e_tot < e_tot_min) and flag == 1:
                struct_out = struct
                e_tot_min = e_tot

    print("e_tot/atom: " + str(e_tot_min))
    # write(filename="POSCAR.best.in", images=struct_out, format="espresso-in")
    # write(filename="POSCAR.best", images=struct_out, format="vasp")

    del model

    return struct_out


def random_structure_on_substrate(symbols, amin, amax, dmin, model_file, Natt=RANDOM_ATTEMPTS):
    # returns random structure (ase Atoms) on substrate with lowest e_tot according to Megnet model
    substrate = read_vasp("POSCAR.substrate")
    adapt = AseAtomsAdaptor()
    model = MEGNetModel.from_file(model_file)
    e_tot_min = 1000.

    for i in range(Natt):
        s = surface(substrate, (0, 0, 1), 1, vacuum=0., tol=1e-10)
        cell = s.get_cell()
        cell[2][2] = CELL_Z
        s.set_cell(cell)
        amin = cell[0][0]
        amax = cell[0][0]
        struct = random_structure(symbols, amin, amax, dmin, iwrite=0)

        j = 0
        atoms = struct.get_chemical_symbols()
        positions = struct.get_positions()
        for atom in atoms:
            at = Atom(atom)
            positions[j][2] = positions[j][2] + SURF_DIST
            pos = positions[j]
            at.position = pos
            s.append(at)
            j = j + 1

        struct_pymatgen = adapt.get_structure(s)
        try:
            e_tot = model.predict_structure(struct_pymatgen)
            # print(e_tot)
        except:
            e_tot = 0.
            print("isolated molecule exception handled")
        if e_tot < e_tot_min:
            struct_out = s
            e_tot_min = e_tot

    print("e_tot min: ", e_tot_min)
    write(filename='best.in', images=struct_out, format="espresso-in")

    del model

    return struct_out
