# TODO: remove hardcoded parameters

# genetic operations

from NGOptIndividual import *
from NGOptRandom import *

# from memory_profiler import profile


def swap_positions(list, pos1, pos2):
    # swaps two positions in a list, returns swapped list
    list_new = list[:]
    l1 = list_new[pos1]
    l2 = list_new[pos2]
    list_new[pos2] = l1
    list_new[pos1] = l2
    return list_new


def minz(zmax, zs):
    # finds minimum but ignores all zs < zmax
    zmin = 100.
    for z in zs:
        if z > zmax and z < zmin:
            zmin = z
    return zmin


def get_layers(struct):
    pos = struct.get_positions()
    chem = struct.get_chemical_symbols()
    debug = 0
    #
    zs = []
    for p in pos:
        z = p[2]
        zs.append(z)
    #
    zmin = min(zs)
    zmaxx = max(zs) + 10.
    dlay = 2.5
    zmax = zmin + dlay
    n = len(zs)
    layers = []
    #
    nl = 1
    while 1:
        tmp = []
        layer = []
        n = len(zs)
        for i in range(n):
            z = zs[i]
            if z >= zmin and z <= zmax:
                if debug:
                    print("atom in layer " + str(nl) + " " + str(z) + " indx: " + str(i) + " " + str(chem[i]))
                tmp.append(z)
                layer.append(i)
        layers.append(layer)
        maxlay = max(tmp)
        zmin = minz(maxlay, zs)
        zmax = zmin + dlay
        nl = nl + 1
        if zmax > zmaxx:
            break

    return layers


# @profile
def new_by_switch_atoms(ind_in: Individual):
    # returns individual made by switch of two random atoms
    struct = ind_in.get_relaxed_structure().copy()
    # print(struct)
    chem_sym = struct.get_chemical_symbols()
    Nat = len(chem_sym)
    if Nat > 2:
        indx1 = randrange(Nat)
        indx2 = indx1
        flag = 0
        while flag == 0:
            indx2 = randrange(Nat)
            if indx2 != indx1 and chem_sym[indx1] != chem_sym[indx2]:
                # print(indx1,indx2)
                chem_sym_new = swap_positions(chem_sym, indx1, indx2)
                struct.set_chemical_symbols(chem_sym_new)
                flag = 1
    ind = Individual(struct, -1)
    return ind


# @profile
def new_by_switch_atoms_on_substrate(ind_in: Individual):
    # returns individual made by switch of two random atoms on particular substrate
    struct = ind_in.get_relaxed_structure().copy()
    substrate = read_vasp("POSCAR.substrate")
    s = surface(substrate, (0, 0, 1), 1, vacuum=0., tol=1e-10)
    cell = s.get_cell()
    cell[2][2] = 25.
    s.set_cell(cell)
    z_max = -100.
    positions = s.get_positions()
    for pos in positions:
        z = pos[2]
        if z > z_max:
            z_max = z
    # print("z_max", z_max)
    chem_sym = struct.get_chemical_symbols()
    Nat = len(chem_sym)
    flag = 0
    while flag == 0:
        indx1 = randrange(Nat)
        indx2 = randrange(Nat)
        z1 = struct.get_positions()[indx1][2]
        z2 = struct.get_positions()[indx2][2]
        # print(z1,z2)
        if indx2 != indx1 and chem_sym[indx1] != chem_sym[indx2] and z1 > z_max and z2 > z_max:
            # print(indx1,indx2)
            chem_sym_new = swap_positions(chem_sym, indx1, indx2)
            struct.set_chemical_symbols(chem_sym_new)
            flag = 1
    ind = Individual(struct, -1)
    return ind


# @profile
def new_by_switch_layers(ind_in: Individual):
    # returns individual made by switch of two random layers
    struct = ind_in.get_relaxed_structure().copy()
    write(filename="POSCAR.in.in", images=struct, format="espresso-in")
    layers = get_layers(struct)
    nlay = len(layers)
    if nlay > 1:
        # random indexes
        indx1 = randrange(nlay)
        flag = 0
        while flag == 0:
            indx2 = randrange(nlay)
            if indx2 != indx1:
                flag = 1
        layer1 = layers[indx1]
        layer2 = layers[indx2]
        print("switching layers " + str(indx1 + 1) + " and " + str(indx2 + 1))
        pos = struct.get_positions()
        chem = struct.get_chemical_symbols()
        cell = struct.get_cell()
        new_pos = []
        new_chem = []
        # 1 - setting all other atoms
        n = len(chem)
        for i in range(n):
            if i not in layer1 and i not in layer2:
                # print(i)
                new_pos.append(pos[i])
                new_chem.append(chem[i])
        # 2 - calculate shift between layers
        z1 = []
        z2 = []
        for i in layer1:
            z = pos[i][2]
            z1.append(z)
        for i in layer2:
            z = pos[i][2]
            z2.append(z)
        min1 = max(z1)
        min2 = max(z2)
        zshift = min1 - min2
        # 3 - setting layer1 and layer2 atoms with shift
        for i in range(n):
            if i in layer1:
                # print(i)
                if min1 > min2:
                    pos[i][2] = pos[i][2] + zshift
                if min1 < min2:
                    pos[i][2] = pos[i][2] - zshift
                new_pos.append(pos[i])
                new_chem.append(chem[i])
            if i in layer2:
                # print(i)
                if min1 > min2:
                    pos[i][2] = pos[i][2] - zshift
                if min1 < min2:
                    pos[i][2] = pos[i][2] + zshift
                new_pos.append(pos[i])
                new_chem.append(chem[i])
        st = ""
        for chem in new_chem:
            st = st + chem
        new_struct = Atoms(st)
        new_struct.set_cell(cell)
        new_pos2 = np.asarray(new_pos)
        new_struct.set_positions(new_pos2)
        ind = Individual(new_struct)
        write(filename="POSCAR.out.in", images=new_struct, format="espresso-in")
    else:
        ind = Individual(struct)
        print("warning: only one layer, no change in structure!")
    return ind


# @profile
def new_by_random_coordinates(ind_in: Individual):
    # returns individual made by changing all coordinates randomly (cell is intact)
    struct = ind_in.get_relaxed_structure().copy()
    Nat = len(struct.get_chemical_symbols())
    flag = 0
    while flag == 0:
        positions = []
        #
        for i in range(Nat):
            pos = [uniform(0., 0.90), uniform(0., 0.90), uniform(-0.1, 0.1)]
            positions.append(pos)
        #    
        struct.set_scaled_positions(positions)
        positions_ang = struct.get_positions()
        # print(positions_ang)
        flag = check_dist(Nat, positions_ang, dmin=1.4)
        # print(flag)
    # print(struct)
    ind = Individual(struct, -1)
    return ind


def new_by_shift_coordinates(ind_in: Individual):
    # returns individual made by random change of x,y coordinates by random delta
    struct = ind_in.get_relaxed_structure().copy()
    Nat = len(struct.get_chemical_symbols())
    positions = struct.get_scaled_positions()
    delx = uniform(0., 0.2)
    dely = uniform(0., 0.2)
    indx1 = randrange(Nat)
    indx2 = indx1
    while indx2 == indx1:
        indx2 = randrange(Nat)
    positions[indx1][0] = positions[indx1][0] + delx
    positions[indx1][1] = positions[indx1][1] + dely
    positions[indx2][0] = positions[indx2][0] - delx
    positions[indx2][1] = positions[indx2][1] - dely
    struct.set_scaled_positions(positions)
    ind = Individual(struct, -1)
    return ind


def new_by_shift_coordinate(ind_in: Individual):
    # returns individual made by random change of one coordinate
    struct = ind_in.get_relaxed_structure().copy()
    Nat = len(struct.get_chemical_symbols())
    positions = struct.get_scaled_positions()
    delta = uniform(-0.1, 0.1)
    indx_at = randrange(Nat)
    indx_coord = randrange(3)
    positions[indx_at][indx_coord] = positions[indx_at][indx_coord] + delta
    struct.set_scaled_positions(positions)
    ind = Individual(struct)
    return ind



