
# i/o tools

import numpy as np
import os
import json

from ase import Atoms
from ase.io import write
from ase.io.vasp import *

from ase.calculators.siesta.import_functions import xv_to_atoms
from ase.calculators.siesta import Siesta

from NGOptConfig import *

class IOTools:

    def __init__(self, calculator):
        self.calculator = calculator  # "VASP" or "siesta" or "openmx" for now
        if self.calculator == "VASP":
            self.vasp_pseudo_path = VASP_PSEUDO_PATH
        if self.calculator == "openmx":
            self.openmx_bases_file = "openmx_bases.json"  # these files need to be present in current dir
            self.openmx_pseudo_file = "pseudos.json"
            self.halfmom = 0.0

    def read_structure(self, structfile):
        # returns ase structure from proper structure file
        structure = ""
        if self.calculator == "VASP":
            structure = read_vasp(structfile)  # structure - instance of ASE.Atom class
        if self.calculator == "siesta":
            structure = xv_to_atoms(structfile)
        if self.calculator == "openmx":  # hashfile of openmx.dat (openmx.dat#)
            pos, species, cell = read_openmx_hashfile(structfile)
            structure = Atoms(species)
            structure.set_positions(pos)
            structure.set_cell(cell)
        return structure

    def write_structure(self, struct, filename="POSCAR", format="vasp"):
        # writes ase struct (Atoms) in one of formats implemented in ase lib
        write(filename=filename, images=struct, format=format)

    def prepare_input(self, struct):
        if self.calculator == "VASP":
            # POTCAR and POSCAR file generation, INCAR and KPOINTS must be provided in current dir
            write(filename="POSCAR", images=struct, format="vasp")
            if os.path.isfile('POTCAR'):
                os.system("rm POTCAR")
            tmp = open("POSCAR", "r")
            first_line = tmp.readline()
            line_sp = first_line.split()
            cat = "cat "
            # path from NGOptConfig
            for s in line_sp:
                cat = cat + self.vasp_pseudo_path+"/POTCAR_" + s + " "
            cat = cat + " > POTCAR"
            tmp.close()
            os.system(cat)
        if self.calculator == "siesta":
            # head.fdf file generation, tail.fdf must be provided in current dir
            calc = Siesta(label='siesta', xc='PBE', basis_set='SZP', )
            f = open("head.fdf", "w")
            calc._write_species(f, struct)
            calc._write_structure(f, struct)
            f.close()
            os.system("cat head.fdf tail.fdf > siesta.fdf")
        if self.calculator == "openmx":
            # center of file generation, head and tail must be provided in current dir
            label = "openmx"
            input = open("openmx.in.dat", "w")
            #
            input.write("\n")
            input.write("System.Name                      " + label + "\n")
            input.write("\n")
            input.write("# \n")
            input.write("# Definition of Atomic Species \n")
            input.write("# \n")
            input.write("\n")
            # print openmx species
            l = len(struct.get_chemical_symbols())
            symbols = struct.get_chemical_symbols()
            symbols_dist = set(symbols)
            input.write("Species.Number       " + str(len(symbols_dist)) + "\n")
            input.write("<Definition.of.Atomic.Species \n")
            for symbol in symbols_dist:
                basis = get_basis_info_openmx(self.openmx_bases_file, symbol)
                pseudo = get_pseudos_openmx(self.openmx_pseudo_file, symbol)
                size = str(basis[2])
                input.write(str(symbol) + "  " + str(basis[0] + "-" + size) + "  " + str(pseudo) + " \n")
            input.write("Definition.of.Atomic.Species> \n")
            input.write("\n")
            input.write("# \n")
            input.write("# Atoms \n")
            input.write("# \n")
            input.write("\n")
            # print openmx positions
            positions = struct.get_positions()
            input.write("Atoms.Number      " + str(l) + "\n")
            input.write("Atoms.SpeciesAndCoordinates.Unit   Ang # Ang|AU \n")
            input.write("<Atoms.SpeciesAndCoordinates \n")
            for i in range(l):
                symbol = symbols[i]
                basis = get_basis_info_openmx(self.openmx_bases_file, symbol)
                valence = float(basis[1]) / 2.
                if basis[2] == "s2p2d2":
                    valence1 = valence + self.halfmom
                    valence2 = valence - self.halfmom
                else:
                    valence1 = valence2 = valence
                    # valence1 = valence+self.halfmom
                    # valence2 = valence-self.halfmom
                input.write(
                    str(i + 1) + "    " + str(symbols[i]) + "    " + str(format(positions[i][0], '.12f')) + "   " + str(
                        format(positions[i][1], '.12f')) + "   " + str(format(positions[i][2], '.12f')) + "    " + str(
                        valence1) + "    " + str(valence2) + "\n")
            input.write("Atoms.SpeciesAndCoordinates> \n")
            # print openmx unit cell
            input.write("Atoms.UnitVectors.Unit             Ang # Ang|AU \n")  # do poprawy!
            input.write("<Atoms.UnitVectors \n")
            cell = struct.get_cell()
            input.write(str(format(cell[0][0], '.12f')) + "    " + str(format(cell[0][1], '.12f')) + "    " + str(
                format(cell[0][2], '.12f')) + "\n")
            input.write(str(format(cell[1][0], '.12f')) + "    " + str(format(cell[1][1], '.12f')) + "    " + str(
                format(cell[1][2], '.12f')) + "\n")
            input.write(str(format(cell[2][0], '.12f')) + "    " + str(format(cell[2][1], '.12f')) + "    " + str(
                format(cell[2][2], '.12f')) + "\n")
            input.write("Atoms.UnitVectors> \n")
            input.close()
            os.system("cat head openmx.in.dat tail > openmx.dat")

    def read_tot_energy(self, out):  # returns e_tot from proper output file
        outfile = open(out, "r")
        e_tot = 0.
        # VASP case
        if self.calculator == "VASP":  # in VASP case - OSZICAR file
            for line in outfile:
                test = line.find("F=")
                if test > -1:
                    i = 0
                    line_sp = line.split()
                    for l in line_sp:
                        if l == "F=":
                            e_tot = float(line_sp[i + 1])
                        i = i + 1
        # siesta case
        if self.calculator == "siesta":
            # in siesta case - proper out text file
            outfile = open(out, "r")
            e_tot = 0.
            for line in outfile:
                test = line.find("siesta:         Total =")
                if test > -1:
                    line_sp = line.split()
                    e_tot = float(line_sp[3])
        # openmx case
        if self.calculator == "openmx":
            # in openmx case - label.out text file
            outfile = open(out, "r")
            e_tot = 0.
            for line in outfile:
                test = line.find("Utot.")
                if test > -1:
                    line_sp = line.split()
                    e_tot = float(line_sp[1])
            outfile.close()
        return e_tot

    def read_tot_mm(self, out):
        # returns mmtot from proper out file, VASP case only implemented
        outfile = open(out, "r")
        mmtot = 0.
        # VASP case
        if self.calculator == "VASP":  # in VASP case - OSZICAR file
            for line in outfile:
                test = line.find("mag=")
                if test > -1:
                    i = 0
                    line_sp = line.split()
                    for l in line_sp:
                        if l == "mag=":
                            mmtot = float(line_sp[i + 1])
                        i = i + 1
        return mmtot

    def get_best_individual(self, out="Individuals", calculator="VASP",  mode="2D"):
        # returns e_tot/atom of best individual, mode="2D" or "all", uses Individuals file
        # deprecated
        outfile = open(out, "r")
        e_tot_min = 0.
        indx_min = 0
        indx_pop = 1
        indx_pop_min = 1
        for line in outfile:
            test = line.find("data of pop:")
            if test > -1:
                line_sp = line.split()
                pop = str(line_sp[3])
                pop = pop.replace('POP', '')
                indx_pop = int(pop)
                # print(indx_pop)
            test = line.find("individual ")
            if test > -1:
                line_sp = line.split()
                e_tot = float(line_sp[4])
                tmp = str(line_sp[2])
                tmp = tmp.replace(':', '')
                indx = int(tmp)
                if mode == "all":
                    if e_tot < e_tot_min:
                        e_tot_min = e_tot
                        indx_min = int(tmp)
                        indx_pop_min = indx_pop
                if mode == "2D":
                    poscar = "POSCAR_POP" + str(indx_pop) + "_" + str(indx)
                    struct = read_vasp(poscar)
                    cell = struct.get_cell()
                    z = cell[2][2]
                    if e_tot < e_tot_min and z > 8. and e_tot > -25.: # z>8. -> 2D, e_tot/at > -25. to avoid not converged results
                        e_tot_min = e_tot
                        indx_min = int(tmp)
                        indx_pop_min = indx_pop
        # print(indx_pop_min,indx_min)

        if calculator == "VASP":
            poscar = "POSCAR_POP" + str(indx_pop_min) + "_" + str(indx_min)
            struct = read_vasp(poscar)
            write(filename='POSCAR.best', images=struct, format='vasp')
            write(filename='POSCAR.best.in', images=struct, format='espresso-in')

        if calculator == "siesta":
            outfile = "siesta.XV." + str(indx_min) + ".POP" + str(indx_pop_min)
            struct = xv_to_atoms(outfile)
            write(filename='POSCAR.best', images=struct, format='vasp')
            write(filename='POSCAR.best.in', images=struct, format='espresso-in')

        if calculator == "openmx":
            # raise ValueError('not yet implemented for openmx case')
            outfile = "openmx.dat#." + str(indx_min) + ".POP" + str(indx_pop_min)
            pos, species, cell = read_openmx_hashfile(outfile)
            struct = Atoms(species)
            struct.set_positions(pos)
            struct.set_cell(cell)
            write(filename='POSCAR.best', images=struct, format='vasp')
            write(filename='POSCAR.best.in', images=struct, format='espresso-in')

        return e_tot_min

    def simple_raport(self, out="Individuals"):
        # short raport, uses Individuals file
        # deprecated
        outfile = open(out, "r")
        datafile = open("generations.dat", "w")
        indx_pop = 1
        for line in outfile:
            test = line.find("data of pop:")
            if test > -1:
                line_sp = line.split()
                pop = str(line_sp[3])
                pop = pop.replace('POP', '')
                indx_pop = int(pop)
                # print(indx_pop)
            test = line.find("individual ")
            if test > -1:
                line_sp = line.split()
                e_tot = float(line_sp[4])
                tmp = str(line_sp[2])
                tmp = tmp.replace(':', '')
                datafile.write(str(indx_pop) + "  " + str(e_tot) + "\n")
        datafile.close()


# helper functions for openmx input generation

def get_basis_info_openmx(basisfile, index):
    with open(basisfile, 'r') as f:
        bases = json.load(f)
    out = bases[index]
    return out


def get_pseudos_openmx(pseudofile, index):
    with open(pseudofile, 'r') as f:
        pseudos = json.load(f)
    out = pseudos[index]
    return out


def get_natoms_openmx(outfile):
    out = open(outfile, "r")
    natoms = 0
    for line in out:
        test = line.find("Atoms.Number")
        if test > -1:
            line_sp = line.split()
            natoms = int(line_sp[1])
    out.close()
    return natoms


def read_openmx_hashfile(hashfile):
    nat = get_natoms_openmx(hashfile)
    out = open(hashfile, "r")
    pos = np.zeros((nat, 3))
    species = []
    #
    #
    i = 0
    for line in out:
        i = i + 1
        test = line.find("<Atoms.SpeciesAndCoordinates")
        if test > -1:
            pline = i
    out.close()
    out = open(hashfile, "r").read()
    lines = out.split("\n")
    i = 0
    j = 0
    for line in lines:
        i = i + 1
        if i > pline and i < pline + nat + 1:
            line_sp = line.split()
            pos[j][0] = float(line_sp[2])
            pos[j][1] = float(line_sp[3])
            pos[j][2] = float(line_sp[4])
            species.append(str(line_sp[1]))
            j = j + 1
    # cell
    cell = np.zeros((3, 3))
    out = open(hashfile, "r")
    #
    i = 0
    for line in out:
        i = i + 1
        test = line.find("<Atoms.UnitVectors")
        if test > -1:
            pline = i
    out.close()
    out = open(hashfile, "r").read()
    lines = out.split("\n")
    i = 0
    j = 0
    for line in lines:
        i = i + 1
        if i > pline and i < pline + 3 + 1:
            line_sp = line.split()
            cell[j][0] = float(line_sp[0]) * 0.529177
            cell[j][1] = float(line_sp[1]) * 0.529177
            cell[j][2] = float(line_sp[2]) * 0.529177
            j = j + 1
    return pos, species, cell


# other helper functions
def get_symbols(chem_sym, stoich):
    Nat = len(chem_sym)
    chem_form =""
    for i in range(Nat):
        chem_form = chem_form + chem_sym[i] + str(stoich[i])

    atm = Atoms(chem_sym)
    symbols = atm.get_chemical_symbols()

    return symbols

