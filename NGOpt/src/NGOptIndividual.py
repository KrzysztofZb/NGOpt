# main class describing individual
# uses ase structures (i.e. ase class - Atoms)

# import numpy as np
# from random import randrange

import uuid
import pickle
import sqlite3 as sqlite


class Individual:

    def __init__(self, structure, id_symm_group=-1):
        self.id = uuid.uuid1()
        self.calculator = "VASP"  # default
        self.init_structure = structure  # ase Atoms class object
        self.relaxed_structure = None
        self.id_symm_group = id_symm_group  # id of symm. group
        self.e_tot = 0.0  # total energy/atom
        self.mmtot = 0.0  # total magnetic moment/cell
        self.hform = 0.0  # formation energy/atom
        self.origin = ""  # origin of an individual
        self.ngen = 0  # generation number

    def save_to_db(self, dbfile, ngen):
        # inserts Individual to database
        iid = str(self.id)
        calculator = self.calculator
        init_structure = pickle.dumps(self.init_structure)
        relaxed_structure = pickle.dumps(self.relaxed_structure)
        id_symm_group = self.id_symm_group
        e_tot = self.e_tot
        mmtot = self.mmtot
        hform = self.hform
        origin = self.origin
        ngen = ngen
        # insert data to db
        db = sqlite.connect(dbfile)
        cur = db.cursor()

        cur.execute('INSERT INTO individuals (iid, calculator, init_structure, relaxed_structure, id_symm_group, e_tot, mmtot, hform, origin, ngen) \
                     VALUES(?,?,?,?,?,?,?,?,?,?)', \
                    (iid, calculator, init_structure, relaxed_structure, id_symm_group, e_tot, mmtot, hform, origin, ngen,))
        db.commit()
        # db.close()

    def __lt__(self, other):
        return self.e_tot < other.e_tot

    def get_init_structure(self):
        return self.init_structure.copy()

    def get_relaxed_structure(self):
        return self.relaxed_structure.copy()

    def get_e_tot(self):
        return self.e_tot

    def set_relaxed_structure(self, structure):
        self.relaxed_structure = structure

    def set_id_symm_group(self, id_symm_group):
        self.id_symm_group = id_symm_group

    def set_e_tot(self, energy):
        self.e_tot = energy

    def set_mmtot(self, mmtot):
        self.mmtot = mmtot

    def set_hform(self, hform):
        self.hform = hform

    def set_origin(self, origin):
        self.origin = origin

    def set_ngen(self, ngen):
        self.ngen = ngen
