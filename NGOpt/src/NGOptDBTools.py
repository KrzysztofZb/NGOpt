# NGOpt Copyright(C) 2020 Krzysztof Zberecki
# This program comes with ABSOLUTELY NO WARRANTY; for details type 'show w'.
# This is free software, and you are welcome to redistribute it under certain
# conditions; type 'show c' for details.

# tools for data extraction from various databases
# now works with c2db database

import os
import numpy as np
import pickle
import sqlite3

import ase.db
from ase import Atoms
from ase.io import write, read

from NGOptIndividual import Individual


def get_c2db_prototypes(dbfile, prototype, stability=False, smagstate=False, sgap=0):
    # returns list of Individuals from c2db
    # Connect to the database - ase interface
    db = ase.db.connect(dbfile)
    # Connect to the database - sqlite3 interface
    conn = sqlite3.connect(dbfile)
    cur = conn.cursor()
    cur.execute('SELECT SQLITE_VERSION()')
    data = cur.fetchone()
    print("SQLite version: " + str(format(data)))

    # nstruct = 0
    # magstate = 'FM'
    dyn_key = 'dynamic_stability_level'
    th_key = 'thermodynamic_stability_level'
    # gap_key = 'gap'

    # rows = db.select('magstate=FM')
    if prototype != "":
        query = 'class=' + prototype
        # print(query)
        rows = db.select(query)
    if prototype == "":
        # rows = db.select('calculator=gpaw')
        rows = db.select('gap>-1')
    print(rows)

    individuals = []

    for row in rows:
        id = row.get("id")
        e_tot = row.get("energy")
        label = row.get('formula')
        gap = row.get('gap')
        #print(gap)
        magstate = row.get('magstate')
        magmom = row.get('magmom')
        positions = row.get('positions')
        cell = row.get('cell')
        hform = row.get('hform')
        e_tot = row.get('energy')

        # Stability info - not implemented is ase interface
        cur.execute("SELECT * FROM number_key_values WHERE id=:id AND key=:dyn_key", {"id": id, "dyn_key": dyn_key})
        data = cur.fetchone()
        if data is not None:
            stability_dyn = data[1]
        else:
            stability_dyn = 0.
        cur.execute("SELECT * FROM number_key_values WHERE id=:id AND key=:th_key", {"id": id, "th_key": th_key})
        data = cur.fetchone()
        if data is not None:
            stability_th = data[1]
        else:
            stability_th = 0.

        if stability == False:
            struct_ase = Atoms(label, positions=positions, cell=cell, pbc=np.array([True, True, True], dtype=bool))
            Nat = len(struct_ase.get_chemical_symbols())
            ind = Individual(struct_ase)
            ind.set_e_tot(e_tot / float(Nat))
            ind.set_mmtot(magmom)
            ind.set_hform(hform)
            individuals.append(ind)

        if stability == True and smagstate == False and sgap==0:
            if stability_th == 3.0:
                struct_ase = Atoms(label, positions=positions, cell=cell, pbc=np.array([True, True, True], dtype=bool))
                Nat = len(struct_ase.get_chemical_symbols())
                ind = Individual(struct_ase)
                ind.set_e_tot(e_tot / float(Nat))
                ind.set_mmtot(magmom)
                ind.set_hform(hform)
                individuals.append(ind)
        
        if stability == True and smagstate == False and sgap>0:
            if stability_th == 3.0 and gap>0.:
                struct_ase = Atoms(label, positions=positions, cell=cell, pbc=np.array([True, True, True], dtype=bool))
                Nat = len(struct_ase.get_chemical_symbols())
                ind = Individual(struct_ase)
                ind.set_e_tot(e_tot / float(Nat))
                ind.set_mmtot(magmom)
                ind.set_hform(hform)
                individuals.append(ind)
        
        if stability == True and smagstate == True and smagstate == "FM":
            if stability_dyn == 3.0 and stability_th == 3.0 and abs(float(magmom)) > 0.1:
                struct_ase = Atoms(label, positions=positions, cell=cell, pbc=np.array([True, True, True], dtype=bool))
                Nat = len(struct_ase.get_chemical_symbols())
                ind = Individual(struct_ase)
                ind.set_e_tot(e_tot / float(Nat))
                ind.set_mmtot(magmom)
                ind.set_hform(hform)
                individuals.append(ind)

    return individuals


def create_populations_db(dbfile):
    os.system("rm " + dbfile)
    conn = sqlite3.connect(dbfile)
    cur = conn.cursor()

    cur.execute('''CREATE TABLE individuals
             ([id]	INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
	          [iid]	TEXT,
	          [calculator]	TEXT,
	          [init_structure]	BLOB,
	          [relaxed_structure]	BLOB,
              [id_symm_group]    INTEGER,
	          [e_tot]	REAL,
	          [mmtot]	REAL,
	          [hform]	REAL,
	          [origin]	TEXT,
	          [ngen]	INTEGER
             )''')

    conn.close()


def report_populations_db(dbfile):
    conn = sqlite3.connect(dbfile)
    cur = conn.cursor()

    cur.execute("SELECT * FROM individuals ORDER BY e_tot ASC")
    rows = cur.fetchall()

    nstr = 1
    for row in rows:
        print("id: " + str(row[0]))
        print("e_tot: " + str(row[6]))
        print("origin: " + str(row[9]))
        print("ngen: " + str(row[10]))
        if nstr == 1:
            struct = pickle.loads(row[4])
            write(filename='best.in', images=struct, format='espresso-in')
        if nstr == 2:
            struct = pickle.loads(row[4])
            write(filename='best_second.in', images=struct, format='espresso-in')
        if nstr == 3:
            struct = pickle.loads(row[4])
            write(filename='best_third.in', images=struct, format='espresso-in')
        nstr = nstr + 1
        print("----------------------")


def get_population_db(dbfile, ngen):
    population = []
    conn = sqlite3.connect(dbfile)
    cur = conn.cursor()

    cur.execute("SELECT * FROM individuals WHERE ngen=?", [ngen])
    rows = cur.fetchall()

    for row in rows:
        struct = pickle.loads(row[4])
        ind = Individual(struct)
        population.append(ind)

    return population


def get_generations_number(dbfile):
    conn = sqlite3.connect(dbfile)
    cur = conn.cursor()

    cur.execute("SELECT MAX(ngen) FROM individuals ")
    rows = cur.fetchall()

    return rows[0][0]










