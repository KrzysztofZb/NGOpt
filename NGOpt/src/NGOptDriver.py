# NGOpt Copyright(C) 2020 Krzysztof Zberecki
# This program comes with ABSOLUTELY NO WARRANTY; for details type 'show w'.
# This is free software, and you are welcome to redistribute it under certain
# conditions; type 'show c' for details.

# TODO: add 'stop' condition

# main driver for the genetic algorithm

from NGOptOperations import *
from NGOptConfig import *
from NGOptDBTools import *


# from memory_profiler import profile


def get_best_individual(population):
    e_tot_min = 0.
    for ind in population:
        try:
            e_tot = ind.e_tot
            struct = ind.get_relaxed_structure()
        except:
            e_tot = 1000.
        if e_tot < e_tot_min:
            struct_best = struct
            e_tot_min = e_tot
    print("e_tot_min: " + str(e_tot_min))
    individual = Individual(struct_best)
    return individual


class Driver:

    def __init__(self, Npop, chem_sym, stoichio, model_file, calculator):
        self.Npop = Npop  # starting no. of individuals
        self.chem_sym = chem_sym  # chemical symbols - table
        self.stoichio = stoichio  # stoichiometry - table
        self.model_file = model_file
        self.chem_form = ""  # chemical formula
        self.calculator = calculator
        self.amin = 0.
        self.amax = 0.
        self.dmin = DMIN
        self.exec = ""
        self.dbfile = "populations.db"

        if self.calculator == "VASP":
            self.exec = VASP_EXEC
        if self.calculator == "siesta":
            self.exec = SIESTA_EXEC

        Nat = len(chem_sym)
        for i in range(Nat):
            self.chem_form = self.chem_form + self.chem_sym[i] + str(self.stoichio[i])

        outfile = open("Individuals", "w")
        outfile.close()

        create_populations_db(self.dbfile)

    def init_population_random(self, model=1):
        # inits random population with model verification, returns list of individuals
        print("generating initial population randomly ...")
        population = []
        i = 1
        while i <= self.Npop:
            if model == 0:
                struct = random_structure(self.chem_form, self.amin, self.amax, self.dmin)
            if model == 1:
                struct = random_structure_model(self.chem_form, self.amin, self.amax, self.dmin, self.model_file)
            individual = Individual(struct)
            individual.origin = "new ind. random from scratch"
            population.append(individual)
            i = i + 1
        return population

    def init_population_random_group(self):
        # inits pyxtal random population with model verification, returns list of individuals
        print("generating initial population randomly but using group theory ...")
        population = []
        i = 1
        while i <= self.Npop:
            struct = random_structure_group(self.chem_sym, self.stoichio, thickness=3.0, tol_factor=0.9,
                                            model_file=self.model_file, dmin=self.dmin, Natt=100)
            individual = Individual(struct)
            individual.origin = "new ind. random from scratch"
            population.append(individual)
            i = i + 1
        return population

    def init_population_random_on_substrate(self):
        print("generating initial population on substrate randomly ...")
        population = []
        i = 1
        while i <= self.Npop:
            struct = random_structure_on_substrate(self.chem_form, self.amin, self.amax, self.dmin, self.model_file)
            individual = Individual(struct)
            individual.origin = "new ind. random on substrate from scratch"
            population.append(individual)
            i = i + 1
        return population

    def init_population_from_file(self, dbfile):
        # reads population from db file, returns list of individuals
        ngen = get_generations_number(dbfile)
        population = get_population_db(dbfile, ngen)
        return population

    # @profile
    def evaluate_population(self, population, postfix):
        # evaluates population using self.calculator
        print("")
        print("evaluating population ...")
        print("")
        # new_population = [] # new - means evaluated...
        iot = IOTools(calculator=self.calculator)
        outfile = open("Individuals", "a")
        outfile.write("data of pop: " + postfix + " \n")
        i = 1
        for ind in population:  # main population loop
            struct = ind.get_init_structure()
            # generation of input
            iot.prepare_input(struct)
            print("evaluating ind. no. " + str(i) + " ...")
            # start calculation
            os.system(self.exec)
            # handle output
            if self.calculator == "VASP":
                structfile = "CONTCAR"
                r_outfile = "OSZICAR"
            if self.calculator == "siesta":
                structfile = "siesta.XV"
                r_outfile = "out.tmp"
            if self.calculator == "openmx":
                structfile = "openmx.dat#"
                r_outfile = "openmx.out"
            try:
                struct_relaxed = iot.read_structure(structfile)
                symbols = struct_relaxed.get_chemical_symbols()
                nat = len(symbols)
                e_tot = iot.read_tot_energy(r_outfile) / float(nat)
                ind.set_relaxed_structure(struct_relaxed)
                ind.set_e_tot(e_tot)
                info = "individual no. " + str(i) + ": e_tot/atom: " + str(e_tot) + " \n"
                print(info)
                outfile.write(info)
            except:
                print("warning, individual no. " + str(i) + " failed!")
                ind.set_e_tot(0.)
            # if self.calculator == "VASP":
            #     os.system("mv CONTCAR POSCAR_" + postfix + "_" + str(i))
            #     os.system("mv OSZICAR OSZICAR_" + postfix + "_" + str(i))
            #     os.system("mv OUTCAR OUTCAR_" + postfix + "_" + str(i))
            if self.calculator == "siesta":
                os.system("mv siesta.XV siesta.XV." + str(i) + "." + postfix)
                os.system("mv out.tmp out." + str(i) + "." + postfix)
            if self.calculator == "openmx":
                os.system("mv openmx.dat# openmx.dat#." + str(i) + "." + postfix)
                os.system("mv openmx.out openmx.out." + str(i) + "." + postfix)
            i = i + 1

    # @profile
    def generate_new_population_alg1(self, population):
        # generates new population using algorithm no. 1
        print("")
        print("generating new population with softmutation ...")
        print("")
        new_population = []
        population.sort()
        Npop = self.Npop
        adapt = AseAtomsAdaptor()
        model = MEGNetModel.from_file(self.model_file)
        for i in range(int(Npop / 4)):  # 1/4 of population by atoms switch
            indx = randrange(Npop / 4)  # +1 # using only best 25% of population
            print("-------------------------------")
            nlay = len(get_layers(population[indx].get_relaxed_structure()))
            if nlay > 1:
                print("new ind. from no. " + str(indx) + " by layers switch")
                individual = new_by_switch_layers(population[indx])
                individual.origin = "new ind. from no. " + str(indx) + " by layers switch"
            else:
                print("new ind. from no. " + str(indx) + " by atoms switch")
                individual = new_by_switch_atoms(population[indx])
                individual.origin = "new ind. from no. " + str(indx) + " by atoms switch"
            new_population.append(individual)
        for i in range(int(Npop / 4)):  # 1/4 of population by kind of softmutation
            indx = randrange(Npop / 4)  # +1 # using only best 25% of population
            print("-------------------------------")
            print("new ind. from no. " + str(indx) + " by softmutation")
            ind_in = population[indx]
            structure_ase = ind_in.get_relaxed_structure().copy()
            structure_pymatgen = adapt.get_structure(structure_ase)
            e_tot_in = model.predict_structure(structure_pymatgen)
            e_tot_min = e_tot_in
            flag = 0
            for i in range(100):
                individual = new_by_shift_coordinate(population[indx])
                structure_ase = individual.get_init_structure().copy()
                structure_pymatgen = adapt.get_structure(structure_ase)
                try:
                    e_tot = model.predict_structure(structure_pymatgen)
                except:
                    e_tot = 0.
                    print("isolated molecule exception handled")
                if e_tot < e_tot_min:
                    flag = 1
                    e_tot_min = e_tot
                    ind_out = individual
                    ind_out.origin = "new ind. from no. " + str(indx) + " by softmutation"
            if flag == 1:
                new_population.append(ind_out)
                print("softmutate energy gain: " + str(e_tot_min - e_tot_in))
            if flag == 0:
                new_population.append(ind_in)
        for i in range(int((Npop / 2) - 1)):
            # 1/2 of population (minus one) by new random structures, pyxtal+model
            print("-------------------------------")
            print("new ind. random from scratch")
            struct = random_structure_group(self.chem_sym, self.stoichio, thickness=3.0, tol_factor=0.9,
                                            model_file=self.model_file, dmin=self.dmin, Natt=RANDOM_ATTEMPTS)
            individual = Individual(struct)
            individual.origin = "new ind. random from scratch"
            new_population.append(individual)
        best_ind = get_best_individual(population)  # last one is the best in the current population         
        best_ind.origin = "kept best"
        new_population.append(best_ind)
        del model
        return new_population

    def generate_new_population_alg2_substrate(self, population):
        # generates new population using algorithm no. 2
        print("")
        print("generating new population with softmutation ...")
        print("")
        new_population = []
        population.sort()
        Npop = self.Npop
        adapt = AseAtomsAdaptor()
        model = MEGNetModel.from_file(self.model_file)
        for i in range(int(Npop / 4)):  # 1/4 of population by atoms switch
            indx = randrange(Npop / 4)  # +1 # using only best 25% of population
            print("-------------------------------")
            print("new ind. from no. " + str(indx) + " by atoms switch")
            individual = new_by_switch_atoms_on_substrate(population[indx])
            individual.origin = "new ind. from no. " + str(indx) + " by atoms switch"
            new_population.append(individual)
        for i in range(int(Npop / 4)):  # 1/4 of population by kind of softmutation
            indx = randrange(Npop / 4)  # +1 # using only best 25% of population
            print("-------------------------------")
            print("new ind. from no. " + str(indx) + " by softmutation")
            ind_in = population[indx]
            structure_ase = ind_in.get_relaxed_structure().copy()
            structure_pymatgen = adapt.get_structure(structure_ase)
            e_tot_in = model.predict_structure(structure_pymatgen)
            e_tot_min = e_tot_in
            flag = 0
            for i in range(100):
                individual = new_by_shift_coordinate(population[indx])
                structure_ase = individual.get_init_structure().copy()
                structure_pymatgen = adapt.get_structure(structure_ase)
                try:
                    e_tot = model.predict_structure(structure_pymatgen)
                except:
                    e_tot = 0.
                    print("isolated molecule exception handled")
                if e_tot < e_tot_min:
                    flag = 1
                    e_tot_min = e_tot
                    ind_out = individual
            ind_out.origin = "new ind. from no. " + str(indx) + " by softmutation"
            if flag == 1:
                new_population.append(ind_out)
                print("softmutate energy gain: " + str(e_tot_min - e_tot_in))
            if flag == 0:
                new_population.append(ind_in)
        for i in range(int((Npop / 2) - 1)):
            # 1/2 of population (minus one) by new random structures, pyxtal+model
            # indx = randrange(Npop/2)+1
            print("-------------------------------")
            print("new ind. random from scratch")
            struct = random_structure_group(self.chem_sym, self.stoichio, thickness=3.0, tol_factor=0.9,
                                            model_file=self.model_file, dmin=self.dmin, Natt=RANDOM_ATTEMPTS)
            individual = Individual(struct)
            individual.origin = "new ind. random from scratch"
            new_population.append(individual)
        best_ind = get_best_individual(population)  # last one is the best in the current population         
        best_ind.origin = "kept best"
        new_population.append(best_ind)
        return new_population

    def population_report(self, population, ngen):
        print("---------------------")
        print("population report:")
        i = 1
        for ind in population:
            e_tot = ind.e_tot
            origin = ind.origin
            try:
                struct = ind.get_relaxed_structure()
                info = "ind. no. " + str(i) + " : e_tot=" + str(e_tot) + " origin: " + origin
            except:
                info = "ind. no. " + str(i) + " failed..."
            i = i + 1
            ind.save_to_db(self.dbfile, ngen)
            print(info)
        print("---------------------")
