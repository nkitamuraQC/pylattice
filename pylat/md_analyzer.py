import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

class MD_Analyzer:
    def __init__(self, prefix, natoms=2):
        self.prefix = prefix
        self.input_name = self.prefix + ".in"
        self.log_name = self.prefix + ".out"
        self.natoms = natoms
        self.energies = None
        self.stresses_au = None
        self.stresses_kbar = None
        self.kinetics = None
        self.temps = None
        self.forces = None
        self.volumes = None
        
        rf = open(self.log_name, "r")
        self.lines = rf.readlines()
        rf.close()


    def get_energies(self):
        self.energies = []
        for i, line in enumerate(self.lines):
            if "!    total energy" in line:
                target_line = line
                self.energies.append(float(target_line.split("=")[1].split("Ry")[0]))
        print(len(self.energies))
        return
    

    def get_stress_tensors(self):
        self.stresses_au = []
        self.stresses_kbar = []
        for i, line in enumerate(self.lines):
            if "total   stress  (Ry/bohr**3)" in line:
                tmp_au = []
                tmp_kbar = []
                for j in range(3):
                    tmp_au.append(list(map(float, self.lines[i+j+1].split()))[:3])
                    tmp_kbar.append(list(map(float, self.lines[i+j+1].split()))[3:])
                self.stresses_au.append(np.array(tmp_au))
                self.stresses_kbar.append(np.array(tmp_kbar))
        print(len(self.stresses_au))       
        return
    

    def get_kinetic_energies(self):
        self.kinetics = []
        for i, line in enumerate(self.lines):
            if "kinetic energy (Ekin)" in line:
                target_line = line
                self.kinetics.append(float(target_line.split("=")[1].split("Ry")[0]))
        print(len(self.kinetics))     
        return
    

    def get_temperatures(self):
        self.temps = []
        for i, line in enumerate(self.lines):
            if "temperature           =" in line:
                target_line = line
                self.temps.append(float(target_line.split("=")[1].split("K")[0]))
        print(len(self.temps))     
        return
    

    def get_forces(self):
        self.forces = []
        for i, line in enumerate(self.lines):
            if "Forces acting on atoms (cartesian axes, Ry/au):" in line:
                count = 2
                tmp = []
                while True:
                    if "force" in self.lines[i+count]:
                        target_line = self.lines[i+count]
                        tmp.append(list(map(float, target_line.split("=")[1].split())))
                    else:
                        self.forces.append(np.array(tmp))
                        break
                    count += 1
        print(len(self.forces))     
        return
    

    def get_coords(self):
        self.coords = []
        self.elems = []
        for i, line in enumerate(self.lines):
            if "ATOMIC_POSITIONS (angstrom)" in line:
                tmp = []
                for j in range(self.natoms):
                    data = list(map(float, self.lines[i+j+1].split()[1:]))
                tmp.append(data)
                self.coords.append(np.array(tmp))
        print(len(self.coords))   
        return
    

    def get_volumes(self):
        raise NotImplementedError
    

    def parse(self):
        self.get_energies()
        self.get_stress_tensors()
        self.get_kinetic_energies()
        self.get_temperatures()
        self.get_forces()
        return
    

    def to_csv(self, csv_name="save.csv"):
        df = pd.DataFrame({
            "energies": self.energies,
            "temperatures": self.temps,
            "kinetics": self.kinetics,
            "stresses_au": self.stresses_au,
            "stresses_kbar": self.stresses_kbar,
            "forces": self.forces
        })

        df.to_csv(csv_name)
        return
    

    def plot(self, target=""):
        if "stress" in target:
            raise NotImplementedError
        elif "forces" in target:
            raise NotImplementedError
        else:
            steps = [i for i, _ in enumerate(self.energies)]
            plt.scatter(steps, self[target])
            plt.savefig(target+".png")
        return
    

    def plot_stress(self, index=None, typ="au"):
        if isinstance(index, list):
            a, b = index
            plot_data = []
            if typ == "au":
                for i, val in enumerate(self.stresses_au):
                    plot_data.append(val[a, b])
            if typ == "kbar":
                for i, val in enumerate(self.stresses_kbar):
                    plot_data.append(val[a, b])
            steps = [i for i, _ in enumerate(self.energies)]
            plt.scatter(steps, plot_data)
            plt.savefig(f"stress_{a}_{b}"+".png")
        else:
            raise NotImplementedError
        return
    

    def plot_forces(self, atom_index=0, index=0):
        plot_data = []
        for i, val in enumerate(self.forces):
            plot_data.append(val[atom_index, index])
        steps = [i for i in range(self.energies)]
        plt.scatter(steps, plot_data)
        plt.savefig(f"stress_{atom_index}_{index}"+".png")
        return

    

    def __len__(self):
        return len(self.__dict__)

    def __repr__(self):
        return str(self.__dict__)

    def __str__(self):
        return str(self.__dict__)

    def __iter__(self):
        return self.__dict__.iteritems()

    def __getitem__(self, key):
        return self.__dict__[key]

    def __setitem__(self, key, value):
        self.__dict__[key] = value
    


    
if __name__ == "__main__":
    prefix = "calculation_try_0"
    mda = MD_Analyzer(prefix)
    mda.parse()
    mda.plot(target="energies")
    mda.plot_stress(index=[0, 0])


