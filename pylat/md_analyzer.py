import matplotlib
from matplotlib import pyplot as plt

class MD_Analyzer:
    def __init__(self, prefix):
        self.prefix = prefix
        self.input_name = self.prefix + ".in"
        self.log_name = self.prefix + ".out"
        self.energies = None
        self.stress = None
        self.kinetic = None
        self.temp = None
        self.force = None
        self.volume = None


    def get_energies(self):
        return
    

    def get_stress_tensors(self):
        return
    

    def get_kinetic_energies(self):
        return
    

    def get_temperatures(self):
        return
    

    def get_forces(self):
        return
    

    def get_volumes(self):
        return
    

    def analyse(self):
        return
    

    def plot(self, target=""):
        return
    

    def plot_stress(self, index=None):
        return
    

    def plot_forces(self, atom_index=0):
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
    


    


