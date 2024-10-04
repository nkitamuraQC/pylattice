from ase.spacegroup import crystal

class SpaceGroupSupplement:
    def __init__(self, ciffile, xyzfile, lattfile):
        self.ciffile = ciffile
        self.xyzfile = xyzfile
        self.lattfile = lattfile

        self.xyz = None
        self.elements = None
        self.latticevec = None
        self.lattice_param = None

    def read_file(self):
        return

    def get_frac(self):
        retutn

    def get_xyz(self):
        return

    def get_lattice_param(self):
        return

    def supplement_sp(self, spacegroup):
        atoms = crystal(self.elements, [self.xyz], 
                        spacegroup=spacegroup,
                        cellpar=self.lattice_param)
        return atoms
        
