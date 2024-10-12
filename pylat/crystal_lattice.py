from ase.spacegroup import crystal
from pymatgen.core import Structure
from ase.io import write

class SpaceGroupSupplement:
    def __init__(self, ciffile, xyzfile=None, lattfile=None):
        self.ciffile = ciffile
        self.xyzfile = xyzfile
        self.lattfile = lattfile

        self.xyz = None
        self.elements = None
        self.latticevec = None
        self.lattice_param = None
        self.spacegroup = None

    def read_file(self):
        structure = Structure.from_file(self.ciffile)
        lattice = structure.lattice
        a, b, c = lattice.a, lattice.b, lattice.c  # 長さ
        alpha, beta, gamma = lattice.alpha, lattice.beta, lattice.gamma  # 角度
        
        # 空間群を取得
        space_group = structure.get_space_group_info()
        
        # 原子の座標を取得
        atomic_coordinates = structure.frac_coords  # 分数座標
        atomic_species = [site.species_string for site in structure] 

        self.xyz = atomic_coordinates
        self.elements = atomic_species
        self.spacegroup = space_group
        self.lattice_param = [a, b, c, alpha, beta, gamma]
        return

    def get_frac(self, cartesian_coords):
        structure = Structure.from_file(self.ciffile)
        lattice = structure.lattice
        fractional_coords = lattice.get_fractional_coords(cartesian_coords)
        return fractional_coords

    def get_xyz(self, fractional_coords):
        structure = Structure.from_file(self.ciffile)
        lattice = structure.lattice
        cartesian_coords = lattice.get_cartesian_coords(fractional_coords)
        return cartesian_coords

    def supplement(self):
        atoms = crystal(self.elements, [self.xyz], 
                        spacegroup=self.spacegroup,
                        cellpar=self.lattice_param)
        self.atoms = atoms
        return atoms
    
    def write_cif(self, cifname):
        write(cifname, self.atoms)
        return
        
