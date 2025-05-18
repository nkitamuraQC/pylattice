from ase.spacegroup import crystal
from pymatgen.core import Structure
from ase.io import write

from ase.io import read
from ase.spacegroup import get_spacegroup

from ase.io import read
import spglib

def get_spacegroup_from_cif_spglib(cif_file_path, symprec=1e-5):
    """
    CIFファイルからspglibを使って空間群を取得する関数

    Parameters:
        cif_file_path (str): CIFファイルへのパス
        symprec (float): 対称性の検出に使う精度（デフォルト: 1e-5）

    Returns:
        str: 空間群の記号（例: 'Fm-3m (225)'）
    """
    atoms = read(cif_file_path)

    # spglibのセル形式に変換
    lattice = atoms.get_cell()
    positions = atoms.get_scaled_positions()
    numbers = atoms.get_atomic_numbers()
    cell = (lattice, positions, numbers)

    # 空間群情報の取得
    spacegroup = spglib.get_spacegroup(cell, symprec=symprec)

    return spacegroup


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


    def get_sg(self):
        sg = get_spacegroup_from_cif_spglib(self.ciffile)
        print(f"Space Group: {sg}")
        return
        
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
        atoms = crystal(
            self.elements,
            self.xyz,
            spacegroup=self.spacegroup[0],
            cellpar=self.lattice_param,
        )
        self.atoms = atoms
        return atoms

    def write_cif(self, cifname):
        write(cifname, self.atoms)
        return
