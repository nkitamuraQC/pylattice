import os
from pymatgen.core import Structure
import periodictable
import numpy as np
import pathlib


def get_section(occ):
    if "tetrahedra" in occ:
        section = {
            "&control": [
                "prefix",
                "calculation",
                "outdir",
                "pseudo_dir",
                "tstress",
                "tprnfor",
                "wf_collect",
            ],
            "&system": [
                "ibrav",
                "nat",
                "ntyp",
                "ecutwfc",
                "ecutrho",
                "occupations",
                "nosym",
                "noinv",
            ],
            "&electrons": ["mixing_beta", "conv_thr", "electron_maxstep"],
        }

    elif occ == "smearing":
        section = {
            "&control": [
                "prefix",
                "calculation",
                "outdir",
                "pseudo_dir",
                "tstress",
                "tprnfor",
                "wf_collect",
            ],
            "&system": [
                "ibrav",
                "nat",
                "ntyp",
                "ecutwfc",
                "ecutrho",
                "occupations",
                "nosym",
                "noinv",
                "smearing",
                "degauss",
            ],
            "&electrons": ["mixing_beta", "conv_thr", "electron_maxstep"],
        }
    return section


class QEController:
    def __init__(self, cif, pseudo_dict):
        self.filename = "calculation.in"
        self.log = "calculation.out"
        self.prefix = "calculation"
        self.calculation = "scf"
        self.outdir = "./work"
        self.pseudo_dir = "./pseudo"
        self.charge = 0
        self.ibrav = 0
        self.ecutwfc = 30.0
        self.ecutrho = 150.0
        self.occupations = "tetrahedra"
        self.mixing_beta = 0.7
        self.conv_thr = "1.0d-8"
        self.tstress = True
        self.tprnfor = True
        self.wf_collect = True
        self.nosym = False
        self.noinv = False
        self.kpoints = [4, 4, 4]
        self.crystal_kpoint = False
        self.smearing = "fd"
        self.degauss = 0.07
        self.offset = [0, 0, 0]
        self.electron_maxstep = 200

        self.pseudo_dict = pseudo_dict
        self.nbnd = None

        self.geoms = None  # [["H",[0,0,0]],..]
        self.atoms = None  # [["H",mass,psudo]]
        self.lattice = None  # [[x,y,z],[x,y,z],[x,y,z]]

        self.read_file(cif)

        self.nat = len(self.geoms)
        self.ntyp = len(self.atoms)
        self.gauss = [
            [0, "dxy", 0.2],
            [0, "dyz", 0.2],
            [0, "dxz", 0.2],
            [0, "dz2", 0.2],
            [0, "dx2-y2", 0.2],
        ]
        self.N_initial_guess = 5

    def get_kpoint(self):
        dx = 1 / self.kpoints[0]
        dy = 1 / self.kpoints[1]
        dz = 1 / self.kpoints[2]

        kxs = [dx * i for i in range(self.kpoints[0])]
        kys = [dy * i for i in range(self.kpoints[1])]
        kzs = [dz * i for i in range(self.kpoints[2])]
        weight = dx * dy * dz
        return kxs, kys, kzs, weight

    def read_file(self, ciffile):
        structure = Structure.from_file(ciffile)
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
        self.lattice = lattice.matrix
        self.xyz_cart = lattice.get_cartesian_coords(atomic_coordinates)

        self.geoms = []
        for i in range(len(self.xyz)):
            self.geoms.append([self.elements[i], self.xyz[i]])

        elems = list(set(self.elements))
        self.atoms = []
        save = []
        for elem in elems:
            if elem in save:
                continue
            element = periodictable.elements.symbol(elem)
            p = self.pseudo_dict[elem]
            weight = element.mass
            self.atoms.append([elem, weight, p])
            save.append(elem)

        return

    def make_input(self, txt=""):
        self.section = get_section(self.occupations)
        for k in self.section:
            print(k)
            txt += k + "\n"
            for a in self.section[k]:
                if isinstance(self[a], bool) and self[a] == True:
                    val = ".true."
                elif isinstance(self[a], bool) and self[a] == False:
                    val = ".false."
                else:
                    val = self[a]
                if isinstance(val, str):
                    if a != "conv_thr" and val != ".true." and val != ".false.":
                        txt += f"{a} = '{val}', \n"
                    else:
                        txt += f"{a} = {val}, \n"
                else:
                    txt += f"{a} = {val}, \n"
                if a == "nbnd" and self.calculation == "nscf":
                    txt += f"{a} = {val}, \n"

            txt += "/" + "\n"
        txt += "ATOMIC_SPECIES\n"
        for at in self.atoms:
            txt += f"{at[0]}  {at[1]}  {at[2]} \n"
        txt += "ATOMIC_POSITIONS\n"
        for at in self.geoms:
            txt += f"{at[0]}  {at[1][0]:.10f}  {at[1][1]:.10f}  {at[1][2]:.10f} \n"
        txt += "CELL_PARAMETERS angstrom \n"
        for l in self.lattice:
            txt += f"{l[0]:.10f}   {l[1]:.10f}   {l[2]:.10f} \n"

        if not self.crystal_kpoint:
            txt += "K_POINTS {" + "automatic" + "}\n"
            txt += f"{self.kpoints[0]} {self.kpoints[1]} {self.kpoints[2]} {self.offset[0]} {self.offset[1]} {self.offset[2]}"
        else:
            txt += "K_POINTS crystal \n"
            nk = self.kpoints[0] * self.kpoints[1] * self.kpoints[2]
            txt += f"{nk} \n"
            kxs, kys, kzs, w = self.get_kpoint()
            for kx in kxs:
                for ky in kys:
                    for kz in kzs:
                        txt += f"{kx} {ky} {kz} {w}\n"

        return txt

    def write_input(self, inp):
        wf = open(self.filename, "w")
        wf.write(inp)
        wf.close()
        return

    def exec(self):
        txt = self.make_input()
        self.write_input(txt)
        os.system("pw.x < {} > {}".format(self.filename, self.log))
        return

    def exec_dos(self):
        txt = f"""&dos
    outdir = '{str(self.outdir)}',
    prefix= '{self.prefix}',
    fildos= '{self.prefix}.dos',
/
        """
        # dos = open(str(self.input_folder.joinpath("dos.in")), "w")
        # dos_out = str(self.input_folder.joinpath("dos.out"))
        dos = open("dos.in", "w")
        dos_out = "dos.out"
        dos.write(txt)
        dos.close()
        os.system(f"dos.x < {dos} > {dos_out}")
        return

    def exec_pdos(self):
        txt = f"""&projwfc
    outdir = '{str(self.outdir)}',
    prefix= '{self.prefix}',
/
        """
        # dos = open(str(self.input_folder.joinpath("pdos.in")), "w")
        # dos_out = str(self.input_folder.joinpath("pdos.out"))
        dos = open("pdos.in", "w")
        dos_out = "pdos.out"
        dos.write(txt)
        dos.close()
        os.system(f"projwfc.x < {dos} > {dos_out}")
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
