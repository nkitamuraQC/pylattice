import os
from pymatgen.core import Structure
import periodictable
import numpy as np
import pathlib
from pylat.get_section_cryspy import get_section_for_cryspy

def get_section(myclass, occ):
    if myclass.calculation == "vc-md":
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
                    "dt",
                    "nstep",
                    "restart_mode",
                ],
                "&system": [
                    "ibrav",
                    "nat",
                    "ntyp",
                    "nosym",
                    "ecutwfc",
                    "ecutrho",
                    "occupations",
                ],
                "&electrons": ["conv_thr"],
                "&ions": ["ion_temperature", "tempw"],
                "&cell": ["press"],
            }

        if "smearing" in occ:
            section = {
                "&control": [
                    "prefix",
                    "calculation",
                    "outdir",
                    "pseudo_dir",
                    "tstress",
                    "tprnfor",
                    "wf_collect",
                    "dt",
                    "nstep",
                    "restart_mode",
                ],
                "&system": [
                    "ibrav",
                    "nat",
                    "ntyp",
                    "nosym",
                    "ecutwfc",
                    "ecutrho",
                    "occupations",
                    "degauss",
                ],
                "&electrons": ["conv_thr"],
                "&ions": ["ion_temperature", "tempw"],
                "&cell": ["press"],
            }

    elif myclass.calculation == "md":
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
                    "dt",
                    "nstep",
                    "restart_mode",
                ],
                "&system": [
                    "ibrav",
                    "nat",
                    "ntyp",
                    "nosym",
                    "ecutwfc",
                    "ecutrho",
                    "occupations",
                ],
                "&electrons": ["conv_thr"],
                "&ions": ["ion_temperature", "tempw"],
                "&cell": [],
            }

        elif "smearing" in occ:
            section = {
                "&control": [
                    "prefix",
                    "calculation",
                    "outdir",
                    "pseudo_dir",
                    "tstress",
                    "tprnfor",
                    "wf_collect",
                    "dt",
                    "nstep",
                    "restart_mode",
                ],
                "&system": [
                    "ibrav",
                    "nat",
                    "ntyp",
                    "nosym",
                    "ecutwfc",
                    "ecutrho",
                    "occupations",
                    "degauss",
                ],
                "&electrons": ["conv_thr"],
                "&ions": ["ion_temperature", "tempw"],
                "&cell": [],
            }
    elif myclass.calculation == "vc-relax":
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
                    "dt",
                    "nstep",
                    "restart_mode",
                ],
                "&system": [
                    "ibrav",
                    "nat",
                    "ntyp",
                    "nosym",
                    "ecutwfc",
                    "ecutrho",
                    "occupations",
                ],
                "&electrons": ["conv_thr"],
                "&ions": ["ion_temperature", "tempw"],
                "&cell": ["press"],
            }

        if "smearing" in occ:
            section = {
                "&control": [
                    "prefix",
                    "calculation",
                    "outdir",
                    "pseudo_dir",
                    "tstress",
                    "tprnfor",
                    "wf_collect",
                    "dt",
                    "nstep",
                    "restart_mode",
                ],
                "&system": [
                    "ibrav",
                    "nat",
                    "ntyp",
                    "nosym",
                    "ecutwfc",
                    "ecutrho",
                    "occupations",
                    "degauss",
                ],
                "&electrons": ["conv_thr"],
                "&ions": ["ion_temperature", "tempw"],
                "&cell": ["press"],
            }

    if myclass.calculation == "relax":
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
                    "dt",
                    "nstep",
                    "restart_mode",
                ],
                "&system": [
                    "ibrav",
                    "nat",
                    "ntyp",
                    "nosym",
                    "ecutwfc",
                    "ecutrho",
                    "occupations",
                ],
                "&electrons": ["conv_thr"],
                "&ions": ["ion_temperature", "tempw"],
                "&cell": ["press"],
            }

        if "smearing" in occ:
            section = {
                "&control": [
                    "prefix",
                    "calculation",
                    "outdir",
                    "pseudo_dir",
                    "tstress",
                    "tprnfor",
                    "wf_collect",
                    "dt",
                    "nstep",
                    "restart_mode",
                ],
                "&system": [
                    "ibrav",
                    "nat",
                    "ntyp",
                    "nosym",
                    "ecutwfc",
                    "ecutrho",
                    "occupations",
                    "degauss",
                ],
                "&electrons": ["conv_thr"],
                "&ions": ["ion_temperature", "tempw"],
                "&cell": ["press"],
            }

    else:
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
                    "restart_mode",
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
                    "restart_mode",
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
                    "degauss",
                ],
                "&electrons": ["mixing_beta", "conv_thr", "electron_maxstep"],
            }
    return section


class QEController:
    def __init__(self, cif, pseudo_dict, supercell=None, coord_type="crystal"):
        self.filename = None
        self.log = None
        self.prefix = "calculation"
        self.calculation = "scf"
        self.outdir = "./work"
        self.restart_mode = "from_scratch"
        self.pseudo_dir = "./pseudo"
        self.charge = 0
        self.ibrav = 0
        self.ecutwfc = 60.0
        self.ecutrho = 240.0
        self.etot_conv_thr = "1.0d-8"
        self.occupations = "tetrahedra"
        self.mixing_beta = 0.7
        self.press = "0.d0"
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
        self.dt = 20
        self.nstep = 50
        self.conv_thr = "1.0d-8"
        self.ion_temperature = "initial"
        self.tempw = 300.0
        self.coord_type = coord_type
        self.do_cryspy = False

        self.pseudo_dict = pseudo_dict
        self.nbnd = None

        self.geoms = None  # [["H",[0,0,0]],..]
        self.atoms = None  # [["H",mass,psudo]]
        self.lattice = None  # [[x,y,z],[x,y,z],[x,y,z]]

        self.read_file(cif, supercell=supercell)

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

    def read_file(self, ciffile, supercell=None):
        structure = Structure.from_file(ciffile)
        if supercell is not None:
            scaling_matrix = np.diag(np.array(supercell))
            supercell_struct = (
                structure.copy()
            )  # To keep the original structure unchanged
            supercell_struct.make_supercell(scaling_matrix)
            structure = supercell_struct
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
        
        if self.coord_type == "crystal":
            self.coord = self.xyz
        elif self.coord_type == "angstrom":
            self.coord = self.xyz_cart
        else:
            raise NotImplementedError
        self.geoms = []
        for i in range(len(self.xyz)):
            if ":" in self.elements[i]:
                self.elements[i] = self.elements[i].split(":")[0]
            self.geoms.append([self.elements[i], self.coord[i]])

        elems = list(set(self.elements))
        self.atoms = []
        save = []
        for elem in elems:
            if ":" in elem:
                elem = elem.split(":")[0]
            if elem in save:
                continue
            element = periodictable.elements.symbol(elem)
            p = self.pseudo_dict[elem]
            weight = element.mass
            self.atoms.append([elem, weight, p])
            save.append(elem)

        return
    

    def make_input(self, txt=""):
        if self.do_cryspy:
            self.section = get_section_for_cryspy(self, self.occupations)
        else:
            self.section = get_section(self, self.occupations)
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
        txt += f"ATOMIC_POSITIONS {self.coord_type} \n"
        for at in self.geoms:
            txt += f"{at[0]}  {at[1][0]:.10f}  {at[1][1]:.10f}  {at[1][2]:.10f} \n"
        txt += "CELL_PARAMETERS angstrom \n"
        for l in self.lattice:
            txt += f"{l[0]:.10f}   {l[1]:.10f}   {l[2]:.10f} \n"

        if not self.crystal_kpoint:
            txt += "K_POINTS {" + "automatic" + "}\n"
            txt += f"{self.kpoints[0]} {self.kpoints[1]} {self.kpoints[2]} {self.offset[0]} {self.offset[1]} {self.offset[2]}"
        else:
            txt += "K_POINTS \n"
            nk = self.kpoints[0] * self.kpoints[1] * self.kpoints[2]
            txt += f"{nk} \n"
            kxs, kys, kzs, w = self.get_kpoint()
            for kx in kxs:
                for ky in kys:
                    for kz in kzs:
                        txt += f"{kx} {ky} {kz} {w}\n"

        return txt
    
    def make_input_for_cryspy(self, txt=""):
        self.section = get_section(self, self.occupations)
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

        return txt

    def write_input(self, inp):
        self.filename = self.prefix+".in"
        self.log = self.prefix+".out"
        wf = open(self.prefix+".in", "w")
        wf.write(inp)
        wf.close()
        return


    def exec(self):
        txt = self.make_input()
        self.write_input(txt)
        os.system("pw.x < {} > {}".format(self.filename, self.log))
        return

    def para_exec(self, n_para=4):
        txt = self.make_input()
        self.write_input(txt)
        # mpirun -np 8 pw.x -in input.in > output.out
        os.system("mpirun -np {} pw.x < {} > {}".format(n_para, self.filename, self.log))
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
        dos_in = "dos.in"
        dos_out = "dos.out"
        dos.write(txt)
        dos.close()
        os.system(f"dos.x < {dos_in} > {dos_out}")
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
        dos_in = "pdos.in"
        dos_out = "pdos.out"
        dos.write(txt)
        dos.close()
        os.system(f"projwfc.x < {dos_in} > {dos_out}")
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
