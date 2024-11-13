import os
from pymatgen.core import Structure
import periodictable
from pylat.wan90_template import *
import numpy as np


class QEController:
    def __init__(self, cif, pseudo_dict):
        self.filename = "calculation.in"
        self.prefix = "calculation"
        self.calculation = "scf"
        self.outdir = "./work"
        self.pseudo_dir = "./pseudo"
        self.charge = 0
        self.ibrav = 0
        self.ecutwfc = 30.0
        self.ecutrho = 150.0
        self.occupations = 'tetrahedra'
        self.mixing_beta = 0.7
        self.conv_thr = "1.0d-8"
        self.tstress = True
        self.tprnfor = True
        self.wf_collect = True
        self.kpoints = [4, 4, 4]
        self.section = {"&control":["prefix","calculation","outdir","pseudo_dir","tstress","tprnfor","wf_collect"],\
                        "&system":["ibrav","nat","ntyp","ecutwfc","ecutrho","occupations","nbnd"],\
                        "&electrons":["mixing_beta","conv_thr"]}
        self.pseudo_dict = pseudo_dict
        self.nbnd = None

        self.geoms = None # [["H",[0,0,0]],..]
        self.atoms = None # [["H",mass,psudo]]
        self.lattice = None # [[x,y,z],[x,y,z],[x,y,z]]

        self.read_file(cif)

        self.nat = len(self.geoms)
        self.ntyp = len(self.atoms)
        self.gauss = [[0, "dxy", 0.2], [0, "dyz", 0.2], [0, "dzx", 0.2],  [0, "dz2", 0.2],  [0, "dx2", 0.2]]
        self.N_initial_guess = 4

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

            txt += "/"+"\n"
        txt += "ATOMIC_SPECIES\n"
        for at in self.atoms:
            txt += f"{at[0]}  {at[1]}  {at[2]} \n"
        txt += "ATOMIC_POSITIONS\n"
        for at in self.geoms:
            txt += f"{at[0]}  {at[1][0]:.10f}  {at[1][1]:.10f}  {at[1][2]:.10f} \n"
        txt += "CELL_PARAMETERS angstrom \n"
        for l in self.lattice:
            txt += f"{l[0]:.10f}   {l[1]:.10f}   {l[2]:.10f} \n"
        
        txt += "K_POINTS {" + "automatic" + "}\n"
        txt += f"{self.kpoints[0]} {self.kpoints[1]} {self.kpoints[2]} 0 0 0"
        return txt

    def write_input(self, inp):
        wf = open(self.filename, "w")
        wf.write(inp)
        wf.close()
        return

    def exec(self):
        txt = self.make_input()
        self.write_input(txt)
        os.system("pw.x {}".format(self.filename))
        return

    def exec_dos(self):
        txt = f"""&dos
    outdir = '{self.outdir}',
    prefix= '{self.prefix}',
    fildos= '{self.prefix}.dos',
/
        """
        dos = open("dos.in", "w")
        dos.write(txt)
        dos.close()
        os.system("dos.x < dos.in > dos.out")
        return
    
    def exec_pdos(self):
        txt = f"""&projwfc
    outdir = '{self.outdir}',
    prefix= '{self.prefix}',
/
        """
        dos = open("pdos.in", "w")
        dos.write(txt)
        dos.close()
        os.system("projwfc.x < pdos.in > pdos.out")
        return

    def write_wan90(self, nb, nw, a, b):
        self.parse_gauss()
        txt = ""
        txt += wan90_temp0.format(nb=nb, nw=nw, a=a, b=b)
        txt += "begin projections \n"
        for i in range(self.N_initial_guess):
            type_orb = self.gaussian_orb[i][0]
            exp_orb = self.gaussian_orb[i][1]
            coord_x = self.gauss_center[i][0]
            coord_y = self.gauss_center[i][1]
            coord_z = self.gauss_center[i][2]
            txt += f"f={coord_x}, {coord_y}, {coord_z}: {type_orb}\n"
        txt += "end projections \n"
        txt += wan90_temp1 
        txt += "begin unit_cell_cart \n"
        txt += "bohr \n"
        for l in self.lattice:
            txt += f"{l[0]:.10f}   {l[1]:.10f}   {l[2]:.10f} \n"
        txt += "end unit_cell_cart \n"
        txt += "begin atoms_frac \n"
        for at in self.geoms:
            txt += f"{at[0]}  {at[1][0]:.10f}  {at[1][1]:.10f}  {at[1][2]:.10f} \n"
        txt += "end atoms_frac \n"
        txt += wan90_temp2.format(nkx=self.kpoints[0], nky=self.kpoints[1], nkz=self.kpoints[2])

        kxs = np.linspace(0, 1, self.kpoints[0])
        kys = np.linspace(0, 1, self.kpoints[1])
        kzs = np.linspace(0, 1, self.kpoints[2])

        txt += "begin kpoints \n"
        for kx in kxs:
            for ky in kys:
                for kz in kzs:
                    txt += f"{kx} {ky} {kz}\n"
            
        txt += "end kpoints \n"
        wf = open(f"{self.prefix}.win.in", "w")
        wf.write(txt)
        wf.close()

        txt = wan90_pw2wan.format(prefix=self.prefix, outdir=self.outdir)
        wf = open(f"{self.prefix}.pw2wan.in", "w")
        wf.write(txt)
        wf.close()
        return

    def parse_gauss(self):
        gauss_orb = []
        gauss_center = []
        for g in self.gauss:
            gauss_orb.append([g[1], g[2]])
            gauss_center.append(self.geoms[g[0]][1])
        self.gaussian_orb = gauss_orb
        self.gauss_center = gauss_center
        return

    def do_wan90(self):
        os.system(f"wannier90.x -pp {self.prefix}")
        os.system(f"pw2wannier90.x < {self.prefix}.pw2wan.in > {self.prefix}.pw2wan.out")
        os.system(f"wannier90.x {self.prefix}")
        return

    def read_hessians(self):
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
