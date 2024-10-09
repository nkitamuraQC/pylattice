import os


class QEController:
    def __init__(self, geoms, atoms, lattice):
        self.filename = "calculation.in"
        self.prefix = "calculation"
        self.calculation = "scf"
        self.outdir = "./"
        self.psudo_dir = "./pseudo"
        self.charge = 0
        self.ibrav = 0
        self.nat = len(geoms)
        self.ntyp = len(atoms)
        self.ecutwfc = 30.0,
        self.ecutrho = 150.0,
        self.occupations = 'tetrahedra'
        self.mixing_beta = 0.7
        self.conv_thr = "1.0d-8"
        self.section = {"&control":["prefix","calculation","outdir","psudo_dir"],\
                "&system":["ibrav","nat","ntyp","ecutwfc","ecutrho","occupations"],\
                "&electron":["mixing_beta","conv_thr"]}

        self.geoms = geoms # [["H",[0,0,0]],..]
        self.atoms = atoms # [["H",mass,psudo]]
        self.lattice = lattice # [[x,y,z],[x,y,z],[x,y,z]]



    def make_input(self, txt=""):
        for k in self.section:
            txt += k + "\n"
            for a in self.section[k]:
                txt += a + " = " + self[k] + "\n"
            txt += "/"+"\n"
        txt += "ATOMIC_SPECIES\n"
        for at in self.atoms:
            txt += at[0] + " " + at[1] + " " + at[2] + "\n"
        txt += "ATOMIC_POSITIONS\n"
        for at in self.geoms:
            txt += at[0] + " " + at[1][0] + " " + at[1][1] + " " + at[1][2] + "\n"
        txt += "CELL_PARAMETERS angstrom \n"
        for l in self.lattice:
            txt += l[0] + " " + l[1] + " " + l[2] + "\n"
        return txt

    def write_input(self, inp)
        wf = open(self.filename, "w")
        wf.write(inp)
        wf.close()
        return

    def exec(self):
        txt = self.make_input()
        self.write_input(txt)
        os.system("pw.x {}".format(self.filename))
        return

    def read_geom(self):
        return

    def read_mo_coeff(self):
        return

    def read_integrals(self):
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

if __name__ == "__main__":
    geom = ["O",[0.0, 0.0, 0.0],
            "H",[0.0, 0.0, 1.0],
            "H",[0.0, 1.0, 0.0]]
    atoms = [["O","AAA.UPF",32],["H","AAA.UPF",2]]
    lattice = [[10,0,0],[0,10,0],[0,0,10]]
    myclass = qe_controller(geom, atoms, lattice)
    inp = myclass.make_input()
    myclass.write_input(inp)
    myclass.exec()