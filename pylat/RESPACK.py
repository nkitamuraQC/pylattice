#
import os
import numpy as np

def get_volume(qe_cls):
    lattice = qe_cls.lattice
    x = lattice[0]
    y = lattice[1]
    z = lattice[2]
    vec = np.cross(x, y)
    vol = np.dot(vec, z)
    return vol

def get_fracCoord(qe_cls, xyz):
    lattice = qe_cls.T
    norm = []
    v = get_volume(qe_cls)
    cos = np.cos
    sin = np.sin
    for i in range(3):
        norm_i = np.linalg.norm(lattice[i])
        norm.append(norm_i)
    norm = np.array(norm)
    alpha_ = np.dot(lattice[1], lattice[2])/(norm[1]*norm[2])
    beta_ = np.dot(lattice[0], lattice[2])/(norm[0]*norm[2])
    gamma_ = np.dot(lattice[0], lattice[1])/(norm[0]*norm[1])
    alpha = np.arccos(alpha_)
    beta = np.arccos(beta_)
    gamma = np.arccos(gamma_)
    tr_mat = np.zeros((3,3),dtype=float)
    a = norm[0]
    b = norm[1]
    c = norm[2]
    tr_mat[0, 0] = 1.0/a
    tr_mat[0, 1] = -cos(gamma)/(a*sin(gamma))
    tr_mat[0, 2] = b*c*cos(gamma)*(cos(alpha) - cos(beta)*cos(gamma))/(v*sin(gamma)) - b*c*cos(beta)*sin(gamma)/v
    tr_mat[1, 1] = 1/(b*sin(gamma))
    tr_mat[1, 2] = a*c*(cos(beta)*cos(gamma)-cos(alpha))/(v*sin(gamma))
    tr_mat[2, 2] = a*b*sin(gamma)/v
    frac = np.dot(tr_mat,xyz)
    return frac

class RESPACKController:
    def __init__(self, qe_cls):    
        self.debug = True
        self.workpath = "./"
        self.qe2respack_path = "/home/users/auv/pyscf/pyscf/pbc/downfolding"
        self.filename = "input.in"
        self.calctype = "calc_wannier"
        self.flg_cRPA = 1
        self.MPI_num_qcomm = 12
        self.N_wannier = 2

        self.Lower_energy_window = "-10.0d0"
        self.Upper_energy_window = "10.0d0"
        self.FLG_BMAT = 0
        self.N_initial_guess = 2 # 6
        self.gaussian_orb = None
        self.Bmat = None
        self.N_sym_points = 5
        self.bandpath = None
        self.RESPACK_path = "/home/nkitamura/respack/RESPACK-20190527-dist/"
        self.MPI = 1
        self.qe = qe_cls
        self.num_atom = None
        self.num_mole = None
        self.gauss = [[0, "px", 0.2]] ## [atom_index, orbital_type, orbital_exp]
        return

        
    def parse_gauss(self):
        gauss_orb = []
        gauss_center = []
        for g in self.gauss:
            gauss_orb.append([g[0], g[1]])
            gauss_center.append(self.qe.geoms[g[0]])
        self.gaussian_orb = gauss_orb
        self.gauss_center = gauss_center
        return

    def get_Bmat(self):
        raise NotImplementedError

    def get_bandpath(self):
        # TMTTF
        bandpath = [[0.000, 0.500, 0.000], [0.000, 0.000, 0.000],[0.500, 0.000, 0.000],\
                    [0.500, 0.500, 0.000],[0.000, 0.000, 0.000]]
        self.bandpath = bandpath
        assert(self.N_sym_points == len(bandpath))
        return

    def prepare(self):
        self.parse_gauss()
        if self.FLG_BMAT == 1:
            self.get_Bmat()
        self.get_bandpath()
        self.write_inp()
        return 

    def write_inp(self):
        txt = ""
        filename = self.workpath + "/" + self.filename
        txt += "&param_chiqw"+"\n"
        txt += "flg_cRPA={rpa}".format(rpa=self.flg_cRPA)+", \n"
        txt += "MPI_num_qcomm={mpi}".format(mpi=self.MPI_num_qcomm)+"/ \n"
        txt += "\n"
        txt += "&param_wannier"+"\n"
        txt += "N_wannier={nwan}".format(nwan = self.N_wannier)+", \n"
        txt += "Lower_energy_window={low}".format(low=self.Lower_energy_window)+", \n"
        txt += "Upper_energy_window={up}".format(up=self.Upper_energy_window)+", \n"
        txt += "FLG_BMAT={flag}".format(flag=self.FLG_BMAT)+", \n"
        txt += "N_initial_guess={Nguess}".format(Nguess=self.N_initial_guess)+"/ \n"
        for i in range(self.N_initial_guess):
            basis = self.N_initial_guess // self.N_wannier
            idx = i // basis
            type_orb = self.gaussian_orb[0][i]
            exp_orb = self.gaussian_orb[1][i]
            coord_x = self.gauss_center[idx, 0]
            coord_y = self.gauss_center[idx, 1]
            coord_z = self.gauss_center[idx, 2]
            txt += type_orb + " " + str(exp_orb) + " " + str(coord_x) \
                     + " " + str(coord_y) + " " + str(coord_z) + "\n"
    
        # Bmat
        for i in range(self.N_initial_guess):
            for j in range(self.N_wannier):
                txt += str(self.Bmat[i, j]) + " "
            txt += "\n"
    
        txt += "&param_interpolation"+"\n"
        txt += "N_sym_points={npoint}".format(npoint=self.N_sym_points)+"/ \n"
        bandpath = self.bandpath
        for i in range(self.N_sym_points):
            for j in range(3):
                txt += str(bandpath[i][j])+" "
            txt += "\n"
    
        txt += "&param_visualization"+"\n"
        txt += "/"+"\n"
    
        txt += "&param_calc_int"+"\n"
        txt += "/"+"\n"
        wf = open(filename, "w")
        wf.write(txt)
        wf.close()
        return 

    def execution(self):
        util = self.qe2respack_path
        work = self.workpath
        filename = work + "/" + self.filename
        dir_ = self.qe.prefix + ".save"
        if self.calctype == "calc_wannier":
            os.system("python {util}/qe2respack.py {work}/{dir_}".format(util=util,work=work,dir_=dir_))
            os.system("calc_wannier < {filename} > {filename}_wan.out".format(filename=filename))
        elif self.calctype == "calc_chiqw":
            if self.MPI > 1:
                os.system("OMP_NUM_THREADS=1 mpirun -np {MPI} calc_chiqw < {filename} > {filename}_chi.out".format(filename=filename,MPI=self.MPI))
            else:
                os.system("calc_chiqw < {filename} > {filename}_chi.out".format(filename=filename,MPI=self.MPI))
        elif self.calctype == "calc_w3d":
            os.system("calc_w3d < {filename} > {filename}_w.out".format(filename=filename))
        elif self.calctype == "calc_j3d":
            os.system("calc_j3d < {filename} > {filename}_j.out".format(filename=filename))
        return 