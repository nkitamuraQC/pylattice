#
import os
import numpy as np
from pylat.template_respack import *
import subprocess

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
    alpha_ = np.dot(lattice[1], lattice[2]) / (norm[1] * norm[2])
    beta_ = np.dot(lattice[0], lattice[2]) / (norm[0] * norm[2])
    gamma_ = np.dot(lattice[0], lattice[1]) / (norm[0] * norm[1])
    alpha = np.arccos(alpha_)
    beta = np.arccos(beta_)
    gamma = np.arccos(gamma_)
    tr_mat = np.zeros((3, 3), dtype=float)
    a = norm[0]
    b = norm[1]
    c = norm[2]
    tr_mat[0, 0] = 1.0 / a
    tr_mat[0, 1] = -cos(gamma) / (a * sin(gamma))
    tr_mat[0, 2] = (
        b * c * cos(gamma) * (cos(alpha) - cos(beta) * cos(gamma)) / (v * sin(gamma))
        - b * c * cos(beta) * sin(gamma) / v
    )
    tr_mat[1, 1] = 1 / (b * sin(gamma))
    tr_mat[1, 2] = a * c * (cos(beta) * cos(gamma) - cos(alpha)) / (v * sin(gamma))
    tr_mat[2, 2] = a * b * sin(gamma) / v
    frac = np.dot(tr_mat, xyz)
    return frac


class RESPACKController:
    def __init__(self, qe_cls):
        self.debug = True
        self.workpath = qe_cls.outdir
        self.prefix = qe_cls.prefix
        self.qe2respack_path = "/Users/qclove00/RESPACK/build/bin"
        self.filename = "input.in"
        self.calctype = "calc_wannier"
        self.flg_cRPA = 1
        self.MPI_num_qcomm = 12
        self.N_wannier = 5

        self.Lower_energy_window = "1.1049d01"
        self.Upper_energy_window = "1.8929d01"
        self.FLG_BMAT = 0
        self.N_initial_guess = 5
        self.gaussian_orb = None
        self.Bmat = None
        self.N_sym_points = 5
        self.bandpath = [
            [0.500, 0.500, 0.500],
            [0.000, 0.000, 0.000],
            [0.500, 0.000, 0.500],
            [0.500, 0.250, 0.750],
            [0.500, 0.500, 0.500],
        ]
        self.MPI = 1
        self.qe = qe_cls
        self.num_atom = None
        self.num_mole = None
        self.gauss = [
            [0, "dxy", 0.2],
            [0, "dyz", 0.2],
            [0, "dzx", 0.2],
            [0, "dz2", 0.2],
            [0, "dx2", 0.2],
        ]  ## [atom_index, orbital_type, orbital_exp]
        return

    def parse_gauss(self):
        gauss_orb = []
        gauss_center = []
        for g in self.gauss:
            gauss_orb.append([g[1], g[2]])
            gauss_center.append(self.qe.geoms[g[0]][1])
        self.gaussian_orb = gauss_orb
        self.gauss_center = gauss_center
        return

    def get_Bmat(self):
        raise NotImplementedError

    def prepare(self):
        self.parse_gauss()
        if self.FLG_BMAT == 1:
            self.get_Bmat()
        self.write_inp()
        return

    def write_inp(self):
        txt = template_res1.format(
            N_wannier=self.N_wannier,
            lower_e=self.Lower_energy_window,
            upper_e=self.Upper_energy_window,
            n_init_guess=self.N_initial_guess,
        )
        for i in range(self.N_initial_guess):
            type_orb = self.gaussian_orb[i][0]
            exp_orb = self.gaussian_orb[i][1]
            coord_x = self.gauss_center[i][0]
            coord_y = self.gauss_center[i][1]
            coord_z = self.gauss_center[i][2]
            txt += f"{type_orb} {exp_orb} {coord_x} {coord_y} {coord_z} \n"

        txt += template_res2
        wf = open(self.filename, "w")
        wf.write(txt)
        wf.close()
        return

    def execution(self):
        util = self.qe2respack_path
        work = self.workpath
        filename = self.filename

        if self.calctype == "calc_wannier":
            # 実行コマンド1: python qe2respack.py
            subprocess.run(
                ["python", f"{util}/qe2respack.py", f"{work}/{self.prefix}.save"],
                check=True
            )
            # 実行コマンド2: calc_wannier
            with open(filename, "r") as infile, open(f"{filename}_wan.out", "w") as outfile:
                subprocess.run(
                    [f"{util}/calc_wannier"], 
                    stdin=infile,
                    stdout=outfile,
                    check=True
                )

        elif self.calctype == "calc_chiqw":
            if self.MPI > 1:
                # 実行コマンド: MPIありの場合
                subprocess.run(
                    [
                        "mpirun",
                        "-np",
                        str(self.MPI),
                        "calc_chiqw"
                    ],
                    stdin=open(filename, "r"),
                    stdout=open(f"{filename}_chi.out", "w"),
                    env={"OMP_NUM_THREADS": "1"},
                    check=True
                )
            else:
                # 実行コマンド: MPIなしの場合
                with open(filename, "r") as infile, open(f"{filename}_chi.out", "w") as outfile:
                    subprocess.run(
                        ["calc_chiqw"],
                        stdin=infile,
                        stdout=outfile,
                        check=True
                    ) 

        elif self.calctype == "calc_w3d":
            # 実行コマンド: calc_w3d
            with open(filename, "r") as infile, open(f"{filename}_w.out", "w") as outfile:
                subprocess.run(
                    ["calc_w3d"],
                    stdin=infile,
                    stdout=outfile,
                    check=True
                )

        elif self.calctype == "calc_j3d":
            # 実行コマンド: calc_j3d
            with open(filename, "r") as infile, open(f"{filename}_j.out", "w") as outfile:
                subprocess.run(
                    ["calc_j3d"],
                    stdin=infile,
                    stdout=outfile,
                    check=True
                )

        return

