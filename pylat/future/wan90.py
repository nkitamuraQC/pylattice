from pylat.future.wan90_template import *
import os
import subprocess
from pylat.k_path import SeekPath

class Wannier90:
    def __init__(self, qe_ctrl):
        self.qe_ctrl = qe_ctrl
        self.prefix = qe_ctrl.prefix
        self.target_dicts = {"Fe": ["dxy", "dyz", "dxz", "dz2", "dx2-y2"]}

    def write_wan90(self, win_min, win_max, nw=5, use_seekpath=True):
        with open(self.qe_ctrl.log, "r") as file:
            for line in file:
                if "number of Kohn-Sham states" in line:
                    band_count = int(line.split("=")[1].strip())
                    print(f"Number of Kohn-Sham states (bands): {band_count}")
                    break
        # self.qe_ctrl.parse_gauss()
        txt = ""
        txt += wan90_temp0.format(
            nb=band_count, nw=nw, dis_win_min=win_min, dis_win_max=win_max
        )
        txt += "\n"
        txt += "begin Unit_Cell_Cart \n"
        txt += "Ang \n"
        for l in self.qe_ctrl.lattice:
            txt += f"{l[0]:.10f}   {l[1]:.10f}   {l[2]:.10f} \n"
        txt += "End Unit_Cell_Cart \n"
        txt += "\n"
        txt += "Begin Projections \n"
        sub_txt = ""
        for k, v in self.target_dicts.items():
            sub_txt += f"{k}:"
            for idx, v_ in enumerate(v):
                if idx < len(v) - 1:
                    sub_txt += f"{v_};"
                else:
                    sub_txt += f"{v_}\n"
        txt += sub_txt
        txt += "End Projections\n"
        txt += "\n"
        txt += "Begin  ATOMS_FRAC\n"
        for at in self.qe_ctrl.geoms:
            txt += f"{at[0]}  {at[1][0]:.10f}  {at[1][1]:.10f}  {at[1][2]:.10f} \n"
        txt += "End ATOMS_FRAC \n"
<<<<<<< Updated upstream
        for i in range(self.qe_ctrl.N_initial_guess):
            type_orb = self.qe_ctrl.gaussian_orb[i][0]
            coord_x = self.qe_ctrl.gauss_center[i][0]
            coord_y = self.qe_ctrl.gauss_center[i][1]
            coord_z = self.qe_ctrl.gauss_center[i][2]
            txt += f"f={coord_x}, {coord_y}, {coord_z}: {type_orb}\n"
        txt += "end projections \n"
=======
        txt += "\n"
        txt += "Begin Kpoint_path\n"
>>>>>>> Stashed changes
        if use_seekpath:
            k_txt = ""
            sp = SeekPath(self.qe_ctrl)
            sp.exec()
            path = sp.path
            kcoord = sp.coord
            for p in path:
                kcoord1 = kcoord[p[0]]
                kcoord2 = kcoord[p[1]]
                k_txt += f"{p[0]} {kcoord1[0]} {kcoord1[1]} {kcoord1[2]} {p[1]} {kcoord2[0]} {kcoord2[1]} {kcoord2[2]}\n"
            txt += k_txt
        else:
            txt += wan90_kpath
<<<<<<< Updated upstream
=======
        txt += "End Kpoint_path\n"
>>>>>>> Stashed changes
        txt += "\n"
        
        txt += wan90_temp2.format(
            nkx=self.qe_ctrl.kpoints[0],
            nky=self.qe_ctrl.kpoints[1],
            nkz=self.qe_ctrl.kpoints[2],
        )
        txt += "\n\n"

        kxs, kys, kzs, w = self.qe_ctrl.get_kpoint()

        txt += "begin kpoints \n"
        for kx in kxs:
            for ky in kys:
                for kz in kzs:
                    txt += f"{kx} {ky} {kz}\n"

        txt += "end kpoints \n"
        txt += "\n"
        win_in = f"{self.qe_ctrl.prefix}.win"
        wf = open(win_in, "w")
        wf.write(txt)
        wf.close()

        txt = wan90_pw2wan.format(
            prefix=self.qe_ctrl.prefix,
            outdir=self.qe_ctrl.outdir,
            win=f"{self.qe_ctrl.prefix}.win",
        )
        wf = open(f"{self.qe_ctrl.prefix}.pw2wan.in", "w")
        wf.write(txt)
        wf.close()
        return

    def parse_gauss(self):
        gauss_orb = []
        gauss_center = []
        for g in self.qe_ctrl.gauss:
            gauss_orb.append([g[1], g[2]])
            gauss_center.append(self.qe_ctrl.geoms[g[0]][1])
        self.gaussian_orb = gauss_orb
        self.gauss_center = gauss_center
        return 

    def do_wan90(self):
        # 1. wannier90.x -pp コマンドの実行
        subprocess.run(["wannier90.x", "-pp", self.prefix], check=True)

        # 2. pw2wannier90.x コマンドの実行
        pw2wan_in = f"{self.qe_ctrl.prefix}.pw2wan.in"
        pw2wan_out = f"{self.qe_ctrl.prefix}.pw2wan.out"
        with open(pw2wan_in, "r") as infile, open(pw2wan_out, "w") as outfile:
            subprocess.run(["pw2wannier90.x"], stdin=infile, stdout=outfile, check=True)

        # 3. 実行済みでコメントアウトされた部分の再実装（必要に応じて有効化）
        # subprocess.run(["wannier90.x", self.prefix], check=True)

        return
