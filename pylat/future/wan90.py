from pylat.future.wan90_template import *
import os
import subprocess


class Wannier90:
    def __init__(self, qe_ctrl):
        self.qe_ctrl = qe_ctrl

    def write_wan90(self, win_min, win_max, nw=5):
        with open(self.qe_ctrl.log, "r") as file:
            for line in file:
                if "number of Kohn-Sham states" in line:
                    band_count = int(line.split("=")[1].strip())
                    print(f"Number of Kohn-Sham states (bands): {band_count}")
        self.qe_ctrl.parse_gauss()
        txt = ""
        txt += wan90_temp0.format(
            nb=band_count, nw=nw, dis_win_min=win_min, dis_win_max=win_max
        )
        txt += "begin projections \n"
        txt += ""
        for i in range(self.qe_ctrl.N_initial_guess):
            type_orb = self.qe_ctrl.gaussian_orb[i][0]
            coord_x = self.qe_ctrl.gauss_center[i][0]
            coord_y = self.qe_ctrl.gauss_center[i][1]
            coord_z = self.qe_ctrl.gauss_center[i][2]
            txt += f"f={coord_x}, {coord_y}, {coord_z}: {type_orb}\n"
        txt += "end projections \n"
        txt += wan90_temp1
        txt += "begin unit_cell_cart \n"
        txt += "angstrom \n"
        for l in self.qe_ctrl.lattice:
            txt += f"{l[0]:.10f}   {l[1]:.10f}   {l[2]:.10f} \n"
        txt += "end unit_cell_cart \n"
        txt += "begin atoms_frac \n"
        for at in self.qe_ctrl.geoms:
            txt += f"{at[0]}  {at[1][0]:.10f}  {at[1][1]:.10f}  {at[1][2]:.10f} \n"
        txt += "end atoms_frac \n"
        txt += wan90_temp2.format(
            nkx=self.qe_ctrl.kpoints[0],
            nky=self.qe_ctrl.kpoints[1],
            nkz=self.qe_ctrl.kpoints[2],
        )

        kxs, kys, kzs, w = self.qe_ctrl.get_kpoint()

        txt += "begin kpoints \n"
        for kx in kxs:
            for ky in kys:
                for kz in kzs:
                    txt += f"{kx} {ky} {kz}\n"

        txt += "end kpoints \n"
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
