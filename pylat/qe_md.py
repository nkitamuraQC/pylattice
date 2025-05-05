import copy, os
import numpy as np
import subprocess

class QEMD:
    def __init__(
        self,
        qe_ctrl,
        dt=20,
        micro_step=50,
        prefix="calculation",
        dL=np.zeros((3, 3)),
        temp=300,
        macro_step=10,
        nkpoints = [4, 4, 4],
        press = None, 
        degauss = None,
    ):
        self.qe_ctrl = qe_ctrl
        self.dL = dL
        self.prefix = prefix
        self.macro_step = macro_step
        self.micro_step = micro_step
        self.qe_ctrl.tempw = temp
        self.dt = dt
        self.nkpoints = nkpoints
        self.default_lattice = copy.copy(self.qe_ctrl.lattice)
        self.press = press
        self.degauss = degauss

    def run(self, n_para=None):
        start = 0
        for istep_ in range(self.macro_step):
            if not os.path.exists(f"./{self.prefix}_try_{istep_}.out"):
                break
        start = istep_
        for istep in range(start, self.macro_step):
            self.main_step(istep, n_para)
        return

    def main_step(self, istep, n_para):
        if self.press is None:
            current_lattice = self.default_lattice + self.dL * istep
            self.qe_ctrl.press = None
            self.qe_ctrl.lattice = current_lattice
            self.qe_ctrl.kpoints = self.nkpoints
            self.qe_ctrl.calculation = "md"
        else:
            self.qe_ctrl.nosym = True
            self.qe_ctrl.calculation = "vc-md"
            self.qe_ctrl.press = self.press
        self.qe_ctrl.dt = self.dt
        self.qe_ctrl.nosym = True
        self.qe_ctrl.nstep = self.micro_step
        self.qe_ctrl.prefix = f"{self.prefix}_try_{istep}"
        self.qe_ctrl.outdir = f"./work/try_{istep}"
        if isinstance(self.degauss, float):
            self.qe_ctrl.occupations = "smearing"
            self.qe_ctrl.degauss = self.degauss
        # os.system(f"mkdir -p {self.qe_ctrl.outdir}")
        subprocess.run(["mkdir", "-p", self.qe_ctrl.outdir], check=True)
        if n_para is None:
            self.qe_ctrl.exec()
        else:
            self.qe_ctrl.para_exec(n_para=n_para)
        return
