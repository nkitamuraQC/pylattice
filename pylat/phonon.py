
import subprocess

def run_command(input_file, output_file, kernel="ph.x"):
    try:
        # Open input and output files
        with open(input_file, "r") as input_file_data, open(output_file, "w") as output_file_data:
            # Execute the command
            subprocess.run(
                [kernel],  # Command as a list
                stdin=input_file_data,  # Redirect standard input
                stdout=output_file_data,  # Redirect standard output
                check=True  # Raise an exception if the command fails
            )
        print(f"Phonon calculation completed. Output written to {self.prefix}.ph.out")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running phonon calculation: {e}")

class Phonon:
    """
    Class to handle phonon calculations.
    """

    def __init__(self, qe):
        self.qe = qe
        self.prefix = self.qe.prefix
        self.outdir = self.qe.outdir
        self.dfpt_str_list = [
            "prefix",
            "fildyn",
            "outdir",
        ]

        self.dfpt_section = [
            "tr2_ph",
	        "prefix",
	        "outdir",
	        "fildyn",
	        "ldisp",
            "epsil",
            "nq1",
            "nq2",
            "nq3",
        ]
        self.tr2_ph = "1.0d-14"
        self.fildyn = f'{self.prefix}.dyn'
        self.ldisp = ".true."
        self.flfrc = f'{self.prefix}.fc'
        self.epsil = ".true."
        self.nq1 = 2
        self.nq2 = 2
        self.nq3 = 2
        self.nk1 = 4
        self.nk2 = 4
        self.nk3 = 4

    def check_qpoint(self):
        if self.ldisp == ".true.":
            self.dfpt_section.append("nq1")
            self.dfpt_section.append("nq2")
            self.dfpt_section.append("nq3")
        return

    def make_input_ph(self, txt=""):
        # self.check_qpoint()
        txt += "\n"
        txt += f"&inputph\n"
        for item in self.dfpt_section:
            if hasattr(self, item):
                val = getattr(self, item)
                if item in self.dfpt_str_list:
                    txt += f"  {item} = '{val}',\n"
                else:
                    txt += f"  {item} = {val},\n"
        txt += f"/\n"
        return txt
    
    def make_input_q2r(self):
        txt = "&input \n"
        txt += f"  fildyn = '{self.fildyn}',\n"
        txt += f"  flfrc = '{self.flfrc}',\n"
        txt += "/"
        return txt
    
    def make_input_matdyn(self, txt=""):
        txt = "&input \n"
        txt += f"  flfrc = '{self.flfrc}',\n"
        txt += "  asr = 'crystal' \n"
        txt += "  dos = .true. \n"
        txt += f"  nk1 = {self.nk1},\n"
        txt += f"  nk2 = {self.nk2},\n"
        txt += f"  nk3 = {self.nk3},\n"
        txt += "/"
        return txt
    
    def write_input(self, inp_ph, inp_q2r, inp_matdyn):
        """
        Write the input file for phonon calculation.
        """
        with open(f"{self.prefix}.ph.in", "w") as f:
            f.write(inp_ph)

        with open(f"{self.prefix}.q2r.in", "w") as f:
            f.write(inp_q2r)

        with open(f"{self.prefix}.matdyn.in", "w") as f:
            f.write(inp_matdyn)
        return


    def run(self):
        """
        Run the phonon calculation.
        """
        inp_ph = self.make_input_ph()
        inp_q2r = self.make_input_q2r()
        inp_matdyn = self.make_input_matdyn()
        self.write_input(inp_ph, inp_q2r, inp_matdyn)
        run_command(self.prefix + ".ph.in", self.prefix + ".ph.out", kernel="ph.x")
        run_command(self.prefix + ".q2r.in", self.prefix + ".q2r.out", kernel="q2r.x")
        run_command(self.prefix + ".matdyn.in", self.prefix + ".matdyn.out", kernel="matdyn.x")
        return
