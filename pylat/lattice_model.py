from pylat.template import *
import pathlib, os

class LatticeModel_gen:
    def __init__(self, lat_type, ints_dict):
        self.ints_dict = ints_dict
        self.lat_type = lat_type
        self.vmcdry_path = ""
        self.tmp_dir = pathlib.Path("")
        self.lattice_dict ={}

    def drive(self):
        temp = self.lattice_dict[self.lat_type]
        temp = temp.format(**self.ints_dict)
        wf = open(str(self.tmp_dir.joinpath("tmp.in")), "w")
        wf.write(temp)
        wf.close()
        os.system(f"{self.vmcdry_path} {
        
            
            
        
