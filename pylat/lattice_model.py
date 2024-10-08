from pylat.template import *
import pathlib, os
import numpy as np

class LatticeModel_gen:
    def __init__(self, lat_type, ints_dict, H, W, nelec):
        self.ints_dict = ints_dict
        self.size_dict = {"H": H, "W": W, "ncond": nelec}
        self.lat_type = lat_type
        self.vmcdry_path = ""
        self.tmp_dir = pathlib.Path("")
        self.lattice_dict = {"chain": chain, "triangular": triangular,
                             "square": square, "honeycomb": honeycomb,
                             "ladder": ladder, "kagome": kagome}

    def drive(self):
        self.ints_dict.update(self.size_dict)
        temp = self.lattice_dict[self.lat_type]
        temp = temp.format(**self.ints_dict)
        name = str(self.tmp_dir.joinpath("tmp.in"))
        wf = open(name, "w")
        wf.write(temp)
        wf.close()
        os.system(f"{self.vmcdry_path} {name}")
        return

    def get_ints(self):
        """
        read results

        Returns:
            tuple(np.array, np.array)
        """
        return
    
        
            
            
        
