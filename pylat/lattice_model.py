from pylat.template import *
import pathlib, os
import numpy as np

class LatticeModel_gen:
    def __init__(self, ints_dict, lat_type="chain", H=4, W=1, nelec=2):
        self.ints_dict = ints_dict
        self.H = H
        self.W = W
        self.size_dict = {"H": H, "W": W, "ncond": nelec}
        self.lat_type = lat_type
        self.vmcdry_path = "/Users/qclove00/mVMC-1.1.0/build/src/mVMC/vmcdry.out"
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

    def get_ints(self, start=5):
        """
        read results

        Returns:
            tuple(np.array, np.array)
        """
        norb = self.H * self.W
        int1e = np.zeros((norb, 2, norb, 2))
        int2e = np.zeros((norb, norb, norb, norb))

        rf = open("./output/trans.def", "r")
        lines = rf.readlines()[start:]
        for line in lines:
            i, si, j, sj, val = int(line[0]), int(line[1]), int(line[2]), int(line[3]), float(line[4])
            int1e[i, si, j, sj] = val

        rf = open("./output/coulombintra.def", "r")
        lines = rf.readlines()[start:]
        for line in lines:
            i, val = int(line[0]), float(line[1])
            int2e[i, i, i, i] = val

        rf = open("./output/coulombinter.def", "r")
        lines = rf.readlines()[start:]
        for line in lines:
            i, j, val = int(line[0]), int(line[1]), float(line[2])
            int2e[i, i, j, j] = val

        rf = open("./output/exchange.def", "r")
        lines = rf.readlines()[start:]
        for line in lines:
            i, j, val = int(line[0]), int(line[1]), float(line[2])
            int2e[i, j, j, i] = val
            int2e[i, j, i, j] = val
        return int1e, int2e
    
    def do(self):
        self.drive()
        int1e, int2e = self.get_ints()
        return int1e, int2e
    
        
            
            
        
