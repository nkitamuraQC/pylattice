#####
import numpy
import matplotlib as mpl
mpl.use('tkagg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import os, copy



class Vibration:
    def __init__(self):
        self.prefix = "TMTTF2AsF6"
        self.inputfile = "TMTTF_qe_cp_{i}.in"
        self.outputfile = "TMTTF_qe_cp_0.out"
        self.workdir = os.getcwd()
        self.basis_vec = numpy.zeros((3,3),dtype=float)
        self.lattice = None
        self.atom_types = []
        self.atom_types_ = []
        self.coordinate = []
        self.mode = []
        self.ntyp = 5
        self.natoms = 59
        self.MD = False
        self.start_elem = "C"
        self.q_list = []
        self.wavenumber = []
        self.frequency = []
        self.nstep = 8000
        self.scale_length = None
        self.large_displace = None
        self.bohr = 0.529177
        self.atom_symbols = []
        self.natoms_element = {"C":20,"H":24,"S":8,"As":1,"F":6}
        self.atom_number = {"C":6,"S":16,"H":1,"As":33,"F":9,"Sb":51,"P":15}
        self.MD_results = {"coord":[],"nfi":[],"time(ps)":[],"ekinc":[],\
                           "T_cell(K)":[],"Tion(K)":[],"etot":[],"enthal":[]}
        return


    
    def read_freq(self):
        rf = open(self.inputfile,"r")
        lines = rf.readlines()
        len_lines = len(lines)
        i =  0
        count_freq = 0
        while i < len_lines:
            if "Diagonalizing the dynamical matrix" in lines[i]:
                j = i + 2
                q_arr = numpy.array(list(map(float,lines[j].split()[3:6])))
                self.q_list.append(q_arr)
            if "freq" in lines[i]:
                self.frequency.append(float(lines[i].split()[4]))
                self.wavenumber.append(float(lines[i].split()[7]))
                tmp = []
                for j in range(1,self.natoms+1):
                    coord_x = float(lines[i+j].split()[1])
                    coord_y = float(lines[i+j].split()[3])
                    coord_z = float(lines[i+j].split()[5])
                    displace = numpy.array([coord_x,coord_y,coord_z])
                    #print("displace = ")
                    #print(displace)
                    tmp.append(displace)
                self.mode.append(tmp)
            i += 1
        return

    def read_coord(self):
        rf = open(self.inputfile,"r")
        lines = rf.readlines()
        len_lines = len(lines)
        i =  0
        while i < len_lines:
            if i == 2:
                self.ntyp = int(lines[i].split()[0])
                self.natoms = int(lines[i].split()[1])
                self.scale_length = float(lines[i].split()[3])
            if "Basis vectors" in lines[i]:
                for j in range(1,4):
                    self.basis_vec[j-1,:] = numpy.array(list(map(float,lines[i+j].split())))
            if self.start_elem in lines[i]:
                for j in range(self.ntyp): 
                    elem_ = lines[i+j].split("'")[1]
                    elem_ = elem_.split()[0]
                    self.atom_types.append(elem_)
            if i == self.ntyp + 7:
                for j in range(self.natoms):
                    elem = int(lines[i+j].split()[1])
                    #print("elem = ")
                    #print(elem)
                    #self.atom_types_.append(self.atom_types[elem-1])
                    arr = numpy.array(list(map(float,lines[i+j].split()[2:5])))  
                    #self.coordinate.append([elem,arr])
                    self.coordinate.append([elem,arr]) 
            if "Diagonalizing the dynamical matrix" in lines[i]:
                break
            i += 1
        self.lattice = self.scale_length * self.basis_vec
        self.get_elem()
        return

    def merge_coord(self, atom_types, coord):
        xyz = []
        natoms = len(atom_types)
        for i in range(natoms):
            xyz.append([atom_types[i],coord[i]])
        print(xyz)
        return xyz

    def read_lattice_from_inp(self):
        rf = open(self.inputfile,"r")
        lines = rf.readlines()
        self.lattice = numpy.zeros((3,3))
        self.atom_symbols = []
        for i, line in enumerate(lines):
            if "CELL_PARAMETERS bohr" in line:
                for j in range(3):
                    self.lattice[j,:] = numpy.array(list(map(float,lines[i+j+1].split())))
        self.atom_types_ = []
        '''
        for k, v in self.natoms_element.items():
            self.atom_symbols += [k for i in range(v)]
            tmp = [self.atom_number[k] for i in range(v)]
            self.atom_types_ += tmp
        '''
        for i, line in enumerate(lines):
            if "ATOMIC_POSITIONS bohr" in line:
                for j in range(self.natoms):
                    elem = lines[i+j+1].split()[0]
                    self.atom_symbols.append(elem)
                    self.atom_types_.append(self.atom_number[elem])
        return

    def get_elem(self):
        self.atom_symbols = []
        self.atom_types_ = []
        self.coordinate_ = []
        for i,coord in enumerate(self.coordinate):
            elem = self.atom_types[coord[0]-1]
            elem = self.atom_number[elem]
            self.atom_types_.append(elem)
            self.coordinate_.append(coord[1])
        return

    def read_MD_conditions(self):
        filename = self.inputfile
        rf = open(filename,"r")
        lines = rf.readlines()
        #label = lines[0].split()
        for i,line in enumerate(lines):
            label = line.split()[0]
            for k in self.MD_conditions.keys():
                if k == label:
                    self.MD_conditions[k].append(float(line.split()[1]))
        return

    
    def read_MD_coord_(self):
        natoms = self.natoms
        nstep = self.nstep // 10
        self.MD_results["coord"] = [[] for i in range(nstep)]
        filename = self.workdir+"/"+self.prefix+".pos"
        rf = open(filename,"r")
        lines = rf.readlines()
        for i,line in enumerate(lines):
            pos = i % (natoms + 1)
            step = i // (natoms + 1)
            #print("step = ")
            #print(step)
            if pos != 0:
                coord = numpy.array(list(map(float,lines[i].split())))
                self.MD_results["coord"][step].append(coord)
        return

    def read_MD_coord(self):
        natoms = self.natoms
        nstep = self.nstep // 10
        self.MD_results["coord"] = [[] for i in range(nstep)]
        #filename = self.workdir+"/"+self.prefix+".pos"
        rf = open(self.outputfile,"r")
        lines = rf.readlines()
        step = -1
        for i,line in enumerate(lines):
            if "ATOMIC_POSITIONS" in line:
                step += 1
                for j in range(self.natoms):             
                    coord = numpy.array(list(map(float,lines[i+j+1].split()[1:4])))
                    self.MD_results["coord"][step].append(coord)
        return

    def read_MD_results(self):
        filename = self.workdir+"/"+self.prefix+".evp"
        rf = open(filename,"r")
        lines = rf.readlines()
        label = lines[0].split()
        for i,line in enumerate(lines):
            if i > 0:
                for k in self.MD_results.keys():
                    #k = "time(ps)"
                    try:
                        idx = label.index(k)
                    except ValueError:
                        #print("Not Found")
                        continue
                    self.MD_results[label[idx]].append(float(line.split()[idx-1]))
                    #print(float(line.split()[idx]))
        return

    def RMSD(self,larger_displace=10):
        rmsd_list = []
        atoms = []
        norm = []
        vec = []
        coord_0 = self.MD_results["coord"][0]
        for i, coord in enumerate(self.MD_results["coord"]):
            rmsd = 0.0
            for j, atom_coord in enumerate(coord):
                val = atom_coord - coord_0[j]
                val_ = numpy.linalg.norm(val)
                norm.append(val_)
                vec.append(val)
                atoms.append(j)
                rmsd += val_ ** 2
            rmsd = numpy.sqrt(rmsd)
            rmsd /= self.natoms
            rmsd_list.append(rmsd)
        norm = numpy.array(norm)
        vec = numpy.array(vec)
        atoms = numpy.array(atoms)
        sorted_idx = numpy.argsort(norm)
        atoms = atoms[sorted_idx]
        norm = norm[sorted_idx]
        vec = vec[sorted_idx]
        self.large_displace = {"atoms":atoms[0:larger_displace].tolist(),"vec":vec[0:larger_displace].tolist()}
        return rmsd_list # plotclass.plot()


    def choose_atomtype(self,MD_idx):
        if MD_idx is not None:
            self.coordinate = self.MD_results["coord"][MD_idx]
        retval = [[] for i in range(self.ntyp)]
        for i,typ in enumerate(self.atom_types):
            for j,coord in enumerate(self.coordinate):
                if coord[0] == typ:
                    retval[i].append(coord)

        return retval


    def average(self,MD_start,MD_end):
        ret = []
        for j in range(self.natoms):
            coord_ave = numpy.zeros(3)
            for i in range(MD_start, MD_end):
                coord = self.MD_results["coord"][i][j]
                coord_ave += coord
            coord_ave /= (MD_end - MD_start)
            ret.append(coord_ave)
        return ret

    

    def write_xsf(self,filename,MD_idx=None, mode_idx=None,tol=0.2):
        wf = open(filename,"w")
        lattice = self.lattice * self.bohr
        lattice = lattice.T
        wf.write("# \n")
        wf.write("CRYSTAL\n\n")
        wf.write("PRIMVEC\n")
        for i in range(3):
            wf.write(str(lattice[0,i])+" "+str(lattice[1,i])+" "+str(lattice[2,i])+"\n")
        wf.write("CONVVEC\n")
        for i in range(3):
            wf.write(str(lattice[0,i])+" "+str(lattice[1,i])+" "+str(lattice[2,i])+"\n")
        wf.write("PRIMCOORD\n")
        wf.write("{nat} 1 \n".format(nat=self.natoms))
        freq = None
        displace = None
        if MD_idx is not None:
            self.coordinate_ = self.MD_results["coord"][MD_idx]
        if mode_idx is not None:
            freq = self.mode[mode_idx]
            #print(freq)
        if self.large_displace is not None and self.MD:
            displace = self.large_displace
        for i, coord in enumerate(self.coordinate_):
            coord_ = copy.deepcopy(coord) 
            coord_ *= self.bohr
            if self.scale_length is not None:
                coord_ *= self.scale_length
            wf.write(str(self.atom_types_[i])+" "+str(coord_[0])+" "+str(coord_[1])+" "+str(coord_[2])+" ")  
            if freq is not None:
                if numpy.linalg.norm(freq[i]) > tol:
                    wf.write(str(freq[i][0])+" "+str(freq[i][1])+" "+str(freq[i][2]))
                else:
                    wf.write(str(0.0)+" "+str(0.0)+" "+str(0.0))
            if displace is not None:
                if i in displace["atoms"]:
                    idx = displace["atoms"].index(i)
                    freq_ = displace["vec"]
                    wf.write(str(freq_[idx][0])+" "+str(freq_[idx][1])+" "+str(freq_[idx][2]))
                else:
                    wf.write(str(0.0)+" "+str(0.0)+" "+str(0.0))
            wf.write("\n")
        wf.close()
        return

###########################################################################################  

if __name__ == "__main__":
    myclass = Vibration()
    myclass.inputfile = "TMTTF2AsF6_ph.dyn1"
    myclass.read_freq()
    myclass.read_coord()
    print(myclass.atom_types)
    myclass.write_xsf("vibration.xsf",mode_idx=0)
    #myclass.plot_mole()
    #myclass.plot_mole()
    #myclass.read_MD_conditions()
    myclass.read_MD_results()
    myclass.read_MD_coord()
    print("len(coord) = ")
    print(len(myclass.MD_results["coord"]))
    rmsd = myclass.RMSD()
    fig = plt.figure()
    ax_new = fig.add_subplot(1,1,1)
    time = myclass.MD_results["time(ps)"]
    #for i in range(len(time)):
    #    print(time[i])
    ax_new.plot(time,rmsd)
    plt.savefig("rmsd.png")
    #myclass.MD_movies()
    #myclass.movie_MD()
