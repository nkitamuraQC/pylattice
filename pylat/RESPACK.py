#
import numpy,os
from pyscf.pbc.TMTTF import processing
from pyscf.pbc.downfolding import QE as QE_mod
#from pyscf.pbc.analysis import analysis_class

def get_midpoint(xyz):
    sum_ = numpy.zeros(3,dtype=float)
    count = 0
    for i in range(len(xyz)):
        if xyz[i][0] == "S":
            count += 1
            sum_ += xyz[i][1]
    mid = sum_/count
    return mid


def _get_gauss_center_Hchain(RESPACK,QE,num_atom=2):
    xyz = QE.atomic_position
    Tx = QE.T[0,0]
    retval = numpy.zeros((num_atom,3),dtype=float)
    for i in range(num_atom):
        retval[i,0] = xyz[i][1][0]/Tx
    #retval[0,0] = xyz[0][1][0]/Tx
    #retval[1,0] = xyz[1][1][0]/Tx
    return retval

def _get_Bmat_Hchain(RESPACK,QE,num_atom=2):
    RESPACK.Bmat = numpy.identity(num_atom,float)
    return RESPACK.Bmat

def _get_Bmat(RESPACK,QE,num_mole=2):
    RESPACK.Bmat = numpy.zeros((12,num_mole),dtype=float)
    xyz = QE.atomic_position
    if num_mole == 2:
        xyz_1 = xyz[0:26]
        xyz_2 = xyz[26:52]
        xyz_mole = [xyz_1,xyz_2]
    if num_mole == 4:
        xyz_1 = xyz[0:26]
        xyz_2 = xyz[26:52]
        #xyz_3 = xyz[59:85]
        xyz_3 = xyz[52:78]    
        #xyz_4 = xyz[85:111]
        xyz_4 = xyz[78:104]
        xyz_mole = [xyz_1,xyz_2,xyz_3,xyz_4]
    for i in range(num_mole):
        xyz_ = xyz_mole[i]
        vec = processing.calc_plane(xyz_)[0]
        mat = numpy.array([vec[0]**2,vec[1]**2,vec[2]**2])
        RESPACK.Bmat[i*3:(i+1)*3,i] = mat 
    return RESPACK.Bmat

def _get_gauss_center(RESPACK,QE,num_mole=2):
    retval = numpy.zeros((num_mole,3),dtype=float)
    xyz = QE.atomic_position
    if num_mole == 2:
        xyz_1 = xyz[0:26]
        xyz_2 = xyz[26:52]
        xyz_mole = [xyz_1,xyz_2]
    if num_mole == 4:
        xyz_1 = xyz[0:26]
        xyz_2 = xyz[26:52]
        #xyz_3 = xyz[59:85]
        xyz_3 = xyz[52:78]    
        #xyz_4 = xyz[85:111]
        xyz_4 = xyz[78:104]
        xyz_mole = [xyz_1,xyz_2,xyz_3,xyz_4]
    for i in range(num_mole):
        xyz_ = xyz_mole[i]
        coord = get_midpoint(xyz_)
        frac = get_fracCoord(RESPACK,QE,coord)
        retval[i,:] = frac
    return retval

def get_volume(RESPACK,QE):
    lattice = QE.T
    x = lattice[0]
    y = lattice[1]
    z = lattice[2]
    vec = numpy.cross(x,y)
    vol = numpy.dot(vec,z)
    return vol

def get_fracCoord(RESPACK,QE,xyz):
    lattice = QE.T
    norm = []
    v = get_volume(RESPACK,QE)
    cos = numpy.cos
    sin = numpy.sin
    for i in range(3):
        norm_i = numpy.linalg.norm(lattice[i])
        norm.append(norm_i)
    norm = numpy.array(norm)
    alpha_ = numpy.dot(lattice[1],lattice[2])/(norm[1]*norm[2])
    beta_ = numpy.dot(lattice[0],lattice[2])/(norm[0]*norm[2])
    gamma_ = numpy.dot(lattice[0],lattice[1])/(norm[0]*norm[1])
    alpha = numpy.arccos(alpha_)
    beta = numpy.arccos(beta_)
    gamma = numpy.arccos(gamma_)
    tr_mat = numpy.zeros((3,3),dtype=float)
    a = norm[0]
    b = norm[1]
    c = norm[2]
    tr_mat[0,0] = 1.0/a
    tr_mat[0,1] = -cos(gamma)/(a*sin(gamma))
    #tr_mat[0,2] = (cos(alpha)*cos(gamma)-cos(beta))/(a*v*sin(gamma))
    tr_mat[0,2] = b*c*cos(gamma)*(cos(alpha) - cos(beta)*cos(gamma))/(v*sin(gamma)) - b*c*cos(beta)*sin(gamma)/v
    tr_mat[1,1] = 1/(b*sin(gamma))
    tr_mat[1,2] = a*c*(cos(beta)*cos(gamma)-cos(alpha))/(v*sin(gamma))
    tr_mat[2,2] = a*b*sin(gamma)/v
    frac = numpy.dot(tr_mat,xyz)
    return frac

def _write_inp(RESPACK):
    filename = RESPACK.workpath + "/" +RESPACK.filename
    wf = open(filename,"w")
    wf.write("&param_chiqw"+"\n")
    wf.write("flg_cRPA={rpa}".format(rpa=RESPACK.flg_cRPA)+", \n")
    wf.write("MPI_num_qcomm={mpi}".format(mpi=RESPACK.MPI_num_qcomm)+"/ \n")
    wf.write("\n")
    wf.write("&param_wannier"+"\n")
    wf.write("N_wannier={nwan}".format(nwan = RESPACK.N_wannier)+", \n")
    wf.write("Lower_energy_window={low}".format(low=RESPACK.Lower_energy_window)+", \n")
    wf.write("Upper_energy_window={up}".format(up=RESPACK.Upper_energy_window)+", \n")
    wf.write("FLG_BMAT={flag}".format(flag=RESPACK.FLG_BMAT)+", \n")
    wf.write("N_initial_guess={Nguess}".format(Nguess=RESPACK.N_initial_guess)+"/ \n")
    # gaussian orbital
    print(RESPACK.gauss_center)
    for i in range(RESPACK.N_initial_guess):
        basis = RESPACK.N_initial_guess // RESPACK.N_wannier
        idx = i // basis
        type_orb = RESPACK.gaussian_orb[0][i]
        exp_orb = RESPACK.gaussian_orb[1][i]
        coord_x = RESPACK.gauss_center[idx,0]
        coord_y = RESPACK.gauss_center[idx,1]
        coord_z = RESPACK.gauss_center[idx,2]
        wf.write(type_orb + " " + str(exp_orb) + " " + str(coord_x) \
                 + " " + str(coord_y) + " " + str(coord_z) + "\n")

    # Bmat
    for i in range(RESPACK.N_initial_guess):
        for j in range(RESPACK.N_wannier):
            wf.write(str(RESPACK.Bmat[i,j]) + " ")
        wf.write("\n")

    wf.write("&param_interpolation"+"\n")
    wf.write("N_sym_points={npoint}".format(npoint=RESPACK.N_sym_points)+"/ \n")
    bandpath = RESPACK.bandpath
    for i in range(RESPACK.N_sym_points):
        for j in range(3):
            wf.write(str(bandpath[i][j])+" ")
        wf.write("\n")

    wf.write("&param_visualization"+"\n")
    wf.write("/"+"\n")

    wf.write("&param_calc_int"+"\n")
    wf.write("/"+"\n")
    return

def _execution(RESPACK,QE):
    #util = RESPACK.RESPACK_path + '/' + 'util/qe2respack'
    util = RESPACK.qe2respack_path
    work = RESPACK.workpath
    filename = work + "/" + RESPACK.filename
    dir_ = QE.prefix + ".save"
    #os.system("python {util}/qe2respack.py {work}".format(util=util,work=work))
    #os.system("bash {util}/qe2respack.sh {work}".format(util=util,work=work))
    if RESPACK.calctype == "calc_wannier":
        os.system("python {util}/qe2respack.py {work}/{dir_}".format(util=util,work=work,dir_=dir_))
        os.system("calc_wannier < {filename} > {filename}_wan.out".format(filename=filename))
    elif RESPACK.calctype == "calc_chiqw":
        if RESPACK.MPI > 1:
            os.system("OMP_NUM_THREADS=1 mpirun -np {MPI} calc_chiqw < {filename} > {filename}_chi.out".format(filename=filename,MPI=RESPACK.MPI))
        else:
            os.system("calc_chiqw < {filename} > {filename}_chi.out".format(filename=filename,MPI=RESPACK.MPI))
    elif RESPACK.calctype == "calc_w3d":
        os.system("calc_w3d < {filename} > {filename}_w.out".format(filename=filename))
    elif RESPACK.calctype == "calc_j3d":
        os.system("calc_j3d < {filename} > {filename}_j.out".format(filename=filename))
    return

def _exec(RESPACK,QE):
    #file_bat = mVMC.work_dir+"/"+mVMC.project_idx+mVMC.projectname+".bat"
    file_bat = RESPACK.workpath+"/"+RESPACK.filename+".bat" 
    work = RESPACK.workpath
    #filename_ = work + "/" + RESPACK.filename
    filename = RESPACK.filename
    #util = RESPACK_path + '/' + 'util/qe2respack'
    util = RESPACK.qe2respack_path
    dir_ = QE.prefix + ".save"
    os.system("python {util}/qe2respack.py {work}/{dir_}".format(util=util,work=work,dir_=dir_))
    wf = open(file_bat,"w")
    wf.write("#!/bin/sh"+"\n")
    wf.write("#PBS -m be"+"\n")
    wf.write("#PBS -l select=1:ncpus={a}:mpiprocs={b}:ompthreads={c}:jobtype={d}".format(a=RESPACK.QE.IMS.ncpus,\
                                                                                         b=RESPACK.QE.IMS.mpiprocs,\
                                                                                         c=RESPACK.QE.IMS.ompthreads,\
                                                                                         d=RESPACK.QE.IMS.jobtype)+"\n")
    wf.write("#PBS -l walltime={e}".format(e=RESPACK.QE.IMS.walltime)+"\n")
    wf.write("if [ ! -z ${PBS_O_WORKDIR}]; then"+"\n")
    wf.write("     cd ${PBS_O_WORKDIR}"+"\n")
    wf.write("fi"+"\n")
    wf.write("bash"+"\n")
    wf.write("source ~/.bash_profile"+"\n")
    wf.write("export HomeDir={path}".format(path=RESPACK.workpath)+"\n")
    wf.write("export WorkDir=/work/users/$USER/${PBS_JOBID}"+"\n")
    wf.write("mkdir -p ${WorkDir}"+"\n")
    wf.write("cd ${WorkDir}"+"\n")
    wf.write("cp $HomeDir/* ."+"\n")
    # modify here
    wf.write("calc_wannier < $WorkDir/{filename} > $HomeDir/{filename}_wan.out".format(filename=filename)+"\n")
    if RESPACK.QE.MPI < 2 or RESPACK.QE.MPI is None:
        wf.write("calc_chiqw < $WorkDir/{filename} > $HomeDir/{filename}_chi.out".format(filename=filename)+"\n")
    else:
        wf.write("OMP_NUM_THREADS=1 mpirun -np {MPI} calc_chiqw < $WorkDir/{filename} > $HomeDir/{filename}_chi.out".format(MPI=RESPACK.MPI,filename=filename)+"\n")
    wf.write("calc_w3d < $WorkDir/{filename} > $HomeDir/{filename}_w.out".format(filename=filename)+"\n")
    wf.write("calc_j3d < $WorkDir/{filename} > $HomeDir/{filename}_j.out".format(filename=filename)+"\n")
    wf.write("mkdir ${HomeDir}/TMTTF_0_output"+"\n")
    wf.write("cp $WorkDir/* ${HomeDir}/TMTTF_0_output"+"\n")
    wf.close()
    filename_sh = file_bat[0:-3] + ".sh"
    wf = open(filename_sh,"w")
    wf.write("#!/bin/bash" + "\n")
    wf.write("jsub -q PN {fn}".format(fn = file_bat)+ "\n")
    wf.close()
    if not RESPACK.debug:
        os.system("bash {sh}".format(sh = filename_sh)+"\n")


def _kernel(RESPACK,QE):
    if QE.flag_Hchain:
        gauss_center = RESPACK.get_gauss_center_Hchain()
        RESPACK.Gauss_Hchain(gauss_center)
    else:
        gauss_center = RESPACK.get_gauss_center()
        RESPACK.Gauss(gauss_center)
    if RESPACK.FLG_BMAT == 1:
        if QE.flag_Hchain:
            RESPACK.get_Bmat_Hchain()
        else:
            RESPACK.get_Bmat()
    RESPACK.get_bandpath()
    RESPACK.write_inp()
    #self.execution()
    return

class RESPACK:
    def __init__(self,QE,path,filename,calc,computer):    
        self.debug = True
        self.workpath = path
        self.qe2respack_path = "/home/users/auv/pyscf/pyscf/pbc/downfolding"
        self.filename = filename
        self.calctype = calc
        self.flg_cRPA = 1
        self.MPI_num_qcomm = 12
        self.N_wannier = 2
        #self.Lower_energy_window = "2.7d0" ## TMTSF
        #self.Upper_energy_window = "4.5d0" ## TMTSF
        self.Lower_energy_window = "-10.0d0"
        self.Upper_energy_window = "10.0d0"
        self.FLG_BMAT = 1
        self.N_initial_guess = 2 # 6
        self.gaussian_orb = None
        self.Bmat = None
        self.N_sym_points = 5
        self.bandpath = None
        if computer == "qcl":
            self.RESPACK_path = "/home/nkitamura/respack/RESPACK-20190527-dist/"
        if computer == "ims":
            self.RESPACK_PATH = ""
        self.MPI = 1
        self.QE = QE
        self.gauss_center = None
        self.flag_Hchain = True
        self.num_atom = None
        self.num_mole = None
        self.computer = computer 
        return

    def Gauss_Hchain(self,gauss_center):
        num_atom = self.num_atom
        #gauss_type = ["s","s"]
        gauss_type = ["s" for i in range(num_atom)]
        #gauss_exp = [0.20,0.20]
        gauss_exp = [0.20 for i in range(num_atom)]
        gauss_center = gauss_center
        gauss = [gauss_type,gauss_exp]
        self.gaussian_orb = gauss
        self.gauss_center = gauss_center
        return
        
    def Gauss(self,gauss_center):
        gauss_type = ["px","py","pz","px","py","pz"]
        gauss_exp = [0.20,0.20,0.20,0.20,0.20,0.20]
        gauss_center = gauss_center
        gauss = [gauss_type,gauss_exp]
        self.gaussian_orb = gauss
        self.gauss_center = gauss_center
        return

    def get_gauss_center(self):
        return _get_gauss_center(self,self.QE)


    def get_gauss_center_Hchain(self):
        return _get_gauss_center_Hchain(self,self.QE,num_atom=self.num_atom)

    
    def get_Bmat(self):
        # self.BMAT = BMAT  
        return _get_Bmat(self,self.QE)

    def get_Bmat_Hchain(self):
        # self.BMAT = BMAT  
        return _get_Bmat_Hchain(self,self.QE,num_atom=self.num_atom)

    def get_bandpath(self):
        # TMTTF
        bandpath = [[0.000,0.500,0.000],[0.000,0.000,0.000],[0.500,0.000,0.000],\
                    [0.500,0.500,0.000],[0.000,0.000,0.000]]
        self.bandpath = bandpath
        assert(self.N_sym_points == len(bandpath))
        return

    def kernel(self):
        return _kernel(self,self.QE)

    def write_inp(self):
        return _write_inp(self)

    def execution(self,QE):
        if self.computer == "ims":
            return _exec(self,QE)
        else:
            return _execution(self,QE)
    
if __name__ == "__main__":
    filename = "TMTTF"
    num_atoms = 4
    path = os.getcwd()
    prefix = "TMTTF2AsF6"
    calc = "scf"
    qe = QE_mod.QuantumEspresso(filename,path,prefix,calc)
    qe.nat = num_atoms
    qe.ntyp = 1
    qe.ibrav = 0
    qe.atoms = QE_mod.Atoms(path)
    qe.atoms.H_dist = 0.4
    qe.atoms.num_atoms = num_atoms
    qe.pseudo_dir = ""
    qe.MPI = 10
    qe.kernel()
    #QE.QE_exec()
    filename_res = "TMTTF_res"
    calc = "calc_wannier"
    respack = RESPACK(qe,path,filename_res,calc,"qcl")
    respack.num_atom = num_atoms
    respack.N_wannier = num_atoms
    respack.N_initial_guess = num_atoms
    respack.kernel()
    #respack.execution(computer="qcl")
