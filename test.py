import pytest
from pylat.crystal_lattice import SpaceGroupSupplement
from pylat.lattice_model import LatticeModel_gen
from pylat.qe_controller import QEController
from pylat.respack import RESPACKController

def test_crystal_lattice():
    #cifname = "cif/TMTSF2AsF6.cif"
    #cifout = "cif/TMTSF2AsF6_sup.cif"

    #cifname = "cif/FeSe_mp-20120_primitive.cif"
    #cifout = "cif/FeSe_mp-20120_primitive_out.cif"

    cifname = "cif/Ac2CuGe_sym.cif"
    cifout = "cif/Ac2CuGe_sym_out.cif"
    sgs = SpaceGroupSupplement(cifname)
    sgs.read_file()
    cart = sgs.get_xyz(sgs.xyz)
    sgs.get_frac(cart)
    print(sgs.spacegroup)
    sgs.supplement()
    sgs.write_cif(cifout)
    return

def test_gen_lat():
    int_dict = {"t1": 1, "t2": 0, "U": 10, "V1": 2, "V2": 0}
    lg = LatticeModel_gen(int_dict)
    lg.drive()
    int1e, int2e = lg.get_ints()
    print(int1e[:, 0, :, 0])
    print(int2e)
    assert(int2e[0, 0, 0, 0] == 10)
    assert(int2e[0, 0, 1, 1] == 2)

    int_dict = {"t1": 1, "t2": 0, "U": 10, "V1": 2, "V2": 0}
    lg = LatticeModel_gen(int_dict, lat_type="triangular", H=4, W=4)
    lg.drive()
    int1e, int2e = lg.get_ints()
    return

def test_qe():
    cifname = "cif/FeSe_mp-20120_primitive.cif"
    pseudo_dict = {"Fe":"Fe.pbe-spn-kjpaw_psl.1.0.0.UPF", "Se":"Se.pbe-dn-kjpaw_psl.1.0.0.UPF"}
    qe = QEController(cifname, pseudo_dict)
    inp = qe.make_input()
    qe.write_input(inp)
    ## qe.exec()
    qe.exec_dos()
    qe.exec_pdos()

    res = RESPACKController(qe)
    res.prepare()
    return

if __name__ == "__main__":
    test_qe()

