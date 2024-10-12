import pytest
from pylat.crystal_lattice import SpaceGroupSupplement
from pylat.lattice_model import LatticeModel_gen
from pylat.qe_controller import QEController
from pylat.respack import RESPACKController

def test_crystal_lattice():
    cifname = "cif/[3](TMTSF)2AsF6.cif"
    cifout = "cif/[3](TMTSF)2AsF6_sup.cif"
    sgs = SpaceGroupSupplement(cifname)
    sgs.read_file()
    sgs.get_frac()
    sgs.get_xyz()
    sgs.supplement()
    sgs.write_cif(cifout)
    return

def test_gen_lat():
    int_dict = {"t1": 1, "U": 10, "V1": 2}
    lg = LatticeModel_gen(int_dict, 4)
    lg.drive()
    lg.get_ints()

    int_dict = {"t1": 1, "U": 10, "V1": 2}
    lg = LatticeModel_gen(int_dict, lat_type="triangular", H=4, W=4)
    lg.drive()
    lg.get_ints()
    return

def test_qe():
    cifname = "cif/FeSe_mp-20120_primitive.cif"
    pseudo_dict = {"Fe":"Fe.pbe-spn-kjpaw_psl.1.0.0.UPF", "Se":"Se.pbe-dn-kjpaw_psl.1.0.0.UPF"}
    qe = QEController(cifname, pseudo_dict)
    inp = qe.make_input()
    qe.write_input(inp)
    ## qe.exec()

    res = RESPACKController(qe)
    res.prepare()
    return

