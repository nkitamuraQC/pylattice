import pytest
from pylat.crystal_lattice import SpaceGroupSupplement
from pylat.lattice_model import LatticeModel_gen
from pylat.qe_controller import QEController
from pylat.respack import RESPACKController
from pylat.qe_md import QEMD
import numpy as np
from pylat.md_analyzer import MD_Analyzer
from pylat.md_properties import MD_Properties
from pylat.future.wan90 import Wannier90

def _test_crystal_lattice():
    # cifname = "cif/TMTSF2AsF6.cif"
    # cifout = "cif/TMTSF2AsF6_sup.cif"

    # cifname = "cif/FeSe_mp-20120_primitive.cif"
    # cifout = "cif/FeSe_mp-20120_primitive_out.cif"

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


def _test_gen_lat():
    int_dict = {"t1": 1, "t2": 0, "U": 10, "V1": 2, "V2": 0}
    lg = LatticeModel_gen(int_dict)
    lg.drive()
    int1e, int2e = lg.get_ints()
    print(int1e[:, 0, :, 0])
    print(int2e)
    assert int2e[0, 0, 0, 0] == 10
    assert int2e[0, 0, 1, 1] == 2

    int_dict = {"t1": 1, "t2": 0, "U": 10, "V1": 2, "V2": 0}
    lg = LatticeModel_gen(int_dict, lat_type="triangular", H=4, W=4)
    lg.drive()
    int1e, int2e = lg.get_ints()
    return


def test_qe():
    cifname = "cif/FeSe_mp-20120_primitive.cif"
    pseudo_dict = {
        "Fe": "Fe.pbe-spn-kjpaw_psl.1.0.0.UPF",
        "Se": "Se.pbe-dn-kjpaw_psl.1.0.0.UPF",
    }
    qe = QEController(cifname, pseudo_dict)
    inp = qe.make_input()
    qe.write_input(inp)
    qe.exec()
    qe.exec_dos()
    qe.exec_pdos()

    res = RESPACKController(qe)
    res.prepare()
    res.execution()
    return


def test_qe_only():
    cifname = "cif/FeSe_mp-20120_primitive.cif"
    pseudo_dict = {
        "Fe": "Fe.pbe-spn-kjpaw_psl.1.0.0.UPF",
        "Se": "Se.pbe-dn-kjpaw_psl.1.0.0.UPF",
    }
    qe = QEController(cifname, pseudo_dict)
    inp = qe.make_input()
    qe.write_input(inp)
    qe.exec()
    return


def test_qemd():
    cifname = "cif/FeSe_mp-20120_primitive.cif"
    pseudo_dict = {
        "Fe": "Fe.pbe-spn-kjpaw_psl.1.0.0.UPF",
        "Se": "Se.pbe-dn-kjpaw_psl.1.0.0.UPF",
    }
    qe = QEController(cifname, pseudo_dict, supercell=None)
    qe.conv_thr = "1d-5"
    qe.ecutwfc = 60.0
    qe.ecutrho = 240.0
    dL = np.array([[0.05, 0, 0], [0, 0, 0], [0, 0, 0]])
    qemd = QEMD(qe, dL=dL, macro_step=8, temp=10, nkpoints=[8, 8, 8], micro_step=30)
    qemd.run()

    return

def test_analyzer():
    prefix = "calculation_try_0"
    mda = MD_Analyzer(prefix)
    mda.parse()
    mda.plot(target="energies")
    mda.plot_stress(index=[0, 0])
    return


def test_young():
    ndata = 8
    mdas = [MD_Analyzer("calculation_try_{i}".format(i=i)) for i in range(ndata)]
    for mda in mdas:
        mda.parse()
    mdp = MD_Properties(mdas)
    dlat = np.array([[0.05, 0, 0], [0, 0, 0], [0, 0, 0]])
    targets_list = [i for i in range(ndata)]
    mdp.get_young(dlat, targets_list=targets_list)
    mdp.plot()
    return


def test_wan90():
    cifname = "./cif/FeSe_mp-20120_primitive.cif"
    pseudo_dict = {
        "Fe": "Fe.pbe-spn-kjpaw_psl.1.0.0.UPF",
        "Se": "Se.pbe-dn-kjpaw_psl.1.0.0.UPF",
    }
    qe = QEController(cifname, pseudo_dict)
    inp = qe.make_input()
    qe.write_input(inp)
    qe.exec()
    # qe.exec_dos()
    # qe.exec_pdos()
    qe.calculation = "nscf"
    qe.nosym = True
    qe.noinv = True
    qe.exec()
    w90 = Wannier90(qe)

    w90.write_wan90(win_min="1.1049d01", win_max="1.8929d01", nw=5)
    w90.do_wan90()
    return



if __name__ == "__main__":
    #test_qemd()
    #test_qe_only()
    #test_analyzer()
    #test_young()
    test_wan90()
