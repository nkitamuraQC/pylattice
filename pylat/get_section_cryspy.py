

def get_section_for_cryspy(myclass, occ):
    if myclass.calculation == "vc-md":
        if "tetrahedra" in occ:
            section = {
                "&control": [
                    "prefix",
                    "calculation",
                    "outdir",
                    "pseudo_dir",
                    "tstress",
                    "tprnfor",
                    "wf_collect",
                    "dt",
                    "nstep",
                    "restart_mode",
                ],
                "&system": [
                    "ibrav",
                    "nat",
                    "ntyp",
                    "nosym",
                    "ecutwfc",
                    "ecutrho",
                    "occupations",
                ],
                "&electrons": ["conv_thr"],
                "&ions": [],
                "&cell": [],
            }

        if "smearing" in occ:
            section = {
                "&control": [
                    "prefix",
                    "calculation",
                    "outdir",
                    "pseudo_dir",
                    "tstress",
                    "tprnfor",
                    "wf_collect",
                    "dt",
                    "nstep",
                    "restart_mode",
                ],
                "&system": [
                    "ibrav",
                    "nat",
                    "ntyp",
                    "nosym",
                    "ecutwfc",
                    "ecutrho",
                    "occupations",
                    "degauss",
                ],
                "&electrons": ["conv_thr"],
                "&ions": [],
                "&cell": [],
            }

    elif myclass.calculation == "md":
        if "tetrahedra" in occ:
            section = {
                "&control": [
                    "prefix",
                    "calculation",
                    "outdir",
                    "pseudo_dir",
                    "tstress",
                    "tprnfor",
                    "wf_collect",
                    "dt",
                    "nstep",
                    "restart_mode",
                ],
                "&system": [
                    "ibrav",
                    "nat",
                    "ntyp",
                    "nosym",
                    "ecutwfc",
                    "ecutrho",
                    "occupations",
                ],
                "&electrons": ["conv_thr"],
                "&ions": [],
                "&cell": [],
            }

        elif "smearing" in occ:
            section = {
                "&control": [
                    "prefix",
                    "calculation",
                    "outdir",
                    "pseudo_dir",
                    "tstress",
                    "tprnfor",
                    "wf_collect",
                    "dt",
                    "nstep",
                    "restart_mode",
                ],
                "&system": [
                    "ibrav",
                    "nat",
                    "ntyp",
                    "nosym",
                    "ecutwfc",
                    "ecutrho",
                    "occupations",
                    "degauss",
                ],
                "&electrons": ["conv_thr"],
                "&ions": [],
                "&cell": [],
            }
    elif myclass.calculation == "vc-relax":
        if "tetrahedra" in occ:
            section = {
                "&control": [
                    "prefix",
                    "calculation",
                    "outdir",
                    "pseudo_dir",
                    "tstress",
                    "tprnfor",
                    "wf_collect",
                    "dt",
                    "nstep",
                    "restart_mode",
                ],
                "&system": [
                    "ibrav",
                    "nat",
                    "ntyp",
                    "nosym",
                    "ecutwfc",
                    "ecutrho",
                    "occupations",
                ],
                "&electrons": ["conv_thr"],
                "&ions": [],
                "&cell": [],
            }

        if "smearing" in occ:
            section = {
                "&control": [
                    "prefix",
                    "calculation",
                    "outdir",
                    "pseudo_dir",
                    "tstress",
                    "tprnfor",
                    "wf_collect",
                    "dt",
                    "nstep",
                    "restart_mode",
                ],
                "&system": [
                    "ibrav",
                    "nat",
                    "ntyp",
                    "nosym",
                    "ecutwfc",
                    "ecutrho",
                    "occupations",
                    "degauss",
                ],
                "&electrons": ["conv_thr"],
                "&ions": [],
                "&cell": [],
            }

    if myclass.calculation == "relax":
        if "tetrahedra" in occ:
            section = {
                "&control": [
                    "prefix",
                    "calculation",
                    "outdir",
                    "pseudo_dir",
                    "tstress",
                    "tprnfor",
                    "wf_collect",
                    "dt",
                    "nstep",
                    "restart_mode",
                ],
                "&system": [
                    "ibrav",
                    "nat",
                    "ntyp",
                    "nosym",
                    "ecutwfc",
                    "ecutrho",
                    "occupations",
                ],
                "&electrons": ["conv_thr"],
                "&ions": [],
                "&cell": [],
            }

        if "smearing" in occ:
            section = {
                "&control": [
                    "prefix",
                    "calculation",
                    "outdir",
                    "pseudo_dir",
                    "tstress",
                    "tprnfor",
                    "wf_collect",
                    "dt",
                    "nstep",
                    "restart_mode",
                ],
                "&system": [
                    "ibrav",
                    "nat",
                    "ntyp",
                    "nosym",
                    "ecutwfc",
                    "ecutrho",
                    "occupations",
                    "degauss",
                ],
                "&electrons": ["conv_thr"],
                "&ions": [],
                "&cell": [],
            }

    else:
        if "tetrahedra" in occ:
            section = {
                "&control": [
                    "prefix",
                    "calculation",
                    "outdir",
                    "pseudo_dir",
                    "tstress",
                    "tprnfor",
                    "wf_collect",
                    "restart_mode",
                ],
                "&system": [
                    "ibrav",
                    "nat",
                    "ntyp",
                    "ecutwfc",
                    "ecutrho",
                    "occupations",
                    "nosym",
                    "noinv",
                ],
                "&electrons": ["mixing_beta", "conv_thr", "electron_maxstep"],
            }

        elif occ == "smearing":
            section = {
                "&control": [
                    "prefix",
                    "calculation",
                    "outdir",
                    "pseudo_dir",
                    "tstress",
                    "tprnfor",
                    "wf_collect",
                    "restart_mode",
                ],
                "&system": [
                    "ibrav",
                    "nat",
                    "ntyp",
                    "ecutwfc",
                    "ecutrho",
                    "occupations",
                    "nosym",
                    "noinv",
                    "degauss",
                ],
                "&electrons": ["mixing_beta", "conv_thr", "electron_maxstep"],
            }
    return section