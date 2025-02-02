wan90_temp0 = """num_bands         =   {nb}
num_wann          =   {nw}
dis_win_max =  {dis_win_max}
dis_win_min =  {dis_win_min}
write_hr = .true.
bands_plot = .true.
bands_num_points = 100
wannier_plot = .true.
"""

wan90_kpath = """
Begin Kpoint_path
L   0.50000     0.00000     0.00000   G   0.00000     0.00000     0.00000
G   0.00000     0.00000     0.00000   X   0.50000     0.50000     0.00000
End Kpoint_path
"""

wan90_temp2 = """mp_grid = {nkx} {nky} {nkz}"""

wan90_pw2wan = """&inputpp
    outdir='{outdir}'
    prefix='{prefix}',
    seedname='{prefix}',
/
"""
