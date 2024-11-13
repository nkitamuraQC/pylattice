wan90_temp0 = """num_bands         =   {nb}
num_wann          =   {nw}
exclude_bands={a}-{b}
"""

wan90_temp1 = """
write_hr = .true.
bands_plot = .true.

begin kpoint_path
L   0.50000     0.00000     0.00000   G   0.00000     0.00000     0.00000
G   0.00000     0.00000     0.00000   X   0.50000     0.50000     0.00000
end kpoint_path
"""

wan90_temp2 = """
mp_grid = {nkx} {nky} {nkz}
"""

wan90_pw2wan = """&inputpp
    outdir='{outdir}'
    prefix='{prefix}',
    seedname = '{prefix}'
/
"""
