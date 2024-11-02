template_res1 = """&param_chiqw
/
&param_wannier
N_wannier={N_wannier}, 
Lower_energy_window={lower_e}, 
Upper_energy_window={upper_e},
N_initial_guess={n_init_guess},
/ 
"""

template_res2 = """&param_interpolation
N_sym_points = 5,!The total number of symmetry points
dense = 16, 16, 16
/
0.50d0 0.50d0 0.50d0 !L
0.00d0 0.00d0 0.00d0 !G
0.50d0 0.00d0 0.50d0 !X
0.50d0 0.25d0 0.75d0 !W
0.50d0 0.50d0 0.50d0 !L
&param_visualization
flg_vis_wannier=1,
!N_write_wannier=5,
ix_vis_min=-1,
ix_vis_max= 1,
iy_vis_min=-1,
iy_vis_max= 1,
iz_vis_min=-1,
iz_vis_max= 1
/
&param_calc_int
/
"""
