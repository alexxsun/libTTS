import libtts
tls_veg_file = "close_stems_3.pts"
th_alpha_sq = 0.01
as_file = libtts.generate_alpha_shape(tls_veg_file, th_alpha_sq)
print("Alpha shape file:", as_file)