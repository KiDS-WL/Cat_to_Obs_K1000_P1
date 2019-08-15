README for numida
=================
Chieh-An Lin
Version 2019.05.07


Cosmological parameters
-----------------------

  - class/parameters.ini    | Parameters read by Class
  - class/unused_parameters | Parameters unread by Class

### Class input parameters

  - h         = 0.6898
  - Omega_b   = 0.0473
  - N_ur      = 3.046
  - Omega_cdm = 0.2432
  - Omega_k   = 0.
  - w0_fld    = -1.0
  - wa_fld    = 0.
  - z_reio    = 10.5
  - A_s       = 2.190290e-09
  - n_s       = 0.969

### Derived parameters

  - Omega_m   = 0.2905
  - sigma_8   = 0.826


Redshift distribution
---------------------
  
  - lens/nOfZ_hist_BOSSA_tomo0.dat     | n(z) as histogram for bin 1 (lens bin 1)
  - lens/nOfZ_hist_BOSSA_tomo1.dat     | n(z) as histogram for bin 2 (lens bin 2)
  - source/nOfZ_hist_KiDSVDA_tomo0.dat | n(z) as histogram for bin 3 (source bin 1)
  - source/nOfZ_hist_KiDSVDA_tomo1.dat | n(z) as histogram for bin 4 (source bin 2)
  - source/nOfZ_hist_KiDSVDA_tomo2.dat | n(z) as histogram for bin 5 (source bin 3)
  - source/nOfZ_hist_KiDSVDA_tomo3.dat | n(z) as histogram for bin 6 (source bin 4)
  - source/nOfZ_hist_KiDSVDA_tomo4.dat | n(z) as histogram for bin 7 (source bin 5)
  
### Columns

  0 = lower edge of the redshift bin
  1 = normalized number counts
  
  
Galaxy bias
-----------

  - b_gal = 2.0 (optional)

  
Linc's wish list
----------------
  
I wish to obtain:
  - All 3x2pt mean signal prediction:
    - w: all lens pairs (11, 12, 22)
    - gamma_t: all lens x source pairs (13, 14, ..., 17, 23, 24, ..., 27)
    - xi_+/-: all source pairs (33, 34, ..., 77)
  - C_l:
    - delta-delta: all pairs
    - delta-kappa: all full x source pairs (13, 14, ..., 17, 23, 24, ..., 27, 33, 34, ... 77)
    - kappa-kappa: all source pairs (33, 34, ..., 77)
    - This is basically the upper triangle of a sqaured matrix of which the side is (delta_1, ..., delta_7, kappa_3, ... kappa_7). 
Thanks!

