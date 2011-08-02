# Copyright (C) 2008 by Jes Frellsen, Ida Moltke and Martin Thiim
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

from math import pi
from Bio.PDB import Vector

def deg2rad(deg):
    return deg * (pi/180)

def rad2deg(rad):
    return rad * (180/pi)


## Distances (all in amino acids)

# Distances in backbone
d_PO5p   = 1.593                  # std. 0.010
d_O5pC5p = 1.440                  # std. 0.016 
d_C5pC4p = 1.510                  # std. 0.013
d_C4pC3p = 1.524                  # std. 0.011
d_C3pO3p = 1.423                  # std. 0.014
d_O3pP   = 1.607                  # std. 0.012

d_PO1P   = 1.485                  # std. 0.017
d_PO2P   = 1.485                  # std. 0.017


# Distances in sugar
# General
d_C1pC2p = 1.528                  # std. 0.010
d_C2pC3p = 1.525                  # std. 0.011
d_C3pC4p = 1.524                  # std. 0.011
d_C4pO4p = 1.453                  # std. 0.012
d_O4pC1p = 1.414                  # std. 0.012
d_C2pO2p = 1.413                  # std. 0.013
d_C1pN   = 1.471                  # std. 0.017  


# C2 endo specific
d_2e_C1pC2p = 1.526                  # std. 0.008
d_2e_C2pC3p = 1.525                  # std. 0.011
d_2e_C3pC4p = 1.527                  # std. 0.011
d_2e_C4pO4p = 1.454                  # std. 0.010
d_2e_O4pC1p = 1.415                  # std. 0.012
d_2e_C3pO3p = 1.427                  # std. 0.012
d_2e_C5pC4p = 1.509                  # std. 0.012
d_2e_C2pO2p = 1.412                  # std. 0.013
d_2e_C1pN   = 1.464                  # std. 0.014
d_2e_O5pC5p = 1.424                  # std. 0.016


# C3 endo specific
d_3e_C1pC2p = 1.529                  # std. 0.011
d_3e_C2pC3p = 1.523                  # std. 0.011
d_3e_C3pC4p = 1.521                  # std. 0.010
d_3e_C4pO4p = 1.451                  # std. 0.013
d_3e_O4pC1p = 1.412                  # std. 0.013
d_3e_C3pO3p = 1.417                  # std. 0.014
d_3e_C5pC4p = 1.508                  # std. 0.007
d_3e_C2pO2p = 1.420                  # std. 0.010
d_3e_C1pN   = 1.483                  # std. 0.015
d_3e_O5pC5p = 1.420                  # std. 0.009


# Distances in bases



# In base Adenine
d_A_N1C2 = 1.339                  # std. 0.009
d_A_C2N3 = 1.331                  # std. 0.009
d_A_N3C4 = 1.344                  # std. 0.006
d_A_C4C5 = 1.383                  # std. 0.007
d_A_C5C6 = 1.406                  # std. 0.009
d_A_C6N1 = 1.351                  # std. 0.007
d_A_C5N7 = 1.388                  # std. 0.006
d_A_N7C8 = 1.311                  # std. 0.007
d_A_C8N9 = 1.373                  # std. 0.008
d_A_N9C4 = 1.374209591    # 1.374 # std. 0.006
d_A_C6N6 = 1.335                  # std. 0.008
d_A_N9C1p= 1.462                  # std. 0.010


# In base Cytosine
d_C_N1C2 = 1.3965761705   # 1.397 # std. 0.010
d_C_C2N3 = 1.353                  # std. 0.008
d_C_N3C4 = 1.335                  # std. 0.007
d_C_C4C5 = 1.425                  # std. 0.008
d_C_C5C6 = 1.339                  # std. 0.008
d_C_C6N1 = 1.367                  # std. 0.006
d_C_C2O2 = 1.240                  # std. 0.009
d_C_C4N4 = 1.335                  # std. 0.009
d_C_N1C1p= 1.470                  # std. 0.012


# In base Guanine
d_G_N1C2 = 1.373                  # std. 0.008
d_G_C2N3 = 1.323                  # std. 0.008
d_G_N3C4 = 1.350                  # std. 0.007
d_G_C4C5 = 1.379                  # std. 0.007
d_G_C5C6 = 1.419                  # std. 0.010
d_G_C6N1 = 1.391                  # std. 0.007
d_G_C5N7 = 1.388                  # std. 0.006
d_G_N7C8 = 1.305                  # std. 0.006
d_G_C8N9 = 1.374                  # std. 0.007
d_G_N9C4 = 1.375420951    # 1.375 # std. 0.008
d_G_C2N2 = 1.341                  # std. 0.010
d_G_C6O6 = 1.237                  # std. 0.009
d_G_N9C1p = 1.459                 # std. 0.009


# In base Uracil
d_U_N1C2 = 1.38052345145  # 1.381 # std. 0.009
s_U_C2N3 = 1.373                  # std. 0.007
d_U_N3C4 = 1.380                  # std. 0.009
d_U_C4C5 = 1.431                  # std. 0.009
d_U_C5C6 = 1.337                  # std. 0.009
d_U_C6N1 = 1.375                  # std. 0.009
d_U_C2O2 = 1.219                  # std. 0.009
d_U_C4O4 = 1.232                  # std. 0.008
d_U_N1C1p = 1.469                 # std. 0.014



## Bond angles (all in radians)


# In backbone
a_O3pPO5p   = deg2rad( 104.0 )    # std. 1.9
a_PO5pC5p   = deg2rad( 120.9 )    # std. 1.6
a_C5pC4pC3p = deg2rad( 115.5 )    # std. 1.5
a_C4pC3pO3p = deg2rad( 110.6 )    # std. 2.6
a_C3pO3pP   = deg2rad( 119.7 )    # std. 1.2

a_O5pPO1P   = deg2rad( 108.1 )    # std. 2.9
a_O5pPO2P   = deg2rad( 108.3 )    # std. 2.7
a_O3pPO1P   = deg2rad( 107.4 )    # std. 3.2
a_O3pPO2P   = deg2rad( 108.3 )    # std. 3.2


# In sugar
a_C1pC2pC3p = deg2rad( 101.5 )    # std. 0.9
a_C2pC3pC4p = deg2rad( 102.7 )    # std. 1.0
a_C3pC4pO4p = deg2rad( 105.5 )    # std. 1.4
a_C4pO4pC1p = deg2rad( 109.6 )    # std. 0.9
a_O4pC1pC2p = deg2rad( 106.4 )    # std. 1.4
a_C1pC2pO2p = deg2rad( 110.6 )    # std. 3.0
a_C3pC2pO2p = deg2rad( 113.3 )    # std. 2.9
a_C2pC3pO3p = deg2rad( 111.0 )    # std. 2.8
a_C5pC4pO4p = deg2rad( 109.2 )    # std. 1.4
a_O4pC1pN   = deg2rad( 108.2 )    # std. 1.0
a_C2pC1pN   = deg2rad( 113.4 )    # std. 1.6
a_O5pC5pC4p = deg2rad( 111.6 )    # std. 1.7
a_C1pN9C4   = deg2rad( 127.1 )    # std. 1.8
a_C1pN1C2   = deg2rad( 117.9 )    # std. 1.3


# C2 endo specific
a_2e_C1pC2pC3p = deg2rad( 101.5 )    # std. 0.8
a_2e_C2pC3pC4p = deg2rad( 102.6 )    # std. 1.0
a_2e_C3pC4pO4p = deg2rad( 106.1 )    # std. 0.8
a_2e_C4pO4pC1p = deg2rad( 109.7 )    # std. 0.7
a_2e_O4pC1pC2p = deg2rad( 105.8 )    # std. 1.0
a_2e_C1pC2pO2p = deg2rad( 111.8 )    # std. 2.6
a_2e_C3pC2pO2p = deg2rad( 114.6 )    # std. 2.2
a_2e_C2pC3pO3p = deg2rad( 109.5 )    # std. 2.2
a_2e_C4pC3pO3p = deg2rad( 109.4 )    # std. 2.1
a_2e_C5pC4pC3p = deg2rad( 115.2 )    # std. 1.4
a_2e_C5pC4pO4p = deg2rad( 109.1 )    # std. 1.2
a_2e_O4pC1pN   = deg2rad( 108.2 )    # std. 0.8
a_2e_C2pC1pN   = deg2rad( 114.0 )    # std. 1.3
a_2e_O5pC5pC4p = deg2rad( 111.7 )    # std. 1.9
a_2e_C1pN9C4   = deg2rad( 127.4 )    # std. 1.2
a_2e_C1pN1C2   = deg2rad( 118.5 )    # std. 1.1


# C3 endo specific
a_3e_C1pC2pC3p = deg2rad( 101.3 )    # std. 0.7
a_3e_C2pC3pC4p = deg2rad( 102.6 )    # std. 1.0
a_3e_C3pC4pO4p = deg2rad( 104.0 )    # std. 1.0
a_3e_C4pO4pC1p = deg2rad( 109.9 )    # std. 0.8
a_3e_O4pC1pC2p = deg2rad( 107.6 )    # std. 0.9
a_3e_C1pC2pO2p = deg2rad( 108.4 )    # std. 2.4
a_3e_C3pC2pO2p = deg2rad( 110.7 )    # std. 2.1
a_3e_C2pC3pO3p = deg2rad( 113.7 )    # std. 1.6
a_3e_C4pC3pO3p = deg2rad( 113.0 )    # std. 2.0
a_3e_C5pC4pC3p = deg2rad( 116.0 )    # std. 1.6
a_3e_C5pC4pO4p = deg2rad( 109.8 )    # std. 0.9
a_3e_O4pC1pN   = deg2rad( 108.5 )    # std. 0.7
a_3e_C2pC1pN   = deg2rad( 112.0 )    # std. 1.1
a_3e_O5pC5pC4p = deg2rad( 111.5 )    # std. 1.6
a_3e_C1pN9C4   = deg2rad( 126.3 )    # std. 2.8
a_3e_C1pN1C2   = deg2rad( 116.7 )    # std. 0.6


# In bases

# In base Adenine 
a_A_C6N1C2  = deg2rad( 118.6 )    # std. 0.6
a_A_N1C2N3  = deg2rad( 129.3 )    # std. 0.5
a_A_C2N3C4  = deg2rad( 110.6 )    # std. 0.5
a_A_N3C4C5  = deg2rad( 126.8 )    # std. 0.7
a_A_C4C5C6  = deg2rad( 117.0 )    # std. 0.5
a_A_C5C6N1  = deg2rad( 117.7 )    # std. 0.5
a_A_C4C5N7  = deg2rad( 110.7 )    # std. 0.5
a_A_C5N7C8  = deg2rad( 103.9 )    # std. 0.5
a_A_N7C8N9  = deg2rad( 113.8 )    # std. 0.5
a_A_C8N9C4  = deg2rad( 105.8 )    # std. 0.4
a_A_N9C4C5  = deg2rad( 105.8 )    # std. 0.4
a_A_N3C4N9  = deg2rad( 127.4 )    # std. 0.8
a_A_C6C5N7  = deg2rad( 132.3 )    # std. 0.7
a_A_N1C6N6  = deg2rad( 118.6 )    # std. 0.6
a_A_C5C6N6  = deg2rad( 123.7 )    # std. 0.8
a_A_C8N9C1p = deg2rad( 127.7 )    # std. 1.8
a_A_C4N9C1p = deg2rad( 126.3 )    # std. 1.8


# In base Cytosine
a_C_C6N1C2  = deg2rad( 120.3 )    # std. 0.4
a_C_N1C2N3  = deg2rad( 119.2 )    # std. 0.7
a_C_C2N3C4  = deg2rad( 119.9 )    # std. 0.5
a_C_N3C4C5  = deg2rad( 121.9 )    # std. 0.4
a_C_C4C5C6  = deg2rad( 117.4 )    # std. 0.5
a_C_C5C6N1  = deg2rad( 121.0 )    # std. 0.5
a_C_N1C2O2  = deg2rad( 118.9 )    # std. 0.6
a_C_N3C2O2  = deg2rad( 121.9 )    # std. 0.7
a_C_N3C4N4  = deg2rad( 118.0 )    # std. 0.7
a_C_C5C4N4  = deg2rad( 120.2 )    # std. 0.7
a_C_C6N1C1p = deg2rad( 120.8 )    # std. 1.2
a_C_C2N1C1p = deg2rad( 118.8 )    # std. 1.1


# In base Guanine
a_G_C6N1C2  = deg2rad( 125.1 )    # std. 0.6
a_G_N1C2N3  = deg2rad( 123.9 )    # std. 0.6
a_G_C2N3C4  = deg2rad( 111.9 )    # std. 0.5
a_G_N3C4C5  = deg2rad( 128.6 )    # std. 0.5
a_G_C4C5C6  = deg2rad( 118.8 )    # std. 0.6
a_G_C5C6N1  = deg2rad( 111.5 )    # std. 0.5
a_G_C4C5N7  = deg2rad( 110.8 )    # std. 0.4
a_G_C5N7C8  = deg2rad( 104.3 )    # std. 0.5
a_G_N7C8N9  = deg2rad( 113.1 )    # std. 0.5
a_G_C8N9C4  = deg2rad( 106.4 )    # std. 0.4
a_G_N9C4C5  = deg2rad( 105.4 )    # std. 0.4
a_G_N3C4N9  = deg2rad( 126.0 )    # std. 0.6
a_G_C6C5N7  = deg2rad( 130.4 )    # std. 0.6
a_G_N1C2N2  = deg2rad( 116.2 )    # std. 0.9
a_G_N3C2N2  = deg2rad( 119.9 )    # std. 0.7
a_G_N1C6O6  = deg2rad( 119.9 )    # std. 0.6
a_G_C5C6O6  = deg2rad( 128.6 )    # std. 0.6
a_G_C8N9C1p = deg2rad( 127.0 )    # std. 1.3
a_G_C4N9C1p = deg2rad( 126.5 )    # std. 1.3


# In base Uracil
a_U_C6N1C2  = deg2rad( 121.0 )    # std. 0.6
a_U_N1C2N3  = deg2rad( 114.9 )    # std. 0.6
a_U_C2N3C4  = deg2rad( 127.0 )    # std. 0.6
a_U_N3C4C5  = deg2rad( 114.6 )    # std. 0.6
a_U_C4C5C6  = deg2rad( 119.7 )    # std. 0.6
a_U_C5C6N1  = deg2rad( 122.7 )    # std. 0.5
a_U_N1C2O2  = deg2rad( 122.8 )    # std. 0.7
a_U_N3C2O2  = deg2rad( 122.2 )    # std. 0.7
a_U_N3C4O4  = deg2rad( 119.4 )    # std. 0.7
a_U_C5C4O4  = deg2rad( 125.9 )    # std. 0.6
a_U_C6N1C1p = deg2rad( 121.2 )    # std. 1.4
a_U_C2N1C1p = deg2rad( 117.7 )    # std. 1.2




## Reference frame coordinates for the bases
# 


# Adenine
A_C1p = Vector(-2.479, 5.346, 0.000)
A_N9  = Vector(-1.291, 4.498, 0.000)
A_C8  = Vector( 0.024, 4.897, 0.000)
A_N7  = Vector( 0.877, 3.902, 0.000)
A_C5  = Vector( 0.071, 2.771, 0.000)
A_C6  = Vector( 0.369, 1.398, 0.000)
A_N6  = Vector( 1.611, 0.909, 0.000)
A_N1  = Vector(-0.668, 0.532, 0.000)
A_C2  = Vector(-1.912, 1.023, 0.000)
A_N3  = Vector(-2.320, 2.290, 0.000)
A_C4  = Vector(-1.267, 3.124, 0.000)

A_coords = [A_C1p,A_N9,A_C8,A_N7,A_C5,A_C6,A_N6,A_N1,A_C2,A_N3,A_C4]
A_names  = ["C1'","N9", "C8", "N7", "C5", "C6", "N6", "N1" , "C2", "N3", "C4"]
A_dict   = dict(zip(A_names, A_coords))

# Cytosine
C_C1p = Vector(-2.477, 5.402, 0.000)
C_N1  = Vector(-1.285, 4.542, 0.000)
C_C2  = Vector(-1.472, 3.158, 0.000)
C_O2  = Vector(-2.628, 2.709, 0.000)
C_N3  = Vector(-0.391, 2.344, 0.000)
C_C4  = Vector( 0.837, 2.868, 0.000)
C_N4  = Vector( 1.875, 2.027, 0.000)
C_C5  = Vector( 1.056, 4.275, 0.000)
C_C6  = Vector(-0.023, 5.068, 0.000)

C_coords = [C_C1p,C_N1,C_C2,C_O2,C_N3,C_C4,C_N4,C_C5,C_C6]
C_names  = ["C1'", "N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"]
C_dict   = dict(zip(C_names, C_coords))

# Guanine
G_C1p = Vector(-2.477, 5.399, 0.000)
G_N9  = Vector(-1.289, 4.551, 0.000)
G_C8  = Vector( 0.023, 4.962, 0.000)
G_N7  = Vector( 0.870, 3.969, 0.000)
G_C5  = Vector( 0.071, 2.833, 0.000)
G_C6  = Vector( 0.424, 1.460, 0.000)
G_O6  = Vector( 1.554, 0.955, 0.000)
G_N1  = Vector(-0.700, 0.641, 0.000)
G_C2  = Vector(-1.999, 1.087, 0.000)
G_N2  = Vector(-2.949, 0.139,-0.001)
G_N3  = Vector(-2.342, 2.364, 0.001)
G_C4  = Vector(-1.265, 3.177, 0.000)

G_coords = [G_C1p,G_N9,G_C8,G_N7,G_C5,G_C6,G_O6,G_N1,G_C2,G_N2,G_N3,G_C4]
G_names  = ["C1'","N9","C8","N7","C5","C6","O6","N1","C2","N2","N3","C4"]
G_dict   = dict(zip(G_names, G_coords))

# Base Uracil
U_C1p = Vector(-2.481, 5.354, 0.000)
U_N1  = Vector(-1.284, 4.500, 0.000)
U_C2  = Vector(-1.462, 3.131, 0.000)
U_O2  = Vector(-2.563, 2.608, 0.000)
U_N3  = Vector(-0.302, 2.397, 0.000)
U_C4  = Vector( 0.989, 2.884, 0.000)
U_O4  = Vector( 1.935, 2.094,-0.001)
U_C5  = Vector( 1.089, 4.311, 0.000)
U_C6  = Vector(-0.024, 5.053, 0.000)

U_coords = [U_C1p,U_N1,U_C2,U_O2,U_N3,U_C4,U_O4,U_C5,U_C6]
U_names  = ["C1'","N1","C2","O2","N3","C4","O4","C5","C6"]
U_dict   = dict(zip(U_names, U_coords))


bb_distances = (d_PO5p, d_O5pC5p, d_C5pC4p, d_C4pC3p, d_C3pO3p, d_O3pP)
bb_angles = (a_O3pPO5p, a_PO5pC5p, a_O5pC5pC4p, a_C5pC4pC3p, a_C4pC3pO3p, a_C3pO3pP)
bb_atomnames = ["P", "O5'", "C5'", "C4'", "C3'", "O3'"]

bb_distances_2e = (d_PO5p, d_O5pC5p, d_C5pC4p, d_2e_C3pC4p, d_C3pO3p, d_O3pP)
bb_angles_2e = (a_O3pPO5p, a_PO5pC5p, a_O5pC5pC4p, a_2e_C5pC4pC3p, a_2e_C4pC3pO3p, a_C3pO3pP)

bb_distances_3e = (d_PO5p, d_O5pC5p, d_C5pC4p, d_3e_C3pC4p, d_C3pO3p, d_O3pP)
bb_angles_3e = (a_O3pPO5p, a_PO5pC5p, a_O5pC5pC4p, a_3e_C5pC4pC3p, a_3e_C4pC3pO3p, a_C3pO3pP)

