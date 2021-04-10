# COSMOL2OpenFOAM
Convert COMSOL mesh files (*.mphtxt) to OpenFOAM mesh files. 

Tested on COMSOL 4.1, 5.2, 5.3, 5.3a. For other versions, please add the version info to the first part in mphtxtToFoam.m, or the program will exit on detection of the COMSOL version.

Tested on OpenFOAM 4.1. Should also work on later versions.

The head lines in COMSOL files are different. Three functions are made for v4.1, v5.2, and v.53, respectively. For newer versions of COMSOL, please compare the head lines and modify the readMphtxt53.m to generate the desired function for reading the files.

If you find this program helpful, please cite the paper Wang, Ningyu, Maša Prodanović, and Hugh Daigle. "Nanopaint application for flow assurance with electromagnetic pig." Journal of Petroleum Science and Engineering 180 (2019): 320-329.

