# Description
BarycentricTurbulence represents an OpenFOAM utility for calculating barycentric coordinates of any volSymmTensorField representing Reynolds stress tensor. Utlity can
operate on any number of fields, whoose names are provided as a word list through the option -fields.
# Compilation
Compilation is done through the wmake wrapper. Tested on version v2012.
# Examples
BarycentricTurbulence -latestTime -fields "(R)"
BarycentricTurbulence -latestTime -fields "(UPrime2Mean RMean)"
