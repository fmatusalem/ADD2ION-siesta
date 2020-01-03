# ADD2ION-siesta
Program to add the self-energy potential to .ion siesta files, used in DFT-1/2 calculations.

This program needs the VTOTAL files calculated with the atomic program and the INP.ae-05 used to calculate the VTOTAL.ae-05.

VTOTAL1.ae Atomic potential for the atom

VTOTAL1.ae-05 Atomic potential for the ionized atom

BY FILIPE MATUSALEM, DEC 2019 filipematus@gmail.com , based on add2ion.f90 by L F Guimaraes

Compile with: $gfortran add2ion_matusa.f90 -o add2ion_matusa

run the program entering the .ion file and CUT value as arguments. The program also accepts a third argument that changes the amplitude.

$add2ion_matusa Si.ion 2.69 #for silicon with CUT=2.69 and amplitude = 1

$add2ion_matusa Si.ion 2.69 -1 #for silicon with CUT=2.69 and amplitude = -1
