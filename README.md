# CarbonNanotube-ForsterEnergyTransfer
Calculate exciton energy transfer rate between different types of carbon nanotubes

-This program uses the output files from the program CarbonNanotube-ExcitonEnergy to calculate the energy transfer rate between two types of carbon nanotubes.
-The current version only consider parallel tubes.
-The quasiparticle wavefunctions are 1st nearest neighbor Bloch functions.
-The code assumes that the mesh size in the k-space is equal for all CNTs. This assumption is first used in the soubroutine to find points that bands cross each other.

Structure of the program:
-physicalConstants.f90 : Contains the physical constants that are either universal or considered universal for all types of carbon nanotubes.
-inputParameters.f90 : Contains the simulation parameters that could vary for each carbon nanotube or simulation i.e. cnt chirality, number of unit cells, number of k-points, subband number, threshold energy that is used as a cutoff to determine the limits of k-points.
-dataClass.f90 : Contains the procedures to write the outputs of the simulation and read the data from the output files of ExcitonEnergy program.
-smallFns.f90 : Contains some small miscellaneous functions.

