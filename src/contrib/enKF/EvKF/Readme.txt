These files illustrate the model used for storing, reading and writing the
ensemble of model states.   It assumes that each member of the ensemble is
stored as a record in a direct unformatted Fortran90 file, in single precision.

Note that if you using the square root scheme with the random 
rotation of anomalies AND the local analysis,
you have to make sure you use the same rotation for each grid point.
You may then compute the random matrix outside the analysis routine, and pass
the same matrix in each local call.


Analysis routines:

analysis.F90:   Routine updated 02/02/2007 including options for perturbed
                observations or symmetric square root schemes

analysis2_EnOI.F90: Ensemble OI version of the perturbed observations analysis

Comments and bug reports to  Geir.Evensen@hydro.com and
Laurent.Bertino@nersc.no
