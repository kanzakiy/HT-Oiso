Hydrothermal circulation + O isotope model

Oiso+HT1_v6_irr.f90 --> calculate steady state O isotopes of rock and porewater; including dispersion 
HT1_v8uw_irr.f90    --> Hydrothermal circulation without considering lateral free flow

Need UMFPACK & STEAM libraries
See https://github.com/PetterS/SuiteSparse/tree/master/UMFPACK for UMFPACK
See MEMO.txt for making STEAM library

First run hydrothermal model: 
    gfortran HT1_v8uw_irr.f90 umf4_f77wrapper.o -lumfpack -lamd -lcholmod -lcolamd -lsuitesparseconfig -lopenblas -lsteam -g -fcheck=all
    ./a
Then run isotope model:
    gfortran Oiso+HT1_v6_irr.f90 umf4_f77wrapper.o -lumfpack -lamd -lcholmod -lcolamd -lsuitesparseconfig -lopenblas -g -fcheck=all
    ./a
    