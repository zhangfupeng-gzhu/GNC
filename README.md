# GNC (A Monte Carlo code for dynamics of Galactic Nuclear star Cluster with a central massive black hole)

  GNC is a Monte-Carlo numerical method to calculate the dynamical evolution of a cluster of particles around a massive black hole, written by Zhang fupeng, in collabration with Amaro Seoane Pau. Here we present the basic version of the method, incorporates two-body relaxation and the effects of the loss cone. GNC can calculate the two-body relaxation of multiple mass components based on a two-dimensional (energy and angular momentum) Fokker-Planck Monte Carlo method, the steady-state solutions for multiple mass components under given boundary conditions can be obtained. 
  
  The Monte-Carlo approach used by GNC is base on Shapiro, S. L., & Marchant, A. B. 1978, ApJ, 225, 603. The method is also an overhaul of the previous version of the code in Zhang, F., Shao, L., & Zhu, W. 2019, ApJ, 877, 87, and Zhang, F., Chen, X., Shao, L., et al. 2021, ApJ, 923, 139.
 
  The code is written in Fortran 2003, and has been tested using gfortran 11.3.0 and open mpi 4.1.2. We haven't tested it using other versions of fortran compiler or open mpi, probabily work, however.
 
  The code and the related documentary will be released soon. 
