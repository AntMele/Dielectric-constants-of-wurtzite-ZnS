# Dielectric-constants-of-wurtzite-ZnS
Here we discuss a program that can compute
(depending by the input) or the imaginary part of dielectric 
constant or the band structure of ZnS with Wurtzite structure 
by the Cohen and Bergstresser method.
Main Reference: Phys. Rev. 164, 1069 (1967).

We discuss the following subroutines:

cb_permittivity.f90        ! The main driver program
cbmod.f90                  ! contains the data structures of the program  
input_cb.f90               ! the input routine
set_cb_parameters.f90 	   ! set the parameters of the pseudopotentials
set_lattice.f90            ! calculate the direct and reciprocal lattice vectors
kgen.f90                   ! generate a mesh of k points used for the sum that appears in 
                           ! the formula for imaginary part of dielectric constant
calculate_radius.f90       ! calculates the radius of the G sphere that
                           ! contains all the spheres |k+G|^2 < Ecut 
ggen.f90                   ! calculates the list of G vectors 
set_hamiltonian.f90        ! set the Hamiltonian in the plane waves basis at each k
diagonalize.f90            ! diagonalize the Hamiltonian using lapack routines 
gaussian.f90               ! return the value of the gaussian used to represent the 
                           ! delta function. 
deallocate_hamiltonian.f90 ! deallocation of the main Hamiltonian arrays
deallocate_all.f90         ! deallocate all the rest
                  
Then we need a few auxiliary routines copied from QE:

recips.f90                 ! To calculate the reciprocal lattice vectors
kind.f90                   ! The kind of the real numbers
constants.f90              ! A few constants (I added a constant to convert eV into Ry)

zheevx.html : the man page of the zheevx routine used to diagonalize the
Hamiltonian.

To compile it you need a fortran90 compiler, such as gfortran, and
compiled blas and lapack libraries.

Modify the file 'src/compile' to link with the blas and lapack of your 
installation or to use a compiler different from gfortran.

If your machine has not the lapack and blas libraries installed,
you can build them by entering the src directory and typing

tar -xzvf lapack-3.2.tar.gz

entering in the lapack-3.2 directory and giving the commands:

cp make.inc.example make.inc
make blaslib
make lapacklib

To compile enter in src and type . compile.
To run the code enter the inputs directory and type
../src/cb.x < input_? 

There are several examples of input:

input_epsx     shows an example of input to obtain the imaginary part 
               of dielectric constant with x-polarization.

input_epsz     shows an example of input to obtain the imaginary part 
               of dielectric constant with z-polarization.

TEST/test.sh   make a simple script to check the convergence of the imaginary part 
               of dielectric constant as a function of k points and sigma for three
               values of the energy. 

input_bands_ZnS    show an input to obtain the band structure of ZnS wurtzite.

There are several plots:

plot_epsx   plot to obtain the imaginary part of dielectric constant 
            with x-polarization.

plot_epsz   plot to obtain the imaginary part of dielectric constant 
            with z-polarization.

plot_bands_ZnS   plot to obtain the band structure of ZnS wurtzite.

In the folder plots_report, you can find the outputs and the plots 
that I put in the report. 
