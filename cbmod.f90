!
! Copyright (C) 2013 A. Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE cbmod
USE kinds, ONLY : DP

IMPLICIT NONE
SAVE
PUBLIC 

CHARACTER(LEN=10) :: crystal_name ! The name of the crystal

REAL(DP) :: alat      ! Size of the unit cell
REAL(DP) :: ecut      ! Kinetic energy cut-off
REAL(DP) :: gcutm2    ! Square of the radius of the sphere in G space 
REAL(DP) :: tpiba2    ! Two pi divided a_0

REAL(DP) :: at(3,3)   ! Direct lattice vectors
REAL(DP) :: bg(3,3)   ! Reciprocal lattice vectors

INTEGER :: ngm                   ! How many G vector we have
REAL(DP), ALLOCATABLE :: g(:,:)  ! Reciprocal lattice vectors
INTEGER,  ALLOCATABLE :: igk(:)  ! index of G vector that we can consider
REAL(DP) :: cb_parameters(12,3)      ! CB parameters (column 1 for |G|^2, column 2 and 3 for sym. and antisym. form factors)
INTEGER  :: npw                      ! The dimension of the Hamiltonian
INTEGER  :: nbnd                     ! The number of calculated bands
COMPLEX(DP), ALLOCATABLE :: h(:,:)   ! The Hamiltonian matrix
COMPLEX(DP), ALLOCATABLE :: evc(:,:) ! The wavefunctions
REAL(DP), ALLOCATABLE :: et(:)       ! The eigenvalues of the Hamiltonian

INTEGER :: nks                   ! The number of k points
REAL(DP), ALLOCATABLE :: xk(:,:) ! The coordinates of the k points
REAL(DP), ALLOCATABLE :: wk(:)   ! The xcoordinate of the k point for output

LOGICAL :: leps                  ! true if computes permittivity
INTEGER :: nk1, nk2, nk3         ! the mesh of k points
INTEGER :: nenergy               ! number of energy points
INTEGER :: pol                   ! the polarization ( it could be 1,2 or 3 (where 1=x, 2=y and 3=z) )

REAL(DP) :: sigma, emin, emax    ! a real number
REAL(DP), ALLOCATABLE :: eps(:)  ! The immaginary part of dielectric constant


END MODULE cbmod
