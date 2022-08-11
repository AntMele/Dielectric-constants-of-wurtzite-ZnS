!
! Copyright (C) 2013 A. Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE set_cb_parameters (crystal_name)
!
! This routine sets the form factors and the size of the conventional 
! cell of the crystal taken by the reference Phys. Rev. 164, 1069 (1967)
! and written here in the appropriate units.
!
USE kinds, ONLY : dp
USE constants, ONLY : aa_to_au, pi
USE cbmod, ONLY : cb_parameters, tpiba2, alat
IMPLICIT NONE
CHARACTER (LEN=*), INTENT(IN)  :: crystal_name
REAL(DP) :: tpiba

cb_parameters=0.0_DP
SELECT CASE (crystal_name)
CASE ('ZnSwu') ! ZnS wurtzite
   cb_parameters(1,1) = (2.00_DP + 2.00_DP/3.00_DP)/2.00_DP ! |G|^2 in units of 2pi/a_wurtzite in the first column
   cb_parameters(1,2) = -0.24_DP !  The symmetric form factors in the second column
   cb_parameters(1,3) = 0.00_DP ! The antisymmetric form factors in the third column
   cb_parameters(2,1) = (3.00_DP)/2.00_DP 
   cb_parameters(2,2) = -0.22_DP
   cb_parameters(2,3) = 0.23_DP 
   cb_parameters(3,1) = (3.00_DP + 5.00_DP/12.00_DP)/2.00_DP
   cb_parameters(3,2) = -0.19_DP
   cb_parameters(3,3) = 0.19_DP
   cb_parameters(4,1) = (5.00_DP + 2.00_DP/3.00_DP)/2.00_DP
   cb_parameters(4,2) = -0.06_DP
   cb_parameters(4,3) = 0.10_DP
   cb_parameters(5,1) = (8.00_DP)/2.00_DP
   cb_parameters(5,2) = 0.03_DP
   cb_parameters(5,3) = 0.00_DP
   cb_parameters(6,1) = (9.00_DP + 5.00_DP/12.00_DP)/2.00_DP
   cb_parameters(6,2) = 0.06_DP
   cb_parameters(6,3) = 0.03_DP
   cb_parameters(7,1) = (10.00_DP + 2.00_DP/3.00_DP)/2.00_DP
   cb_parameters(7,2) = 0.07_DP
   cb_parameters(7,3) = 0.00_DP
   cb_parameters(8,1) = (11.00_DP)/2.00_DP
   cb_parameters(8,2) = 0.07_DP
   cb_parameters(8,3) = 0.02_DP
   cb_parameters(9,1) = (11.00_DP + 5.00_DP/12.00_DP)/2.00_DP
   cb_parameters(9,2) = 0.07_DP
   cb_parameters(9,3) = 0.02_DP
   cb_parameters(10,1) = (12.00_DP)/2.00_DP
   cb_parameters(10,2) = 0.00_DP
   cb_parameters(10,3) = 0.02_DP
   cb_parameters(11,1) = (13.00_DP + 2.00_DP/3.00_DP)/2.00_DP
   cb_parameters(11,2) = 0.04_DP
   cb_parameters(11,3) = 0.01_DP
   cb_parameters(12,1) = (14.00_DP + 2.00_DP/3.00_DP)/2.00_DP
   cb_parameters(12,2) = 0.00_DP
   cb_parameters(12,3) = 0.01_DP
   alat = 3.811_DP * aa_to_au
CASE DEFAULT
   WRITE(6,*) 'set_cb_parameters'
   WRITE(6,*) 'crystal name not known'
   STOP 1
END SELECT

tpiba  = 2.0_DP * pi / alat
tpiba2 = tpiba**2

RETURN
END SUBROUTINE set_cb_parameters
