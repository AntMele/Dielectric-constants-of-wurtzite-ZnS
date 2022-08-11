!
! Copyright (C) 2013 A. Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE calculate_radius(ecut, xk, nks, gcutm2)
!
!  This routine calculates the radius of the sphere that contains all
!  the spheres |k+G|^2 < ecut, for all the k points.
!
USE kinds, ONLY : DP
IMPLICIT NONE
INTEGER, INTENT(IN) :: nks
REAL(DP), INTENT(IN) :: ecut, xk(3,nks)
REAL(DP), INTENT(OUT) :: gcutm2 
REAL(DP) :: kmax, kmod
INTEGER :: ik
!
!  Find the maximum modulus of the k vectors
!
kmax=0.0_DP
DO ik=1,nks
   kmod=xk(1,ik)**2 + xk(2,ik)**2 + xk(3,ik)**2
   IF (kmod > kmax) kmax=kmod
ENDDO
!
!  Set the radius of the sphere
!
gcutm2 = (sqrt(ecut) + sqrt (kmax))**2

write(6,'("Kinetic energy cut-off in units Ry / (2\pi/a)**2  ",f15.10)') ecut
write(6,'("Radius of the G vectors sphere ** 2 ((2\pi/a)**2) ",f15.10)') gcutm2
write(6,'("Maximum k (2\pi/a)",f15.10)') sqrt(kmax)

RETURN
END SUBROUTINE calculate_radius
