!
! Copyright (C) 2013 A. Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE ggen(gcutm2)
!
! This routine generates all the reciprocal lattice vectors contained in
! a sphere of radius sqrt(gcutm2) and saves them in the array g.
!
USE kinds, ONLY : dp
USE cbmod, ONLY : at, bg, g, ngm

IMPLICIT NONE
REAL(DP), INTENT(IN)   :: gcutm2

REAL(DP)               :: t(3), tt
INTEGER                :: i, j, k, nx1, nx2, nx3
!
!  Determine the size of a mesh that contains the sphere
!
nx1=int(sqrt(gcutm2)*sqrt(at(1,1)**2+at(2,1)**2+at(3,1)**2))+1
nx2=int(sqrt(gcutm2)*sqrt(at(1,2)**2+at(2,2)**2+at(3,2)**2))+1
nx3=int(sqrt(gcutm2)*sqrt(at(1,3)**2+at(2,3)**2+at(3,3)**2))+1
!
! count the total number of g-vectors such that |g|^2 <= gcutm2
!
ngm = 0
DO i = -nx1, nx1
   DO j = -nx2, nx2
      DO k = -nx3, nx3
         t(:) = i*bg(:,1) + j*bg(:,2) + k*bg(:,3)
         tt   = t(1)*t(1) + t(2)*t(2) + t(3)*t(3)
         IF ( tt <= gcutm2 ) ngm = ngm + 1
      END DO
   END DO
END DO
!
! allocate memory for g-vectors
!
ALLOCATE(g(3,ngm))
!
! re-calculate the g-vectors and their squared norms and store them
!
ngm = 0
DO i = -nx1, nx1
   DO j = -nx2, nx2
      DO k = -nx3, nx3
         t(:) = i*bg(:,1) + j*bg(:,2) + k*bg(:,3)
         tt   = t(1)*t(1) + t(2)*t(2) + t(3)*t(3)
         IF ( tt <= gcutm2 ) THEN
            ngm      = ngm + 1
            g(:,ngm) = t(:)
         END IF
      END DO
   END DO
END DO

RETURN
END SUBROUTINE ggen
