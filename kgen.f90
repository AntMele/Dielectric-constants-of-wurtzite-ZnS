!
! Copyright (C) 2013 A. Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE kgen(n1, n2, n3, bg)
USE kinds, ONLY : DP
USE cbmod, ONLY : nks, xk
IMPLICIT NONE
INTEGER :: n1, n2, n3    ! the mesh of k points
REAL(DP), INTENT(IN) :: bg(3,3)
INTEGER :: i, j, k, ik

nks = n1 * n2 * n3

ALLOCATE(xk(3,nks))

ik=0
DO i=1, n1
   DO j=1, n2
      DO k=1, n3
         ik = ik + 1
         xk(:,ik) = (i-1) * bg(:,1) / n1 + (j-1) * bg(:,2) / n2 +  &
                    (k-1) * bg(:,3) / n3
      ENDDO
   ENDDO
ENDDO  

RETURN
END SUBROUTINE kgen

