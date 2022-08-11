!
! Copyright (C) 2013 A. Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE deallocate_all
USE cbmod, ONLY : g, et, xk, wk, eps, leps
IMPLICIT NONE

DEALLOCATE(g)
DEALLOCATE(xk)
IF (leps) THEN
   DEALLOCATE(eps)
ELSE
   DEALLOCATE(wk)
ENDIF

RETURN
END SUBROUTINE deallocate_all
