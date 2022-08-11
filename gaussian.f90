!
! Copyright (C) 2013 A. Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
FUNCTION gaussian(x,sigma)
USE kinds, ONLY : DP
IMPLICIT NONE
REAL(DP) :: gaussian
REAL(DP), INTENT(IN) :: sigma, x
REAL(DP) :: arg

arg=x*x/(2.0_DP*sigma*sigma)
IF (arg < 150_DP) THEN
   gaussian=0.39894228_DP/sigma*exp(-arg)
ELSE
   gaussian=0.0_DP
ENDIF

RETURN
END FUNCTION gaussian
