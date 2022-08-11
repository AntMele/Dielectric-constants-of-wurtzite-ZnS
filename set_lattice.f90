!
! Copyright (C) 2013 A. Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE set_lattice(at, bg, lattice_type)
USE kinds, ONLY : dp
IMPLICIT NONE
CHARACTER(LEN=*) :: lattice_type 
REAL(DP), INTENT(OUT) :: at(3,3)
REAL(DP), INTENT(OUT) :: bg(3,3)

at=0.0_DP

IF (lattice_type=='SC' .OR. lattice_type=='sc') THEN
   at(1,1)=1.0_DP
   at(2,2)=1.0_DP
   at(3,3)=1.0_DP
ELSEIF (lattice_type=='FCC' .OR. lattice_type=='fcc') THEN
   at(2,1)=0.5_DP
   at(3,1)=0.5_DP
   at(1,2)=0.5_DP
   at(3,2)=0.5_DP
   at(1,3)=0.5_DP
   at(2,3)=0.5_DP
ELSEIF (lattice_type=='BCC' .OR. lattice_type=='bcc') THEN
   at(1,1)=-0.5_DP
   at(2,1)=0.5_DP
   at(3,1)=0.5_DP
   at(1,2)=0.5_DP
   at(2,2)=-0.5_DP
   at(3,2)=0.5_DP
   at(1,3)=0.5_DP
   at(2,3)=0.5_DP
   at(3,3)=-0.5_DP
ELSEIF (lattice_type=='HEX' .OR. lattice_type=='hex') THEN
   at(1,1)=0.5_DP
   at(2,1)=SQRT(3.00_DP)/2.00_DP
   at(1,2)=-0.5_DP
   at(2,2)=SQRT(3.00_DP)/2.00_DP
   at(3,3)=SQRT(8.00_DP/3.00_DP)
ELSE
   WRITE(6,'("Lattice type not programmed")')
   STOP 1
ENDIF

CALL recips (at, bg)

RETURN
END SUBROUTINE set_lattice

