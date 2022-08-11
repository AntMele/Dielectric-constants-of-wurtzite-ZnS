!
! Copyright (C) 2013 A. Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
SUBROUTINE recips (at, bg)
  !---------------------------------------------------------------------
  !
  !   This routine generates the reciprocal lattice vectors b1, b2, b3
  !   from the principal vectors a1, a2, a3 of the Bravais lattice. 
  !   On input the at(:,i) are a1, a2, a3, in units of a. 
  !   On output the bt(:,i) are b1, b2, b3, in units of 2 pi/a.
  !
  !     first the input variables
  !
  USE kinds, ONLY: DP
  IMPLICIT NONE
  REAL(DP) :: at (3,3), bg (3,3)
  ! input: direct lattice vector
  ! output: reciprocal lattice vector
  !
  !   then the local variables
  !
  REAL(DP) :: den
  ! the denominator
  !
  !    first we compute the denominator
  !
  den=at(1,1)*(at(2,2)*at(3,3)-at(3,2)*at(2,3))-   &
      at(1,2)*(at(2,1)*at(3,3)-at(3,1)*at(2,3))+   &
      at(1,3)*(at(2,1)*at(3,2)-at(3,1)*at(2,2))
!
!  and here the three reciprocal lattice vectors
!
  bg(1,1) = (at(2,2)*at(3,3)-at(2,3)*at(3,2))/den
  bg(2,1) = (at(3,2)*at(1,3)-at(3,3)*at(1,2))/den
  bg(3,1) = (at(1,2)*at(2,3)-at(1,3)*at(2,2))/den

  bg(1,2) = (at(2,3)*at(3,1)-at(2,1)*at(3,3))/den
  bg(2,2) = (at(3,3)*at(1,1)-at(3,1)*at(1,3))/den
  bg(3,2) = (at(1,3)*at(2,1)-at(1,1)*at(2,3))/den

  bg(1,3) = (at(2,1)*at(3,2)-at(2,2)*at(3,1))/den
  bg(2,3) = (at(3,1)*at(1,2)-at(3,2)*at(1,1))/den
  bg(3,3) = (at(1,1)*at(2,2)-at(1,2)*at(2,1))/den
!
!  Print the input and output vectors
!
  write(6,'(/,"Direct lattice vectors")')
  write(6,'(3x,"a1=(", 2(f15.5,","),f15.5, " )")') at(:,1)
  write(6,'(3x,"a2=(", 2(f15.5,","),f15.5, " )")') at(:,2) 
  write(6,'(3x,"a3=(", 2(f15.5,","),f15.5, " )")') at(:,3) 

  write(6,'(/,"Reciprocal lattice vectors")')
  write(6,'(3x,"b1=(", 2(f15.5,","),f15.5, " )")') bg(:,1) 
  write(6,'(3x,"b2=(", 2(f15.5,","),f15.5, " )")') bg(:,2) 
  write(6,'(3x,"b3=(", 2(f15.5,","),f15.5, " )")') bg(:,3) 
  write(6,*)

  RETURN
END SUBROUTINE recips
