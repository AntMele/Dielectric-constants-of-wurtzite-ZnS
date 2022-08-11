!
! Copyright (C) 2013 A. Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This routine is a cleanup of a similar routine contained 
! (or once contained) in the Quantum ESPRESSO package. See 
! http://www.quantum-espresso.org
!
!
SUBROUTINE diagonalize (npw,nbnd,h,e,v)
!
! f90 interface to LAPACK routine ZHEEVX which calculates
! nbnd eigenvalues and eigenvectors of a complex hermitean matrix h. 
! 
USE kinds, ONLY : dp

IMPLICIT NONE

INTEGER, INTENT(IN) :: npw
INTEGER, INTENT(IN) :: nbnd

COMPLEX(DP), INTENT(INOUT) :: h(npw,npw)   ! matrix to be diagonalized
REAL(DP),    INTENT(OUT)   :: e(nbnd)      ! eigenvalues
COMPLEX(DP), INTENT(OUT)   :: v(npw,nbnd)  ! eigenvectors (column-wise)

INTEGER  :: lwork,  &! auxiliary variable
            info,   &! flag saying if LAPACK execution was OK
            m        ! number of eigenvalues

CHARACTER(LEN=1) :: jobz, range, uplo  ! select the task in LAPACK

REAL(DP), ALLOCATABLE    :: rwork(:)     ! work array used by LAPACK
COMPLEX(DP), ALLOCATABLE :: work(:)      ! as above
INTEGER, ALLOCATABLE     :: iwork(:)     !    "
INTEGER, ALLOCATABLE     :: ifail(:)     !    "
COMPLEX(DP)              :: copt_sz      ! size of the work array
REAL(DP)                 :: rdummy, zero ! dummy variable, zero
REAL(DP), ALLOCATABLE    :: ee(:)        ! axiliary space for eigenvalues
!
!   Initialize flags
!
jobz  = 'V' ! compute eigenvalues and eigenvectors
uplo  = 'U' ! LAPACK routines use the upper triangle of the input matrix
range = 'I' ! compute bands from 1 to nbnd

zero = 0.0_DP
v(:,:) = (0.0_DP,0.0_DP)
!
! allocate arrays of known size
!
ALLOCATE( ee(npw) )
ALLOCATE( rwork(7*npw) )
ALLOCATE( iwork(5*npw) )
ALLOCATE( ifail(npw) )
!
! ZHEEVX calculates selected eigenvalues and eigenvectors
!
!
! calculate optimal size of work arrays and allocate them
!
lwork = -1  ! flag saying that the size of work arrays 
            ! has to be estimated by LAPACK 
rdummy=0.0_DP
CALL zheevx( jobz, range, uplo, npw, h, npw, rdummy, rdummy, 1, nbnd, zero, &
            m, ee, v, npw, copt_sz, lwork, rwork, iwork, ifail, info )
       
IF (ABS(info) /= 0) THEN
   WRITE(6,'("Error finding lwork, info= ", i5)') info
   STOP 1
ENDIF
!
! allocate work
!
lwork  = NINT(REAL(copt_sz))
ALLOCATE(work(lwork))
!
! and diagonalize the matrix
!
CALL zheevx(jobz, range, uplo, npw, h, npw, rdummy, rdummy, 1, nbnd, zero, &
            m, ee, v, npw, work, lwork, rwork, iwork, ifail, info)
!
IF (ABS(info) /= 0) THEN
   WRITE(6,'("Error in the diagonalization, info= ", i5)') info
   STOP 1
ENDIF
!
!
!  NB: h is overwritten by this routine. We save the eigenvalues on e
!      for plotting
!
e(1:nbnd)=ee(1:nbnd) 
       
DEALLOCATE(work)
DEALLOCATE(rwork)
DEALLOCATE(iwork)
DEALLOCATE(ifail)
DEALLOCATE(ee)

RETURN
END SUBROUTINE diagonalize
