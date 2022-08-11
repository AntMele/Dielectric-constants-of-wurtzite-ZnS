!
! Copyright (C) 2013 A. Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE input_cb()
USE kinds, ONLY : dp
USE cbmod, ONLY : xk, wk, nks, nbnd, ecut, crystal_name, &
                    nk1, nk2, nk3, nenergy, emin, emax, sigma, leps, pol

IMPLICIT NONE
INTEGER :: nks_dum
REAL(DP), ALLOCATABLE :: xk_dum(:,:)
INTEGER, ALLOCATABLE :: nk_dum(:)
REAL(DP) :: xkmod
REAL(DP) :: delta_k(3)
CHARACTER(LEN=20) :: what

INTEGER :: ik, jk

WRITE(6,'("Crystal name?")')
READ(5,*) crystal_name

WRITE(6,'("Cut-off energy (units Ry (2\pi /a )^2) ?")')
READ(5,*) ecut

WRITE(6,'("Number of bands?")') 
READ(5,*) nbnd

WRITE(6,'("Bands or eps (= immaginary part of dielectric constant) ?")')
READ(5,*) what

IF ( TRIM(what)=='eps' .OR. TRIM(what)=='EPS') THEN
   WRITE(6,'("Number of k points?")')
   READ(5,*) nk1, nk2, nk3

   WRITE(6,'("Minimum and maximum energie for eps? (in eV) ")')
   READ(5,*) emin, emax

   WRITE(6,'("Number of energy points?")')
   READ(5,*) nenergy

   IF (nenergy<1) THEN
      WRITE(6,*) 'Wrong number of energy points'
      STOP 1
   ENDIF

   WRITE(6,'("Sigma?")')
   READ(5,*) sigma
   WRITE(6,'("Polarization? 1(=x), 2(=y) or 3(=z)")')
   READ(5,*) pol

   leps=.TRUE.
ELSE

   WRITE(6,'("Number of k points?")') 
   READ(5,*) nks_dum

   ALLOCATE(xk_dum(3,nks_dum))
   ALLOCATE(nk_dum(nks_dum))

   DO ik=1, nks_dum
      READ(5,*) xk_dum(1,ik), xk_dum(2,ik), xk_dum(3,ik), nk_dum(ik)
   ENDDO

   IF (nks_dum==1) THEN
      nks=nks_dum
   ELSEIF (nks_dum>1) THEN
      nks=SUM(nk_dum(1:nks_dum-1))+1
      DO ik=1,nks_dum-1
         IF (nk_dum(ik)==0) nks=nks+1
      ENDDO
   ELSE
      WRITE(6,*) 'input_cb'
      WRITE(6,*) 'wrong number of k points'
      STOP 1
   ENDIF
   ALLOCATE(xk(3,nks))
   ALLOCATE(wk(nks))

   IF (nks_dum==1) THEN
      xk(:,1) = xk_dum(:,1)
      wk(1) = nk_dum(1)
   ELSE
      nks=1
      xk(:,nks)=xk_dum(:,1)
      wk(nks)= 0.0_DP
      DO ik=1,nks_dum-1
         IF (nk_dum(ik) > 0 ) THEN
            delta_k(:) = (xk_dum(:,ik+1)-xk_dum(:,ik))/nk_dum(ik)
            xkmod=SQRT(delta_k(1)**2+delta_k(2)**2+delta_k(3)**2)
            DO jk=1, nk_dum(ik)
               nks=nks+1
               xk(:,nks) = xk_dum(:,ik) + delta_k(:) * jk
               wk(nks) = wk(nks-1) + xkmod
            ENDDO
         ELSE
            nks=nks+1
            xk(:,nks) = xk_dum(:,ik+1)
            wk(nks) = wk(nks-1)
         ENDIF
      ENDDO
   ENDIF
   leps=.FALSE.

   DEALLOCATE(xk_dum)
   DEALLOCATE(nk_dum)
ENDIF

RETURN
END SUBROUTINE input_cb
