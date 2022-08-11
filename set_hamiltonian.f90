!
! Copyright (C) 2013 A. Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE set_hamiltonian(xk,ecut)
!
!  This routine sets the CB Hamiltonian. 
!
USE kinds, ONLY   : dp
USE constants, ONLY : pi
USE cbmod, ONLY   : g, h, evc, et, ngm, npw, nbnd, tpiba2, &
                    cb_parameters, igk

IMPLICIT NONE

REAL(DP), INTENT(IN) :: xk(3), ecut

INTEGER      :: i, j, ig, indx
REAL(DP)     :: dg(3), g_dot_tau, ssg, ssn, t(3), t2(3), tt, dg2, eps
REAL(DP), ALLOCATABLE :: xkpg2(:)

ALLOCATE(igk(ngm))
ALLOCATE(xkpg2(ngm))


npw=0
DO ig=1,ngm
   t(:) = xk(:) + g(:,ig) 
   tt   = t(1) * t(1) + t(2) * t(2) + t(3) * t(3)
   IF (tt < ecut) THEN
      npw = npw + 1
      xkpg2 (npw) = tt * tpiba2
      igk (npw) = ig
   END IF
END DO

ALLOCATE( h(npw,npw) )
ALLOCATE( evc(npw,nbnd) ) 
ALLOCATE( et(nbnd) )

h(:,:)  = (0.0_DP,0.0_DP)
!
! Only two atomic positions are needed. The origin is in the center of the
! bond.
!
t(1) = 0.5_DP
t(2) = 0.5_DP/SQRT(3.0_DP)
t(3) =-(7.0_DP/16.0_DP)*SQRT(8.0_DP/3.0_DP)
t2(1) = 0.5_DP
t2(2) = 0.5_DP/SQRT(3.0_DP)
t2(3) = -(1.0_DP/16.0_DP)*SQRT(8.0_DP/3.0_DP)
!
!  Loop over G:
!
DO i=1, npw
!
!  The diagonal term: kinetic energy.
!
   h(i,i) = CMPLX(xkpg2(i),0.0_DP)
!
!  There is a double loop, on G and G', but we can use the hermiticity of H.
!
   DO j=i+1,npw
!
! We have non vanishing terms only if |G-G'|^2 is not higher 
! than (14+2/3)/2 (in units (2\pi/a_wurtzite)^2).
      dg(:) = g(:,igk(i)) - g(:,igk(j))
      dg2  = dg(1)**2 + dg(2)**2 + dg(3)**2 
      IF ( dg2 <= 7.0_DP + 1.0_DP/3.0_DP + 0.1_DP) THEN
!
!  tau or tau2 are in units of a=alat, dg is in units (2\pi/a), 
!  their scalar product is adimensional, but we must multiply by 2\pi before
!  calculating the cosine and sine.
!
         g_dot_tau = DOT_PRODUCT(dg(:),t(:)) * pi * 2.0_DP
         ssg       = 0.5_DP*COS(g_dot_tau)
         ssn       = -0.5_DP*SIN(g_dot_tau)
         g_dot_tau = DOT_PRODUCT(dg(:),t2(:)) * pi * 2.0_DP
         ssg       = ssg + 0.5_DP*COS(g_dot_tau)
         ssn       = ssn + 0.5_DP*SIN(g_dot_tau) ! as in the formula that I wrote in the report


!    Here we find the position of |G-G'|^2 in the first column of cb_parameters.

         DO indx=1,12
            IF ( ABS( dg2 - cb_parameters(indx,1) ) < 0.1_DP ) THEN  ! ABS( dg2 - cb_parameters(indx,1) ) > 0.1_DP always, apart when dg2=cb_parameters(indx,1). If you need, you can put a lower precision to generalize.
                h(i,j)   = CMPLX(ssg*cb_parameters(indx,2),ssn*cb_parameters(indx,3)) ! set the hamiltonian using the form factors
                EXIT
            END IF 
         ENDDO

!  
!  H is an hermitean operator, but the diagonalization routine does not
!  need the lower part
!
!         h(j,i)   = CONJG(h(i,j))
      END IF
   END DO
END DO

DEALLOCATE(xkpg2)


RETURN

END SUBROUTINE set_hamiltonian
