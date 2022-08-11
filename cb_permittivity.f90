PROGRAM cb_permittivity
!
!  Depending on the input, this program computes or the band structure or the immaginary part of dielectric constant.
!  In the first case, this program reads from input the name of a crystal (ZnSwu, ..), a path of k points in the
!  hexagonal Brillouin zone and some parameters, and writes on output the Cohen-Bergstresser bands of the crystal
!  along that path. Output is in Ry.
!  In the second case, it reads from input some parameters such as the energy interval in eV and writes on output the immaginary
!  part of dielectric constant.
!
!  Reference: Phys. Rev. 164, 1069 (1967).
!
USE kinds,   ONLY : dp
USE cbmod, ONLY : crystal_name, at, bg, ecut, gcutm2,  &
                    nks, xk, wk, ngm, nbnd, leps, emin, emax, nenergy, &
                    sigma, nk1, nk2, nk3, eps, npw, h, evc, et, igk, g, pol
USE constants, ONLY : rydb

IMPLICIT NONE
INTEGER :: ik, iener, ivalbnd, icondbnd, maxvalbnd, ios, ig
REAL(DP) :: ener, deltae, pvc2
REAL(DP), EXTERNAL :: gaussian
COMPLEX(DP) :: pvc
!
! Read the input
!
CALL input_cb()
!
! set the cb parameters and the size of the unit cell
!
CALL set_cb_parameters(crystal_name)
!
! set the direct and reciprocal lattice vectors of the lattice
!
CALL set_lattice(at, bg, 'hex')
!
! generate the mesh of k point for permittivity
!
IF (leps) CALL kgen(nk1, nk2, nk3, bg)
!
! given ecut and the k-points, calculate the radius of the sphere in
! reciprocal space that contains all the necessary G vectors.
!
CALL calculate_radius(ecut, xk, nks, gcutm2)
!
! compute all the reciprocal lattice vectors within a sphere of radius 
! sqrt(gcutm2)
!
CALL ggen(gcutm2)
!
! open the output file
!
IF (leps) THEN
   maxvalbnd=8 ! max number of valence bands (In the case of ZnS wurtzite we have 8 valence bands).
   ALLOCATE(eps(nenergy)) ! allocate space for immaginary part of dielectric constant
   eps=0.0_DP
   IF (nenergy > 1) THEN 
      deltae=(emax-emin)/(nenergy-1)
   ELSE
      deltae=0.0_DP
      nenergy=1
   ENDIF
   deltae=deltae/rydb ! eV -> Rydberg
   emin=emin/rydb ! eV -> Rydberg
   sigma=sigma/rydb ! eV -> Rydberg
ELSE
!
! open the output file
!
   OPEN(unit=26,file='output',status='unknown',err=100,iostat=ios)
   100 IF (ios /= 0) STOP 'opening output'
ENDIF

!
! Cycling over all the k points
!
DO ik=1, nks

!
!   set the hamiltonian matrix
!
   CALL set_hamiltonian(xk(:,ik), ecut)
!
!   and diagonalize it. nbnd bands are calculated.
!
   CALL diagonalize(npw, nbnd, h, et, evc)

   IF (leps) THEN
!
!  in this case accumulate the immaginary part of permittivity due to this k point.
!
      DO iener=1,nenergy
         ener=(emin+deltae*(iener-1)) ! in Rydberg units
         DO ivalbnd=1, maxvalbnd ! cycle over valence bands
            DO icondbnd=maxvalbnd+1,nbnd ! cycle over conduction bands
               pvc=(0.0_DP,0.0_DP)
               DO ig=1, npw  ! cicle over G vector in order to compute p_vc
                  pvc=pvc + CONJG(evc(ig,ivalbnd))*evc(ig,icondbnd)*( xk(pol,ik) + g(pol,igk(ig)) ) ! pvc formula as in the report
               ENDDO
               pvc2=ABS(pvc)**2
               eps(iener)=eps(iener)+ ((1.0_DP/ener)**2)*gaussian(et(icondbnd)-et(ivalbnd)-ener,sigma)*pvc2 ! as in the formula for epsilon.
            ENDDO
         ENDDO
      ENDDO
   ELSE
!
! or write on output the eigenvalues
!
      WRITE(26,'(50f10.4)') wk(ik), et(1:nbnd)
   ENDIF
!
!  deallocate the Hamiltonian variables
!
   CALL deallocate_hamiltonian()

   WRITE(6,*) (100.0*ik)/nks, ' /100'
ENDDO

!
! write on output the immaginary part of dielectric constant
!
IF (leps) THEN
   OPEN(unit=26,file='output',status='unknown',err=200,iostat=ios)
   200 IF (ios /= 0) STOP 'opening output'
   eps=eps/nks
   DO iener=1,nenergy
      ener= emin*rydb + deltae*rydb*(iener-1)  ! Rydberg -> eV
      WRITE(26,'(50f10.4)') ener, eps(iener) ! immaginary part of dielectric constant in arbitrary units
   ENDDO
ENDIF

CALL deallocate_all()
CLOSE(26)
END PROGRAM cb_permittivity
