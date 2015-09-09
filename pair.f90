subroutine initPAIR
! this subroutine calculate non-Hamaker interactions from pair potentials

use ellipsoid
use transform
use system
use spline
implicit none
real*8 radius ! radius of the particles, must be all equal and particles should be spherical
LOGICAL :: file_exists
integer IOS, ii, i,j
real*8 dist(100), energy(100)
real*8 yp1, ypn

! check radius
radius = Aell(1,1)

do j = 1, NNN
 do i = 1,3
  if(Aell(i,j).ne.radius) then
     print*, 'initPAIR: all particles must be spherical and have same size' 
     stop
  endif
 enddo
enddo

! check if pair.txt is present or not
INQUIRE(FILE="pair.txt", EXIST=file_exists)   ! file_exists will be TRUE if the file

if(file_exists.eqv..false.) then
 print*, 'initPAIR: pair.txt does not exists, pairwise additive non-Hamaker energy will not be calculated'
 return
endif

! load pairwise potential
ii = 1
open(file='pair.txt', unit=10)
do while(.true.)
read(10,*,IOSTAT=IOS)dist(ii),energy(ii)
if(IOS.ne.0)exit
ii=ii+1
enddo
close(10)

! make spline
allocate(x(ii))
allocate(y(ii))
allocate(y2(ii))

do i = 1, ii
x(i)=dist(i)
y(i)=energy(i)
enddo


yp1 = 1E31
ypn = 1E31
call splines(x,y,ii,yp1,ypn,y2)
nspline = ii
end

subroutine calcPAIR
! this subroutine calculate non-Hamaker interactions from pair potentials

use ellipsoid
use transform
use system
use spline
implicit none
integer, parameter :: Nlat = 5 ! number of lattices to scan around central cell
real*8 radius ! radius of the particles, must be all equal and particles should be spherical
integer j,i, ix,iy,iz
real*8 dist
real*8 energy, energypair
real*8 vect1T(3), vect1R(3),vect2T(3),vect2R(3),vectdiff(3)

call initPAIR ! init pair calculation

! loop over cell indexes 
energy = 0.0
do j = 1, NNN ! loop over all particles in central cell

      vect1T(1) = Rellf(1,j)*delta*dfloat(dimx)
      vect1T(2) = Rellf(2,j)*delta*dfloat(dimy)
      vect1T(3) = Rellf(3,j)*delta*dfloat(dimz)
      vect1R = MATMUL(IMAT,vect1T) ! coordinates of particle in real space

 do ix = -Nlat,Nlat 
 do iy = -Nlat,Nlat 
 do iz = -Nlat,Nlat 
   do i = 1, NNN ! loop over all particles in cell ix,iy,iz

      vect2T(1) = Rellf(1,i)*delta*dfloat(dimx)+dfloat(dimx*ix)*delta
      vect2T(2) = Rellf(2,i)*delta*dfloat(dimy)+dfloat(dimy*iy)*delta
      vect2T(3) = Rellf(3,i)*delta*dfloat(dimz)+dfloat(dimz*iz)*delta
      vect2R = MATMUL(IMAT,vect2T) ! coordinates of particle in real space

      vectdiff = vect1R-vect2R
      dist = norm2(vectdiff) 

      if(dist.ne.0.0) then ! not the same particle

         if(dist.lt.x(1)) then ! data not in spline
           print*, 'calcPAIR: particles closer tahn minimal distance in pair.txt!'
           stop
         endif
         if(dist.lt.x(nspline)) then
           call splint(x,y,y2,nspline,dist,energypair) 
           energy=energy+energypair
         endif

    endif ! dist
    enddo ! i
  enddo ! ix
  enddo ! iy
  enddo ! iz
enddo ! j

print*,'calcPAIR: energy is', energy

open(file='F_pair.dat', unit=10)
write(10,*)energy
close(10)

end



