
subroutine readinput
use ellipsoid
use transform
use system
implicit none

! Input related variables
character (len=100)  buffer,label
integer pos
integer, parameter :: fh = 15
integer ios
integer line, linemax
integer i, j
character(len=50) :: filename = 'DEFINITIONS.txt'
character basura
integer ndi
real*8 ndr

! not defined variables, change if any variable can take the value

PBC = 1

ndi = -1e5
ndr = -1.0d10

verbose = 5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Check validity of input
!

dimx = ndi
dimy = ndi
dimz = ndi

delta = ndr
cdiva = ndr
gama0 = ndr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Control file variables

line = 0
ios = 0

open(fh, file=filename)

if(rank.eq.0)print*, 'parser:', 'Reading parameters from ', filename

! ios is negative  if an end of record condition is encountered or if
! an endfile condition was detected.  It is positive  if an error was
! detected.  ios is zero otherwise.

do while (ios == 0)

 read(fh, '(A)', iostat=ios) buffer
 if (ios == 0) then
 line = line + 1

! Find the first instance of whitespace.  Split label and data.

 pos = scan(buffer, ' ')

 label = buffer(1:pos)
 buffer = buffer(pos+1:)


 select case (label)

 case ('PBC')
   read(buffer, *, iostat=ios) PBC(1),PBC(2),PBC(3),PBC(4),PBC(5),PBC(6)
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

   do j = 1,5,2
    if((PBC(j).eq.1).and.(PBC(j+1).ne.1)) then 
      print*, 'parser:', 'Error in PBC'
      stop
    endif
    if((PBC(j+1).eq.1).and.(PBC(j).ne.1)) then
      print*, 'parser:', 'Error in PBC'
      stop
    endif
   enddo

 case ('dimx')
   read(buffer, *, iostat=ios) dimx
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('delta')
   read(buffer, *, iostat=ios) delta
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('cdiva')
   read(buffer, *, iostat=ios) cdiva
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('dimy')
   read(buffer, *, iostat=ios) dimy
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('dimz')
   read(buffer, *, iostat=ios) dimz
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

 case ('gama')
   read(buffer, *, iostat=ios) gama0
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

   do i = 1, nst
   read(fh,*)sts(i)
   enddo 

 case ('kaptype')
   read(buffer, *, iostat=ios) kaptype
   if(rank.eq.0)print*, 'parser:','Set ',trim(label),' = ',trim(buffer)

   select case (kaptype)
    case(1) 
     read(fh, *) basura
     read(fh, *)NNN

     if(NNN.ne.0) then

     call allocateell
     read(fh, *), basura
     do j = 1, NNN
     read(fh, *), Rellf(1,j), Rellf(2,j), Rellf(3,j)
     if(rank.eq.0)print*, 'parser:','Set particle',j,'pos to',  Rellf(1,j), Rellf(2,j), Rellf(3,j)
     enddo
     read(fh, *), basura
     do j = 1, NNN
     read(fh, *), Aell(1,j), Aell(2,j), Aell(3,j)
     if(rank.eq.0)print*, 'parser:','Set particle',j,'axis to',  Aell(1,j), Aell(2,j), Aell(3,j)
     enddo
     read(fh, *), basura
     do j = 1, NNN
     read(fh, *), rotmatrix(1,1,j), rotmatrix(1,2,j), rotmatrix(1,3,j)
     read(fh, *), rotmatrix(2,1,j), rotmatrix(2,2,j), rotmatrix(2,3,j)
     read(fh, *), rotmatrix(3,1,j), rotmatrix(3,2,j), rotmatrix(3,3,j)
     if(rank.eq.0) then
         print*, 'parser:','Set particle',j,'rotation to:'
         print*, 'parser:', rotmatrix(1,1,j), rotmatrix(1,2,j), rotmatrix(1,3,j)
         print*, 'parser:', rotmatrix(2,1,j), rotmatrix(2,2,j), rotmatrix(2,3,j)
         print*, 'parser:', rotmatrix(3,1,j), rotmatrix(3,2,j), rotmatrix(3,3,j)
     endif
     enddo

     read(fh, *), basura
     do j = 1, NNN
     read(fh, *), sigma(j)
     if(rank.eq.0)print*, 'parser:','Set particle',j,'surface coverage to', sigma(j)
     enddo

     read(fh, *), basura
     do j = 1, NNN
     read(fh, *), echarge(j)
     if(rank.eq.0)print*, 'parser:','Set particle',j,'charge to', echarge(j)
     enddo
     read(fh, *), basura
     do j = 1, NNN
     read(fh, *), eeps(j)
     if(rank.eq.0)print*, 'parser:','Set particle',j,'hydrophobicity to', eeps(j)
     enddo

     endif ! NNN
endselect
endselect

endif

enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Check validity of input
! 


if(dimx.eq.ndi)call stopundef('dimx')
if(dimy.eq.ndi)call stopundef('dimy')
if(dimz.eq.ndi)call stopundef('dimz')

if(delta.eq.ndr)call stopundef('delta')
if(cdiva.eq.ndr)call stopundef('cdiva')
if(gama0.eq.ndr)call stopundef('gama')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine

subroutine stopundef(namevar)
character(len=*) :: namevar
print*, 'parser:', 'Variable ', namevar, ' is undefined '
end
