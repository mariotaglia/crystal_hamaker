subroutine allocateellCO

use system 
use ellipsoid

! cuboctahedron
ALLOCATE (Rellf(3,NNN))
ALLOCATE (Rell(3,NNN))
ALLOCATE (Loctall(NNN))
ALLOCATE (Lcubell(NNN))
ALLOCATE (echarge(NNN))
ALLOCATE (sigma(NNN))
ALLOCATE (eeps(NNN))

end subroutine
