module system 
real*8 delta
real*8 cdiva
integer  dimx 
integer  dimy 
integer  dimz 
integer PBC(6)
endmodule

module ellipsoid
integer kaptype
integer NNN
real*8, allocatable :: rotmatrix(:,:,:)
real*8, allocatable :: Aell(:,:)
real*8, allocatable :: AellS(:,:)
real*8, allocatable :: AellL(:,:)
real*8, allocatable :: AellX(:,:)
real*8, allocatable :: AAA(:,:,:)
real*8, allocatable :: AAAS(:,:,:)
real*8, allocatable :: AAAL(:,:,:)
real*8, allocatable :: AAAX(:,:,:)
real*8, allocatable :: Rell(:,:)
real*8, allocatable :: Rellf(:,:)
real*8, allocatable :: orient(:,:)
real*8, allocatable :: echarge(:)
real*8, allocatable :: sigma(:)
real*8, allocatable :: eeps(:)
end module

module transform
real*8 gama0
real*8 MAT(3,3)
real*8 TMAT(3,3)
real*8 IMAT(3,3)
endmodule

