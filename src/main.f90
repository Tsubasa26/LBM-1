module basic
  implicit none
  integer(4) :: IUT6 = 6
  type lbmmodel
    real(8),allocatable :: w(:)
    real(8),allocatable :: e(:,:)
  endtype lbmmodel

  type field
    integer(4)          :: nx,ny,nq
    real(8)             :: Lx,Ly
    real(8)             :: dx,dy
    real(8),allocatable :: x(:)
    real(8),allocatable :: y(:)
    real(8),allocatable :: rho(:,:)
    real(8),allocatable :: u(:,:)
    real(8),allocatable :: v(:,:)
    real(8),allocatable :: f(:,:,:)
    real(8),allocatable :: feq(:,:,:)
  endtype field

  type solver
    integer :: ntime
    real(8) :: usetime
  endtype solver
contains
  function nfield(nx,ny,nq,Lx,Ly)
    implicit none
    integer, intent(in) :: nx,ny,nq
    real(8), intent(in) :: Lx,Ly
    type(field) :: nfield
    integer :: i,j
    nfield%nx = nx
    nfield%ny = ny
    nfield%nq = nq
    nfield%Lx = Lx
    nfield%Ly = Ly
    nfield%dx = Lx / dble(nx)
    nfield%dy = Ly / dble(ny)
    write(IUT6,'(a18, 3i5)') "allocate nx ny nq:",nfield%nx, nfield%ny, nfield%nq
    allocate(nfield%x(nx))
    allocate(nfield%y(ny))
    allocate(nfield%rho(nx,ny))
    allocate(nfield%u(nx,ny))
    allocate(nfield%v(nx,ny))
    ! LBM
    allocate(nfield%f(nx,ny,nq))
    allocate(nfield%feq(nx,ny,nq))
    ! Set Field
    do i = 1, nfield%nx
      nfield%x(i) = nfield%dx*(i-1)
    enddo
    do j = 1, nfield%ny
      nfield%y(j) = nfield%dy*(j-1)
    enddo
  endfunction nfield

  function d2q9(nd, nq)
    implicit none
    type(lbmmodel)     :: d2q9
    integer,intent(in) :: nd,nq
    ! w : weight
    ! e : velocity
    ! LBM
    allocate(d2q9%w(nq))
    allocate(d2q9%e(nd,nq))
    d2q9%w(1)   = 4.d0/9.d0
    d2q9%w(2:5) = 1.d0/9.d0
    d2q9%w(6:9) = 1.d0/36.d0
    d2q9%e(1,:) = (/0.0d0,1.0d0,0.0d0,-1.0d0,0.0d0,1.0d0,-1.0d0,-1.0d0,1.0d0/)
    d2q9%e(2,:) = (/0.0d0,0.0d0,1.0d0,0.0d0,-1.0d0,1.0d0,1.0d0,-1.0d0,-1.0d0/)
    write(IUT6,'(a19)') "set D2Q9 model Done"
  endfunction d2q9

  function collision()
    implicit none


end module basic

program main
  use basic
  implicit none
  integer :: i, j, t
  type(field)    :: v
  type(lbmmodel) :: model
  type(solver)   :: sol

v = nfield(11, 11, 9, 1.0d0, 1.0d0)
model = d2q9(2, 9)
sol%ntime = 10
do t = 1, sol%ntime
  write(IUT6, '(a10, i6)') "Iteration", t
  ! Collision
  ! Stream
  ! BB
  ! Boundary

enddo

! save data
open(10,file="u.dat")
do i = 1, v%nx
  do j = 1, v%ny
    write(10,*) v%x(i), v%y(j), v%u(i, j)
  enddo
    write(10,*)
enddo
close(10)

end program main
