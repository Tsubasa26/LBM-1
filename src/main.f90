module basic
  implicit none
  integer(4) :: IUT6 = 6
  type lbmmodel
    real(8),allocatable :: w(:)
    real(8),allocatable :: e(:,:)
    ! list BB
    integer,allocatable :: lbb(:)
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
    integer(8),allocatable :: ob(:,:)
  endtype field

  type solver
    integer :: ntime
    real(8) :: usetime
    real(8) :: Re,nu,tau
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
    nfield%dx = Lx / dble(nx-1)
    nfield%dy = Ly / dble(ny-1)
    write(IUT6,'(a18, 3i5)') "allocate nx ny nq:",nfield%nx, nfield%ny, nfield%nq
    allocate(nfield%x(nx))
    allocate(nfield%y(ny))
    allocate(nfield%rho(nx,ny))
    allocate(nfield%u(nx,ny))
    allocate(nfield%v(nx,ny))
    ! LBM
    allocate(nfield%f(nx,ny,nq))
    allocate(nfield%feq(nx,ny,nq))
    allocate(nfield%ob(nx,ny))
    ! Set Field
    do i = 1, nfield%nx
      nfield%x(i) = nfield%dx*(i-1)
    enddo
    do j = 1, nfield%ny
      nfield%y(j) = nfield%dy*(j-1)
    enddo
    ! Init
    nfield%u(:,:)     = 0.0d0
    nfield%v(:,:)     = 0.0d0
    nfield%f(:,:,:)   = 0.0d0
    nfield%feq(:,:,:) = 0.0d0
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
    allocate(d2q9%lbb(nq))
    d2q9%w(1)   = 4.d0/9.d0
    d2q9%w(2:5) = 1.d0/9.d0
    d2q9%w(6:9) = 1.d0/36.d0
    d2q9%e(1,:) = (/0.0d0,1.0d0,0.0d0,-1.0d0,0.0d0,1.0d0,-1.0d0,-1.0d0,1.0d0/)
    d2q9%e(2,:) = (/0.0d0,0.0d0,1.0d0,0.0d0,-1.0d0,1.0d0,1.0d0,-1.0d0,-1.0d0/)
    d2q9%lbb = (/1,4,5,2,3,8,9,6,7/)
    write(IUT6,'(a19)') "set D2Q9 model Done"
  endfunction d2q9

  subroutine calcfeq(v, model)
    implicit none
    type(field)    :: v
    type(lbmmodel) :: model
    integer        :: i,j,l
    real(8)        :: eu, uu, s

    v%feq = 0.0d0
    do j =1, v%ny
    do i =1, v%nx
      uu = v%u(i,j)*v%u(i,j) + v%v(i,j)*v%v(i,j)
      do l =1, v%nq
      eu = model%e(l,1) * v%u(i,j) + model%e(l,2) * v%v(i,j)
      s  = model%w(l) * (1.0d0 + 3.0d0 * eu + 4.5d0 * eu**2 -1.5d0 * uu)
      v%feq(i,j,l) = v%rho(i,j) * s
      enddo
    enddo
    enddo
  endsubroutine calcfeq

  subroutine streaming(v, model)
    implicit none
    type(field)    :: v
    type(lbmmodel) :: model
    integer        :: i,j,l

    v%f(:,:,1) = v%f(:,:,1)
    do j =1, v%ny
    do i =1, v%nx
    do l =2, v%nq
      if( (i+int(model%e(l,1)))<1    .or. &
          (i+int(model%e(l,1)))>v%nx .or. &
          (j+int(model%e(l,2)))<1    .or. &
          (j+int(model%e(l,2)))>v%ny )then
          !write(*,*) i,j
      else
        v%f(i+int(model%e(l,1)), j+int(model%e(l,2)),l) = v%f(i,j,l)
      endif
    enddo
    enddo
    enddo
  endsubroutine streaming

  subroutine calc_macroscopic(v,model)
    implicit none
    type(field)    :: v
    type(lbmmodel) :: model
    integer        :: i,j,l
    real(8)        :: tmp1, tmp2
   
    do j =1, v%ny
    do i =1, v%nx
      tmp1 = 0.0d0
      do l =1, v%nq
        tmp1 = tmp1 + v%f(i,j,l)
      enddo
      v%rho(i,j) = tmp1
    enddo
    enddo
    ! Zou he
    do i =2, v%nx-1
      v%rho(i,v%ny) = (1.d0 / 1.d0 + v%v(i,v%ny)) &
                    * (v%f(i,v%ny,1) + v%f(i,v%ny,2) + v%f(i,v%ny,4)) &
                    + 2.d0 * (v%f(i,v%ny,3) + v%f(i,v%ny,6) + v%f(i,v%ny,7))
    enddo
 
    do j =2, v%ny-1
    do i =2, v%nx-1
      tmp1 = 0.0d0
      tmp2 = 0.0d0
      do l =1, v%nq
        tmp1 = tmp1 + model%e(1,l) * v%f(i,j,l)
        tmp2 = tmp2 + model%e(2,l) * v%f(i,j,l)
      enddo
      v%u(i,j) = tmp1 / v%rho(i,j)
      v%v(i,j) = tmp2 / v%rho(i,j)
    enddo
    enddo
    endsubroutine calc_macroscopic

    subroutine set_boundary(v,model)
    implicit none
    type(field)    :: v
    type(lbmmodel) :: model
    integer        :: i,j,l
    integer        :: nx, ny
    nx = v%nx
    ny = v%ny
    !! Bounce Back inter
    !do j =2, ny-1
    !do i =2, nx-1
    !do l=2, v%nq
    !  ! v ob      : object
    !  ! model lbb : list BB
    !  if(v%ob(i+int(model%e(l,0)), j+int(model%e(l,1))) * v%ob(i,j) .eq. 0)then
    !    v%f(i,j,model%lbb(l)) = v%f(i,j,l)
    !  endif
    !enddo
    !enddo
    !enddo
    ! Bounce Back Outer
    ! Velocity
    do i =1, nx
      ! South
      v%f(i,1,model%lbb(5)) = v%f(i,1,5)
      v%f(i,1,model%lbb(8)) = v%f(i,1,8)
      v%f(i,1,model%lbb(9)) = v%f(i,1,9)
      v%u(i,1) = 0.d0
      v%v(i,1) = 0.d0
      !v%f(i,1,5) = 0.0d0 
      !v%f(i,1,8) = 0.0d0
      !v%f(i,1,9) = 0.0d0
    enddo
    do j =1, ny
      ! West
      v%f(1,j,model%lbb(4)) = v%f(1,j,4)
      v%f(1,j,model%lbb(7)) = v%f(1,j,7)
      v%f(1,j,model%lbb(8)) = v%f(1,j,8)
      v%u(1,j) = 0.d0
      v%v(1,j) = 0.d0
      !v%f(1,j,4) = 0.0d0 
      !v%f(1,j,7) = 0.0d0
      !v%f(1,j,8) = 0.0d0
      ! East
      v%f(nx,j,model%lbb(2)) = v%f(nx,j,2)
      v%f(nx,j,model%lbb(6)) = v%f(nx,j,6)
      v%f(nx,j,model%lbb(9)) = v%f(nx,j,9)
      v%u(nx,j) = 0.d0
      v%v(nx,j) = 0.d0
      !v%f(nx,j,2) = 0.0d0 
      !v%f(nx,j,6) = 0.0d0
      !v%f(nx,j,9) = 0.0d0
    enddo
 
    do i =2, nx-1
    ! North
    v%u(i,ny) = 1.d0
    v%v(i,ny) = 0.d0
    !v%rho(i,ny) = (1.d0 / 1.d0 + v%v(i,ny)) * (v%f(i,ny,1) + v%f(i,ny,2) + v%f(i,ny,4)) &
    !                                 + 2.d0 * (v%f(i,ny,3) + v%f(i,ny,6) + v%f(i,ny,7))
 
    v%f(i,ny,model%lbb(3)) = v%f(i,ny,3)
    v%f(i,ny,model%lbb(6)) = v%f(i,ny,6) &
                           - v%rho(i, ny) * v%u(i,ny) / 0.6d0
    v%f(i,ny,model%lbb(7)) = v%f(i,ny,7) &
                           + v%rho(i, ny) * v%u(i,ny) / 0.6d0
    !v%f(i,ny,3) = 0.0d0 
    !v%f(i,ny,6) = 0.0d0
    !v%f(i,ny,7) = 0.0d0
    enddo
    endsubroutine

 
end module basic

program main
  use basic
  implicit none
  integer :: i, j, l, t
  type(field)    :: v
  type(lbmmodel) :: model
  type(solver)   :: sol
  character(20) :: filename
  real(8)       :: tmp1
  real(8)       :: tmp2

v = nfield(121, 121, 9, 1.0d0, 1.0d0)
model = d2q9(2, 9)
sol%ntime = 1000
sol%Re = 100.0d0
sol%nu = 0.064d0
sol%tau= 0.69d0

do l = 1, v%nq
  v%f(:,:,l) = model%w(l)
enddo
v%rho(:,:)  = 1.d0
v%u(2:v%nx-1,v%ny) = 1.d0

do t = 1, sol%ntime
  if(mod(t,100).eq.0)then
    tmp1 = 0.0d0
    tmp2 = 99999.d0
    do j =1, v%ny
    do i =1, v%nx
     tmp1 = max(v%rho(i,j), tmp1)
     tmp2 = min(v%rho(i,j), tmp2)
    enddo
    enddo
    write(IUT6, '(a10, i6)') "Iteration", t
    write(IUT6, *) tmp1,tmp2
  endif
  ! Collision
  call calcfeq(v, model)
  do j =1, v%ny
  do i =1, v%nx
  do l =1, v%nq
    v%f(i,j,l) = v%f(i,j,l) - (1.0d0 / sol%tau) * (v%f(i,j,l) - v%feq(i,j,l))
  enddo
  enddo
  enddo
 ! Streaming
  call streaming(v,model)
  ! BB & Boundary
  call set_boundary(v,model)
 ! macroscopic
  call calc_macroscopic(v, model)

enddo

! for gnuplot data
filename="u.dat"
call push_gnuplot(v, v%u, filename)
filename="v.dat"
call push_gnuplot(v, v%v, filename)
filename="rho.dat"
call push_gnuplot(v, v%rho, filename)

contains
  subroutine push_gnuplot(v, values, filename)
  implicit none
  integer     :: i,j
  type(field) :: v
  real(8)     :: values(:,:)
  character(20) :: filename
  ! save data
  open(10,file=trim(filename))
  do i = 1, v%nx
    do j = 1, v%ny
      write(10,*) v%x(i), v%y(j), values(i, j)
    enddo
      write(10,*)
  enddo
  close(10)
  end subroutine push_gnuplot
end program main
