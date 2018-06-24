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
    real(8)             :: Uwall
    real(8)             :: Vwall
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
    real(8) :: Re,nu,tau,dt
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
    write(IUT6,'(a18, 3i5, f10.6)') "allocate nx ny nq:",nfield%nx, nfield%ny, nfield%nq, nfield%Uwall
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
    d2q9%e(1,:) = (/0.0d0, 1.0d0, 0.0d0, -1.0d0, 0.0d0, 1.0d0, -1.0d0, -1.0d0, 1.0d0/)
    d2q9%e(2,:) = (/0.0d0, 0.0d0, 1.0d0, 0.0d0, -1.0d0, 1.0d0, 1.0d0,  -1.0d0,-1.0d0/)
    !opposite of 1,2,3,4,5,6,7,8,9
    d2q9%lbb = (/1,4,5,2,3,8,9,6,7/)
    write(IUT6,'(a19)') "set D2Q9 model Done"
  endfunction d2q9

  subroutine calcfeq(v, model)
    implicit none
    type(field)    :: v
    type(lbmmodel) :: model
    integer        :: i,j,l
    real(8)        :: eu, uu
    do j =1, v%ny
    do i =1, v%nx
      uu = v%u(i,j)*v%u(i,j) + v%v(i,j)*v%v(i,j)
      do l =1, v%nq
        eu = model%e(l,1) * v%u(i,j) + model%e(l,2) * v%v(i,j)
        v%feq(i,j,l) = model%w(l) * v%rho(i,j) * (1.0d0 + 3.0d0 * eu + 4.5d0 * eu*eu -1.5d0 * uu)
      enddo
    enddo
    enddo
  endsubroutine calcfeq

  subroutine transition(v, model)
    implicit none
    type(field)    :: v
    type(lbmmodel) :: model
    integer        :: i,j,l
    integer        :: ii,jj
    integer        :: id,jd
    real(8),allocatable :: ft(:,:)
    allocate(ft(v%nx, v%ny))
    do l=1,v%nq
      id = model%e(1,l)
      jd = model%e(2,l)

      do j = 1, v%ny
      do i = 1, v%nx
        ii = max(min(i-id,v%nx),1)
        jj = max(min(j-jd,v%ny),1)
        ft(i,j) = v%f(ii,jj,l)
      enddo
      enddo
      do j = 1, v%ny
      do i = 1, v%nx
         v%f(ii,jj,l) = ft(i,j)
      enddo
      enddo
 
    enddo
    deallocate(ft)
  endsubroutine transition

 
  subroutine streaming(v, model)
    implicit none
    type(field)    :: v
    type(lbmmodel) :: model
    integer        :: i,j,l

    ! 2 Right
    do j =1, v%ny
      do i =v%nx, 2, -1
        v%f(i,j,2) = v%f(i-1,j,2)
      enddo
      ! 4 Left
      do i =1, v%nx-1
        v%f(i,j,4) = v%f(i+1,j,4)
      enddo
    enddo

    ! 3 Up
    do j =v%ny, 2, -1
      do i =1, v%nx
        v%f(i,j,3) = v%f(i,j-1,3)
      enddo
      ! 6 UpRight
      do i =v%nx,2, -1
        v%f(i,j,6) = v%f(i-1,j-1,6)
      enddo
      do i =1,v%nx-1
        v%f(i,j,7) = v%f(i+1,j-1,7)
      enddo
    enddo

    ! 5 Down
    do j =1, v%ny-1
      do i =1, v%nx
        v%f(i,j,5) = v%f(i,j+1,5)
      enddo
      ! 8 DowmnLeft
      do i =1, v%nx-1
        v%f(i,j,8) = v%f(i+1,j+1,8)
      enddo
      ! 9 DownRight
      do i =v%nx, 2, -1
        v%f(i,j,9) = v%f(i-1,j+1,9)
      enddo
    enddo

    !! 6 UpRight
    !do j =v%ny,2, -1
    !do i =v%nx,2, -1
    !  v%f(i,j,6) = v%f(i-1,j-1,6)
    !enddo
    !enddo
    !! 7 UpLeft
    !do j =v%ny, 2, -1
    !do i =1,v%nx-1
    !  v%f(i,j,7) = v%f(i+1,j-1,7)
    !enddo
    !enddo

    !! 8 DowmnLeft
    !do j =1, v%ny-1
    !do i =1, v%nx-1
    !  v%f(i,j,8) = v%f(i+1,j+1,8)
    !enddo
    !enddo
    !! 9 DownRight
    !do j =1, v%ny-1
    !do i =v%nx, 2, -1
    !  v%f(i,j,9) = v%f(i-1,j+1,9)
    !enddo
    !enddo

    !do j =1, v%ny
    !do i =1, v%nx
    !do l =2, v%nq
    !  if( (i+int(model%e(l,1)))<1    .or. &
    !      (i+int(model%e(l,1)))>v%nx .or. &
    !      (j+int(model%e(l,2)))<1    .or. &
    !      (j+int(model%e(l,2)))>v%ny )then
    !      !write(*,*) i,j
    !  else
    !    v%f(i+int(model%e(l,1)), j+int(model%e(l,2)),l) = v%f(i,j,l)
    !  endif
    !enddo
    !enddo
    !enddo
    ! Lef
    !do j =1, v%ny
    !do i =1, v%nx
    !      (i+int(model%e(l,1)))>v%nx .or. &
    !      (j+int(model%e(l,2)))<1    .or. &
    !      (j+int(model%e(l,2)))>v%ny )then
    !      !write(*,*) i,j
    !  else
    !    v%f(i+int(model%e(l,1)), j+int(model%e(l,2)),l) = v%f(i,j,l)
    !  endif
    !enddo
    !enddo
    !enddo

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
    ! Zou-He
    do i =1, v%nx
      v%rho(i,v%ny) = (1.d0 / 1.d0 + v%Vwall) &
                           * (v%f(i,v%ny,1) + v%f(i,v%ny,2) + v%f(i,v%ny,4)) &
                    + 2.d0 * (v%f(i,v%ny,3) + v%f(i,v%ny,6) + v%f(i,v%ny,7))
    enddo
 
    do j =1, v%ny
    do i =1, v%nx
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
    ! Bottom
    v%u(:,1) = 0.0d0
    v%v(:,1) = 0.0d0
    ! East
    v%u(v%nx,:) = 0.0d0
    v%v(v%nx,:) = 0.0d0
    ! West
    v%u(1,:) = 0.0d0
    v%v(1,:) = 0.0d0
    ! Top
    ! miki
    do i =1, v%nx
      v%u(i,v%ny) = v%Uwall
      v%v(i,v%ny) = v%Vwall
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

    ! Bounce Back
    ! South
    do i =1, nx
      v%f(i,1,model%lbb(5))  = v%f(i,1,5)
      v%f(i,1,model%lbb(8))  = v%f(i,1,8)
      v%f(i,1,model%lbb(9))  = v%f(i,1,9)
    enddo
    do j =1, ny
    ! West
      v%f(1,j,model%lbb(4))  = v%f(1,j,4)
      v%f(1,j,model%lbb(7))  = v%f(1,j,7)
      v%f(1,j,model%lbb(8))  = v%f(1,j,8)
    ! East
      v%f(nx,j,model%lbb(2)) = v%f(nx,j,2)
      v%f(nx,j,model%lbb(6)) = v%f(nx,j,6)
      v%f(nx,j,model%lbb(9)) = v%f(nx,j,9)
    enddo
 
    ! North
    do i =1, nx
      ! BB
      !v%f(i,ny,model%lbb(3)) = v%f(i,ny,3)
      !v%f(i,ny,model%lbb(6)) = v%f(i,ny,6)
      !v%f(i,ny,model%lbb(7)) = v%f(i,ny,7)
      ! Zou-He
      v%rho(i,ny) = (1.d0 / 1.d0 + v%Vwall) * (v%f(i,ny,1) + v%f(i,ny,2) + v%f(i,ny,4)) &
                                     + 2.d0 * (v%f(i,ny,3) + v%f(i,ny,6) + v%f(i,ny,7))
      v%f(i,ny,model%lbb(3)) = v%f(i,ny,3) - 2.0d0 * v%Vwall / 3.0d0
      v%f(i,ny,model%lbb(6)) = v%f(i,ny,6)  - v%Vwall / 3.0d0 &
                             +(v%f(i, ny,2) - v%f(i,ny,4))/ 2.0d0 &
                             - v%rho(i, ny) * v%Uwall / 2.0d0 &
                             - v%rho(i, ny) * v%Vwall / 2.0d0
      v%f(i,ny,model%lbb(7)) = v%f(i,ny,7)  + v%Vwall / 3.0d0 &
                             -(v%f(i, ny,2) - v%f(i,ny,4))/ 2.0d0 &
                             + v%rho(i, ny) * v%Uwall / 2.0d0 &
                             - v%rho(i, ny) * v%Vwall / 2.0d0
    enddo
    endsubroutine set_boundary
end module basic
