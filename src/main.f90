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

v = nfield(32,32, 9, 1.0d0, 1.0d0)
!v = nfield(128,128, 9, 1.0d0, 1.0d0)
model = d2q9(2, 9)
! Setting
sol%Re  = 1000.0d0
v%Uwall = 0.05d0
v%Vwall = 0.0d0
sol%c  = 1.0d0
!sol%dt = v%dx / sol%c
sol%dt = v%Uwall * v%dx / sol%c
sol%nu = v%Uwall * v%Lx / sol%Re
sol%tau= 0.50d0 + 3.0d0 * sol%nu / (sol%dt * (sol%c**2))
!sol%tau= 0.50d0 + 3.0d0 * sol%nu

write(IUT6, '(a10, f12.6)') "Re",    sol%Re
write(IUT6, '(a10, f12.6)') "dt",    sol%dt
write(IUT6, '(a10, f12.6)') "Uwall", v%Uwall
write(IUT6, '(a10, f12.6)') "nu",    sol%nu
write(IUT6, '(a10, f12.6)') "tau",   sol%tau

!if(sol%tau > 0.8d0)then
!  write(IUT6, '(a17, f12.6)') "tau is too big :",   sol%tau
!  stop
!endif
sol%ntime = 1000
do l = 1, v%nq
  v%f(:,:,l) = model%w(l)
enddo
v%rho(:,:)  = 1.d0
! miki
do i = 2,v%nx-1
  v%u(i, v%ny)  = v%Uwall
  v%rho(i,v%ny) = (1.d0 / 1.d0 + v%Vwall) &
                       * (v%f(i,v%ny,1) + v%f(i,v%ny,2) + v%f(i,v%ny,4)) &
                + 2.d0 * (v%f(i,v%ny,3) + v%f(i,v%ny,6) + v%f(i,v%ny,7))
enddo

do t = 1, sol%ntime
  if(mod(t,10).eq.0)then
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
    call write_vtk2d(v%nx,v%ny,t,v%x,v%y,v%u,v%v,v%rho/3.0d0)
  endif
  ! Collision
  v%feq = 0.0d0
  call calcfeq(v, model)
  do j =1, v%ny
  do i =1, v%nx
  do l =1, v%nq
    v%f(i,j,l) = v%f(i,j,l) - (1.0d0 / sol%tau) * (v%f(i,j,l) - v%feq(i,j,l))
  enddo
  enddo
  enddo
  ! Streaming
  !call streaming(v,model)
  call streaming2(v,model)
  v%f = v%fs
  !call transition(v,model)
  ! BB & Boundary
  call set_boundary(v,model)
  ! macroscopic
  call calc_macroscopic(v, model)

enddo
write(*,*) model%e(1,1)
write(*,*) model%e(2,1)

! Cavity
open(10, file="cavityXV.txt")
do i = 1,v%nx
  write(10,*) v%x(i), v%v(i,int(v%ny/2))
enddo
close(10)
open(10, file="cavityYU.txt")
do j = 1,v%ny
  write(10,*) v%y(j), v%u(int(v%nx/2),j)
enddo
close(10)
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
