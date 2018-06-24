subroutine write_vtk2d(im,jm,ntime,x_plt,y_plt,u_plt,v_plt,p_plt)
  implicit none
  integer(4)    :: im, jm, ntime
  integer(4)    :: i,j
  real(8)       :: x_plt(im)
  real(8)       :: y_plt(jm)
  real(8)       :: u_plt(im,jm)
  real(8)       :: v_plt(im,jm)
  real(8)       :: p_plt(im,jm)
  character(30) :: filename


  !! ファイルの出力
  write(filename,"('data',1i6.6,'.vtk')") ntime !! 出力ファイル番号をファイル名に書き込む
  open(200,file=trim(filename),status="unknown",form="formatted",position="rewind")
  
  write(200,"('# vtk DataFile Version 3.0')")
  write(200,"('2D flow')")
  write(200,"('ASCII ')")
  
  write(200,"('DATASET STRUCTURED_GRID')")
  write(200,"('DIMENSIONS ',3(1x,i3))") im, jm, 1
  
  write(200,"('POINTS ',i9,' float')") im*jm
  do j=1,jm
  do i=1,im
    write(200,"(3(f9.4,1x))") x_plt(i), y_plt(j), 0.0d0
  enddo
  enddo
  
  write(200,"('POINT_DATA ',i9)") im*jm
  !! velocity vector
  write(200,"('VECTORS velocity float')")
  do j=1,jm
  do i=1,im
    write(200,"(3(f9.4,1x))") u_plt(i,j), v_plt(i,j), 0.0d0
  enddo
  enddo
  
  !! pressure
  write(200,"('SCALARS pressure float')")
  write(200,"('LOOKUP_TABLE default')")
  do j=1,jm
  do i=1,im
    write(200,"(f9.4)") p_plt(i,j)
  enddo
  enddo
  
  !!! z_vorticity
  !write(200,"('SCALARS z_vorticity float')")
  !write(200,"('LOOKUP_TABLE default')")
  !do j=1,jm
  !do i=1,im
  !  write(200,"(f9.4)") z_vorticity(i,j)
  !enddo
  !enddo
  close(200)

endsubroutine  write_vtk2d
