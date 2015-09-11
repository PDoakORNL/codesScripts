program spectrum
implicit none

integer :: nqpt,npath,natom,nmode
integer :: i,j,k
real*8,allocatable :: q_x(:),q_y(:),q_z(:)
real*8,allocatable :: frequency(:,:)
character*80 :: string,filename

open(unit=7,file="band.yaml")
read(7,*) string, nqpt
read(7,*) string, npath
read(7,*) string, natom
read(7,*) 
read(7,*) 
read(7,*) 
read(7,*) 
read(7,*) 

nmode=3*natom

allocate(q_x(nqpt),q_y(nqpt),q_z(nqpt))
allocate(frequency(nqpt,nmode))

do i=1,nqpt
  read(7,*) string,string,string,q_x(i),q_y(i),q_z(i)
  read(7,*)
  read(7,*)
  do j=1,nmode
    read(7,*)
    read(7,*) string,frequency(i,j)
  end do
  read(7,*)
end do

open(unit=2,file="modes.txt")
do i=1,nmode
  write(2,"(I,A,F6.2,A)") i," f = ",frequency(1,i)*33.356407889," cm-1"
end do 
close(unit=2)


do k=1,npath
   write(filename,'(a,I0,a)') 'frequency',k,' '
   open(unit=k,file=filename)
end do

do j=1,nmode
  do k=1,npath
    do i=nqpt/npath*(k-1)+1,nqpt/npath*(k)
        write(k,*) i-nqpt/npath*(k-1),frequency(i,j)*33.356407889
    end do
    write(k,*) 
    write(k,*) 
  end do
end do
 


end program spectrum



