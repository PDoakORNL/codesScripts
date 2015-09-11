program toband
implicit none


integer :: nkpt,nband,nspin,nlines,option
integer :: natom,nelectron,d1,d2
integer :: i,j,k
integer :: index_buffer
real*8 :: energy_buffer
real*8 :: a_x,a_y,a_z
real*8 :: b_x,b_y,b_z
real*8 :: c_x,c_y,c_z
real*8,allocatable :: kccordinates_x(:),kccordinates_y(:),kccordinates_z(:)
real*8,allocatable :: eigens_up(:,:),eigens_dw(:,:)
real*8,allocatable :: eigens_tmp(:)
real*8 :: E_f,bandgap,E_c,E_v
character*80 :: filename,line
character*60 ::  char1, char2, char3
logical :: existence, odd_electron_number


open(unit=7,file="KPOINTS")
open(unit=9,file="EIGENVAL")
read(7,*) nlines
read(9,*) natom,d1,d2,nspin

do i=1,4
  read(9,*) 
end do

read(9,*) nelectron,nkpt,nband

allocate(kccordinates_x(nkpt),kccordinates_y(nkpt),kccordinates_z(nkpt))
allocate(eigens_up(nband,nkpt),eigens_dw(nband,nkpt))
allocate(eigens_tmp(2*nband))


do j=1,nkpt
  read(9,*) 
  read(9,*) 
  do i=1,nband
    if(nspin==1) then
      read(9,*) line,eigens_up(i,j)
    else
      read(9,*) line,eigens_up(i,j),eigens_dw(i,j)
    end if
  end do
end do


open(unit=11,file="bandgap")

if(mod(nelectron,2)==0) then
if(nspin==1) then
  inquire(file="option",exist=existence)
  if(existence) then
    open(unit=1,file="option")
    read(1,*) option
    close(1)
  else
    write(*,*) "Please input 0 if this is noncollinear calculations, otherwise 1: "
    read(*,*) option
  end if
  if(option==1) nelectron=nelectron/2
  E_c=100000000
  E_v=-100000000
  do j=1,nkpt
    do i=1,nband
      eigens_tmp(i)=eigens_up(i,j)
    end do
  do i=1,(nband-1)
    energy_buffer=eigens_tmp(i)
    index_buffer=i
    do k=i+1,nband
      if(eigens_tmp(k) .lt. energy_buffer) then
        energy_buffer=eigens_tmp(k)
        index_buffer=k
      end if
    end do
    if(index_buffer .ne. i) then
      eigens_tmp(index_buffer)=eigens_tmp(i)
      eigens_tmp(i)=energy_buffer
    end if
  end do
  if(eigens_tmp(nelectron) .gt. E_v) E_v=eigens_tmp(nelectron)
  if(eigens_tmp(nelectron+1) .lt. E_c) E_c=eigens_tmp(nelectron+1)
  end do

  E_f=(E_c+E_v)/2.0
  bandgap=E_c-E_v
  write(*,*) "The maximum valence band is ",E_v
  write(*,*) "The minimum conduction band is ",E_c
  write(*,*) "The Fermi level is ", E_f
  write(*,*) "The bandgap is ", bandgap
  write(11,"(4F7.2)") E_f, bandgap, E_v, E_c
  
  do i=1,nband
    do j=1,nkpt
       eigens_up(i,j)=eigens_up(i,j)-E_f
    end do
  end do
  do k=1,nlines
     write(filename,'(a,I0,a)') 'eigenvalue',k,' '
     open(unit=k,file=filename)
  end do
  do i=1,nband
    do k=1,nlines
      do j=nkpt/nlines*(k-1)+1,nkpt/nlines*(k)
        write(k,*) j-nkpt/nlines*(k-1),eigens_up(i,j)
      end do
      write(k,*) 
      write(k,*) 
    end do
  end do
  do i=1,nband
    do k=1,nlines
      do j=nkpt/nlines*(k-1)+1,nkpt/nlines*(k)
        write(k,*) j-nkpt/nlines*(k-1),eigens_up(i,j)    !!!!!! bands are degenerate in non spin case
      end do
      write(k,*) 
      write(k,*) 
    end do
  end do
  do k=1,nlines
     close(unit=k)
  end do
else
  E_c=100000000
  E_v=-100000000
  do j=1,nkpt
    do i=1,nband
      eigens_tmp(i)=eigens_up(i,j)
      eigens_tmp(i+nband)=eigens_dw(i,j)
    end do
  do i=1,(2*nband-1)
    energy_buffer=eigens_tmp(i)
    index_buffer=i
    do k=i+1,2*nband
      if(eigens_tmp(k) .lt. energy_buffer) then
        energy_buffer=eigens_tmp(k)
        index_buffer=k
      end if
    end do
    if(index_buffer .ne. i) then
      eigens_tmp(index_buffer)=eigens_tmp(i)
      eigens_tmp(i)=energy_buffer
    end if
  end do
  if(eigens_tmp(nelectron) .gt. E_v) E_v=eigens_tmp(nelectron)
  if(eigens_tmp(nelectron+1) .lt. E_c) E_c=eigens_tmp(nelectron+1)
  end do

  E_f=(E_c+E_v)/2.0
  bandgap=E_c-E_v
  write(*,*) "The maximum valence band is ",E_v
  write(*,*) "The minimum conduction band is ",E_c
  write(*,*) "The Fermi level is ", E_f
  write(*,*) "The bandgap is ", bandgap
  write(11,"(4F7.2)") E_f, bandgap, E_v, E_c
  
  do i=1,nband
    do j=1,nkpt
       eigens_up(i,j)=eigens_up(i,j)-E_f
       eigens_dw(i,j)=eigens_dw(i,j)-E_f
    end do
  end do
  do k=1,nlines
     write(filename,'(a,I0,a)') 'eigenvalue',k,' '
     open(unit=k,file=filename)
  end do
  do i=1,nband
    do k=1,nlines
      do j=nkpt/nlines*(k-1)+1,nkpt/nlines*(k)
        write(k,*) j-nkpt/nlines*(k-1),eigens_up(i,j)
      end do
      write(k,*) 
      write(k,*) 
    end do
  end do
  do i=1,nband
    do k=1,nlines
      do j=nkpt/nlines*(k-1)+1,nkpt/nlines*(k)
        write(k,*) j-nkpt/nlines*(k-1),eigens_dw(i,j)
      end do
      write(k,*) 
      write(k,*) 
    end do
  end do
 
   do k=1,nlines
     close(unit=k)
   end do
end if
else
  write(*,*) "Caution: number of electrons is odd and fermi level used is from VASP itself!!!"
!!!! Fermi level is thus determined from VASP itself
  open(unit=1,file="OUTCAR")

  do
    read(1,*,end=10)
  end do

10 backspace(unit=1)

100 do     
      backspace(unit=1)     
      read(1,"(A)") line 
      line=trim(line)
      line=adjustl(line)  
      backspace(unit=1)
      read(line,*,end=100) char1 !!! in case we have an empty line
      if(char1=="E-fermi") then     
        exit     
      end if     
    end do

  read(line,*) char1,char1, E_f

  close(1)

  write(*,*) "The Fermi level is ", E_f

!   do i=1,nband
!     do j=1,nkpt
!        eigens_up(i,j)=eigens_up(i,j)-E_f
!     end do
!   end do

  do k=1,nlines
     write(filename,'(a,I0,a)') 'eigenvalue',k,' '
     open(unit=k,file=filename)
  end do
  do i=1,nband
    do k=1,nlines
      do j=nkpt/nlines*(k-1)+1,nkpt/nlines*(k)
        write(k,*) j-nkpt/nlines*(k-1),eigens_up(i,j)
      end do
      write(k,*) 
      write(k,*) 
    end do
  end do
  do i=1,nband
    do k=1,nlines
      do j=nkpt/nlines*(k-1)+1,nkpt/nlines*(k)
        write(k,*) j-nkpt/nlines*(k-1),eigens_up(i,j)    !!!!!! bands are degenerate in non spin case
      end do
      write(k,*) 
      write(k,*) 
    end do
  end do
  do k=1,nlines
     close(unit=k)
  end do  
  
end if

close(unit=7)
close(unit=9)
close(11)
 





end program toband

