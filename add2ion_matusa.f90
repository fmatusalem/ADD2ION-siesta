
Program add2ion
   implicit none
   character(20):: Vatomo, Vion, ION, IONnovo, INPUT,rfile,date
   character(40):: formato, arg,str,str1
   character(120)::linha
   character(255) :: dummy,dm,species_array(2550) = '',core,orb,a
   integer:: pula,pula1,pula2,nupon,i,n,nION,ios,nr_lines, nr_words, stats, valence, orbital
   real(8)::r(2001),patomo(2001),pion(2001),VAR(2001),corte,amplitude,ponto,Vna

!Filipe Matusalem adds

CALL get_DDMonYY(date)

if (iargc() < 2) then
   write(*,*) 'ERROR!! Please give the name of .ion file and CUT as arguments!!'
   stop
   end if

call getarg(1,ION)
call getarg(2,str)

if (iargc() > 2) then
   call getarg(3,str1)
   i = scan(trim(ION),".", BACK= .true.)
   IONnovo = trim(ION(1:i-1))//'cut'//trim(str)//'A'//trim(str1)//'.ion'
   read(str1,*)amplitude  !convert string to real
else
   amplitude=1
   i = scan(trim(ION),".", BACK= .true.)
   IONnovo = trim(ION(1:i-1))//'cut'//trim(str)//'.ion'
end if

read(str,*)corte  !convert string to real


write(*,*)ION,IONnovo,corte,amplitude


   open(unit=12,file='add2ion.out',iostat=ios) 
	if (ios/=0) then 
	write(*,*)'Impossible to create add2ion.out'
	stop
        end if

pula=1
n=8
rfile='VTOTAL1.ae'
Vatomo='VTOTAL1.ae'
Vion='VTOTAL1.ae-05'
formato='(2g22.12)'
nION=500
   
   open(unit=1,file=rfile, status='old',iostat=ios) 
	if (ios/=0) then 
	write(*,*)'Impossible to open file  ',rfile
	stop
        end if

nr_lines = 0  ! number of lines in the file
nr_words = 0  ! number of words found

do
   read(1,'(255a)',iostat=ios) dummy
   ! end of file
   if (ios/=0) exit

   dummy=trim(adjustl(dummy)) 
   ! for debugging uncomment the following line
   !write(*,*) dummy
   nr_lines = nr_lines + 1
   call process_line(trim(dummy), nr_lines, nr_words, species_array, n)
   nupon=nr_words-8
   pula1=nr_lines+1
   if (dummy(1:4)=='Down') then
   exit
   end if
end do
   close (1)

!   write(*,*) nupon, pula1

open(unit=1,file=ION, status='old',iostat=ios) 
	if (ios/=0) then 
	write(*,*)'Impossible to open file  ',ION
	stop
        end if

do
   read(1,'(255a)',iostat=ios) dummy
   i=i+1
   ! end of file
   if (ios/=0) exit

   dummy=trim(adjustl(dummy)) 
   ! for debugging uncomment the following line
   !write(*,*) dummy
   if (dummy(1:13)=='<basis_specs>') then
      read(1,*) dm,dm,dm,dm,dm,dm,dm,dm
      exit
   end if
end do
   close (1)

read(dm,*)ponto  !convert string to real
 if (ponto/=0) then
 write(*,*) 'ERROR!! Please use siesta version 2.0.2 to generate .ion file.'
 write(*,*) 'However, the DFT-1/2 corrected .ion file can be used in more recent versions!!'
 stop
 end if


open(unit=1,file=ION, status='old',iostat=ios) 
	if (ios/=0) then 
	write(*,*)'Impossible to open file  ',ION
	stop
        end if

i=0
do
   read(1,'(255a)',iostat=ios) dummy
   i=i+1
   ! end of file
   if (ios/=0) exit

   dummy=trim(adjustl(dummy)) 
   ! for debugging uncomment the following line
   !write(*,*) dummy
   if (dummy(1:6)=='# Vna:') then
   exit
   end if
   pula2=i+2
end do
   close (1)

!write(*,*) pula2

   write(12,'(2a)')rfile,' = nome do arquivo de raios'
   write(12,'(2i5,a)')PULA,nupon,'           = PULA,nupon'
   open(unit=11,file=rfile, status='old',iostat=ios) 
	if (ios/=0) then 
	write(*,*)'Impossible to open file  ',rfile
	stop
        end if
   do i=1,pula;read(11,*);enddo
   read(11,*)(r(i),i=1,nupon)
   close (11)

   write(12,'(2a)')Vatomo,' = nome do arquivo do potencial do atomo'
   open(unit=11,file=Vatomo, status='old',iostat=ios) 
	if (ios/=0) then 
	write(*,*)'Impossible to open file  ',Vatomo
	stop
        end if
   write(12,'(i5,a)')pula1,'                = PULA'
   do i=1,pula1;read(11,*);enddo
   read(11,*)(patomo(i),i=1,nupon)
   close (11)

   write(12,'(2a)')Vion,' = nome do arquivo do potencial do ion'
   open(unit=11,file=Vion, status='old',iostat=ios) 
	if (ios/=0) then 
	write(*,*)'Impossible to open file  ',Vion
	stop
        end if
   write(12,'(i5,a)')pula1,'                = PULA'
   do i=1,pula1;read(11,*);enddo
   read(11,*)(pion(i),i=1,nupon)
   close (11)

   write(12,'(i4,2f10.6,a)')n,corte,amplitude,' = n   corte  amplitude'
   do i=1,nupon
      VAR(i)=0.d0
      if(r(i).lt.corte) &
         VAR(i)=(1.d0-r(i)**n/corte**n)**3*&
          (pion(i)-patomo(i))*amplitude/r(i)
   enddo

   write(12,'(2a)')ION,' = nome do arquivo ION'
   open(unit=11,file=ION, status='old',iostat=ios) 
	if (ios/=0) then 
	write(*,*)'Impossible to open file  ',ION
	stop
        end if
   write(12,'(i5,a)')pula2,'                = PULA'
   write(12,'(2a)')IONnovo,' = arquivo IONnovo'
   open(unit=13,file=IONnovo,iostat=ios) 
	if (ios/=0) then 
	write(*,*)'Impossible to open file  ',IONnovo
	stop
        end if

!---------------------read INP.ae files---------------
   open(unit=19,file='INP.ae-05', status='old',iostat=ios) 
	if (ios/=0) then 
	write(*,*)'Impossible to open file  INP.ae-05'
	stop
        end if

   do i=0,4
   read(19,'(a)')linha
   enddo 

   read(19,*)dm,dm
   read(dm,*)valence  !convert string to int

str1=''
a=''
   do i=0,valence-1
   read(19,*)core
   backspace(19)
   read(19,*)dm,dm;read(dm,*)orbital
   backspace(19)
   read(19,*)dm,dm,dm

   if(orbital==0) orb='s'
   if(orbital==1) orb='p'
   if(orbital==2) orb='d'
   if(orbital==3) orb='f'

   str = '  '//trim(core)//trim(orb)//trim(dm)
   str1 = trim(a)//str
   a= str1
   enddo
   write(*,*)str1
   close(19)
!-----------------------------------------------------


   read(11,'(a)')linha;write(13,'(a)')TRIM(linha)
   write(13,'(a,f5.3,a,a,f6.3,a,a7)')'# DFT-1/2  CUT=',corte, trim(str1),' Amp.=',amplitude,' Date: ',date
   read(11,'(a)')linha
   do i=1,pula2-2;read(11,'(a)')linha;write(13,'(a)')TRIM(linha);enddo

   write(12,'(2a)')formato,' = formato no arquivo ION'
   write(12,'(i5,a)')nION,'                = numero de raios e potenciais Vna no arquivo ION'
   close(12)
   open(unit=21,file='trimmed')
   do i=1,nupon;write(21,'(2f20.10)')r(i),VAR(i);enddo
   close(21)

! Aqui soma a interpolaccao de VAR e escreve
   do i=1,nION;read(11,*)ponto,Vna
      Vna=Vna+Deriva_Interpola(.False.,nupon,r,VAR,6,ponto)
      if(i.ne.1)then
         if(i.eq.2)write(13,formato)0.,Vna
         write(13,formato)ponto,Vna
      endif
   enddo
!
! Escreve o resto
100 read(11,'(a)',end=200,err=200)linha
    write(13,'(a)')TRIM(linha)
    go to 100
200 stop
contains

FUNCTION Deriva_Interpola(deriva,npon,x,y,nproximo,ponto)
! Dados os arrays x(1:npon) onde é conhecida a função y(1:npon), esta
! subrotina interpola ou deriva a funccao y(x) no ponto 'ponto', usando
! os 'nproximo' pontos mais proximos de ponto.
! Se deriva=.TRUE. ==> deriva a funccao.
! Se deriva=.FALSE. ==> interpola a funccao.
   IMPLICIT none
   LOGICAL :: deriva
   INTEGER :: npon,nproximo,i,j,k
   REAL(8) :: ponto,x(npon),y(npon),Deriva_Interpola,array(1:npon),&
            slave(2,1:npon),aux,fator,bux
!INTERFACE
!   subroutine ordena(narray,nslave,array,slave)
!   integer::narray,nslave
!   real(8)::array(:),slave(:,:)
!   end subroutine ordena
!END INTERFACE
   array(1:npon)=(x(1:npon)-ponto)**2
   slave(1,1:npon)=x(1:npon)
   slave(2,1:npon)=y(1:npon)
   call ordena(npon,2,array,slave)
   if(deriva)then
      aux=0.d0
      do i=1,nproximo
         bux=0.d0
         do j=1,nproximo
            if(j.ne.i)then
               fator=1.d0
               do k=1,nproximo
                  if(k.ne.i.and.k.ne.j)then
                     fator=fator*(ponto-slave(1,k))/(slave(1,i)-slave(1,k))
                  endif
               enddo
               bux=bux+fator/(slave(1,i)-slave(1,j))
            endif
         enddo
         aux=aux+bux*slave(2,i)
      enddo
   else
      aux=0.d0
      do i=1,nproximo
         fator=1.d0
         do j=1,nproximo
            if(i.ne.j)then
               fator=fator*(ponto-slave(1,j))/(slave(1,i)-slave(1,j))
            endif
         enddo
         aux=aux+slave(2,i)*fator
      enddo
   endif
   Deriva_Interpola=aux
END FUNCTION Deriva_Interpola

subroutine ordena(narray,nslave,array,slave)
!  Subrotina para ordenar o array array crescentemente. Os arrays slave sao
!  ordenados da mesma maneira.
   implicit none
   integer::narray,nslave,i,j,k
   real(8)::array(1:narray),slave(1:nslave,1:narray),replace
   do i=1,narray-1
      do j=i+1,narray
         if(array(j).ge.array(i))cycle
         replace=array(i)
         array(i)=array(j)
         array(j)=replace
         do k=1,nslave
            replace=slave(k,i)
            slave(k,i)=slave(k,j)
            slave(k,j)=replace
         enddo
      enddo
   enddo
end subroutine ordena
end Program add2ion

subroutine process_line(line, i, k, array, n)
  ! line  = input line
  ! i     = line index
  ! k     = word index index in the words array
  ! array =  array of words
  ! n     = dimension of the words array
  character(*), intent(in) :: line
  integer, intent(in) :: i, n 
  integer, intent(in out) :: k
  character(*), dimension(n), intent(out) :: array 
  ! temp vars
  character(20) :: word
  integer :: j, beg_idx, end_idx

  !write(*,*) i, ".line: '", line, "'"
  beg_idx = 0 ! begin index of the word
  end_idx = 0 ! end index of the word
  j = 1       ! character position in the line
  do while (j <= len(line))
    if ((beg_idx .eq. 0) .and. (line(j:j) .ne. ' ')) then
       beg_idx = j
       !write (*,*) 'beg_idx = ', beg_idx       
    end if
    if (beg_idx .gt. 0) then
      if (line(j:j) .eq. ' ') then
        end_idx = j - 1
      end if
      if (j .eq. len(line)) then 
        end_idx = j
      end if
      if (end_idx .gt. 0) then
         !write(*,*) 'end_idx = ', end_idx
         k = k + 1
         word = line(beg_idx:end_idx)
         !write (*,*) '  ', k,".word = '", trim(word),"'"
         array(k) = word
         ! initialize indices
         beg_idx = 0
         end_idx = 0
      end if
    end if
    j = j + 1  
  end do     
end subroutine process_line

SUBROUTINE get_DDMonYY(date)
    CHARACTER(len=7), INTENT(out) :: date

    CHARACTER(len=2) :: dd
    CHARACTER(len=3) :: mons(12)
    CHARACTER(len=4) :: yyyy
    INTEGER :: values(8)

    mons = ['Jan','Feb','Mar','Apr','May','Jun',&
      'Jul','Aug','Sep','Oct','Nov','Dec']

    CALL DATE_AND_TIME(VALUES=values)

    WRITE(  dd,'(i2)') values(3)
    WRITE(yyyy,'(i4)') values(1)

    date = trim(dd)//trim(mons(values(2)))//trim(yyyy(3:4))
END SUBROUTINE get_DDMonYY


