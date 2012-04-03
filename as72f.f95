	program AS72f
	implicit none

	include 'mpif.h'
	integer status(MPI_STATUS_SIZE), ierror, ReqID
	integer nprocs, mypid, other, tag
	integer nn
	parameter (nn=4)

	real*4 a(nn,nn),b(nn,nn)
	integer N,M,bndry1,bndry2
	integer i, j, iter

	write(*,'(a)')' AS72/P1 33-241, F08 GREER'

	call MPI_Init(ierror)
	call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierror)
	call MPI_Comm_rank(MPI_COMM_WORLD, mypid, ierror)

	if (mypid.eq.0) then
		other=1
		N=2
		M=nn/2
		bndry1=nn/2
		bndry2=nn/2+1
	else
		other=0
		N=(nn/2)+1
		M=nn-1
		bndry1=nn/2+1
		bndry2=nn/2
	endif

! init arrays
	a=0.25 ; a(:,1)=1. ; a(:,nn)=0. ; a(1,:)=0. ; a(nn,:)=0. ; b=a

	tag = 1
	do iter=1,20
		a(2:nn-1,N:M)=(b(3:nn,N:M)+b(1:nn-2,N:M)+b(2:nn-1,N+1:M+1)+b(2:nn-1,N-1:M-1))/4

		call MPI_ISEND(a(1,bndry1), nn, MPI_REAL, other, tag, MPI_COMM_WORLD, ReqID, ierror)
		call MPI_RECV(a(1,bndry2), nn, MPI_REAL, other, tag, MPI_COMM_WORLD, status, ierror)
		call MPI_WAIT(ReqID, status, ierror)

		b(1:nn,N-1:M+1)=a(1:nn,N-1:M+1)
	enddo

	do j = N-1, M+1
		write(10+mypid,'(xa,i2,i3,1x,4f6.3)')   'id',mypid,j,(a(i,j),i=1,nn)
	enddo

	do j=nn/2+1, nn
		if (mypid.ne.0) then
			call MPI_SEND(a(1,j),nn,MPI_REAL,other,tag,MPI_COMM_WORLD,ierror)
		else
			call MPI_RECV(a(1,j),nn,MPI_REAL,other,tag,MPI_COMM_WORLD,status,ierror)
		endif
	enddo
	call MPI_BARRIER(MPI_COMM_WORLD, ierror)

	if (mypid.eq.0) then
		write(10+mypid,'(/x4f7.4)') ((a(i,j),j=1,4),i=1,4)
	endif

	call MPI_FINALIZE(ierror)

	end program AS72f




