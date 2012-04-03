	program AS73f
	implicit none

	include 'mpif.h'
	integer status(MPI_STATUS_SIZE), ierror, ReqID
	integer nprocs, mypid, other, tag
	integer nn
	parameter (nn=20)

	real*4 a(0:20,0:20),b(0:20,0:20), aed(19), aec(19), dad, dac, rms, dad2, dac2, rms2
	integer N,M,bndry1,bndry2
	integer i, j, iter, s

	aed = (/.4972559,.4890261,.4753282,.4562242,.4318678,.4025638,.3688278,.3314301,.2914033,.2499993,.2085953,.1685686,.1311710, & 
		.0974351,.0681313,.0437752,.0246714,.0109737,.0027441/)
	aec = (/.4972645,.4890595,.4753997,.4563411,.4320283,.4027535,.3690196,.3315883,.2914935,.2500000,.2085065,.1684117,.1309804, & 
		.0972465,.0679717,.0436589,.0246003,.0109505,.00273552/)
	dad = 0.
	dac = 0.

	write(*,'(a)')' AS73/P1 33-241, F08 GREER'

	call MPI_Init(ierror)
	call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierror)
	call MPI_Comm_rank(MPI_COMM_WORLD, mypid, ierror)

	s = 10

	if (mypid.eq.0) then
		other=1
		N=1
		M=s
		bndry1=s
		bndry2=s+1
	else
		other=0
		N=s+1
		M=19
		bndry1=s+1
		bndry2=s
	endif

! init arrays
	a=0.25 ; a(:,0)=1. ; a(:,20)=0. ; a(0,:)=0. ; a(20,:)=0.
	b=0.25 ; b(:,0)=1. ; b(:,20)=0. ; b(0,:)=0. ; b(20,:)=0.

	tag = 1
	do iter=1,20
		dad=0.
		dad2=0.
		dac=0.
		dac2=0.
		rms=0.
		rms2=0.
		!write(10+mypid,'(i4,a,i4,a,21f6.3)'),mypid,' iter=',iter,'    ',(a(10,j),j=N-1,M+1)

		!a(1:19,N:M)=(b(2:20,N:M)+b(0:18,N:M)+b(1:19,N+1:M+1)+b(1:19,N-1:M-1))/4
		do i=1,19
			do j=N,M
				a(i,j) = (b(i+1,j) + b(i-1,j) + b(i,j+1) + b(i,j-1))/4.
			enddo
		enddo

		rms = SUM((a(1:19,N:M)-b(1:19,N:M))**2)
		if (mypid.ne.0) then
			do j=s+1,19
				dad = dad + (a(j,j)-aed(j))**2.
				dac = dac + (a(j,j)-aec(j))**2.
				!write(*,'(i4,a,i4,2f9.3)'),mypid,' debug=',j,a(j,j),aec(j)
			enddo
		else
			do j=1,s
				dad = dad + (a(j,j)-aed(j))**2.
				dac = dac + (a(j,j)-aec(j))**2.
				!write(*,'(i4,a,i4,2f9.3)'),mypid,' debug=',j,a(j,j),aec(j)
			enddo
		endif

		call MPI_ISEND(a(0,bndry1), 20, MPI_REAL, other, tag, MPI_COMM_WORLD, ReqID, ierror)
		call MPI_RECV(a(0,bndry2), 20, MPI_REAL, other, tag, MPI_COMM_WORLD, status, ierror)
		call MPI_WAIT(ReqID, status, ierror)

		if (mypid.ne.0) then
			call MPI_SEND(dac,1,MPI_REAL,other,2,MPI_COMM_WORLD,ierror)
			call MPI_SEND(dad,1,MPI_REAL,other,3,MPI_COMM_WORLD,ierror)
			call MPI_SEND(rms,1,MPI_REAL,other,4,MPI_COMM_WORLD,ierror)
		else
			call MPI_RECV(dac2,1,MPI_REAL,other,2,MPI_COMM_WORLD,status,ierror)
			call MPI_RECV(dad2,1,MPI_REAL,other,3,MPI_COMM_WORLD,status,ierror)
			call MPI_RECV(rms2,1,MPI_REAL,other,4,MPI_COMM_WORLD,status,ierror)
			dad = dad + dad2
			dac = dac + dac2
			rms = rms + rms2
			dad = SQRT(dad/19.)
			dac = SQRT(dac/19.)
			rms = SQRT(rms/(19.*19.))
			write(10+mypid,'(xa,3e15.5)'),'rms,dad,dac=',rms,dad,dac
		endif

		b(0:20,N-1:M+1)=a(0:20,N-1:M+1)
	enddo

	do j = N-1, M+1
		write(10+mypid,'(xa,i2,i3,1x,21f6.3)')   'id',mypid,j,(a(i,j),i=0,20)
	enddo

	do j=s+1, 20
		if (mypid.ne.0) then
			call MPI_SEND(a(0,j),20,MPI_REAL,other,tag,MPI_COMM_WORLD,ierror)
		else
			call MPI_RECV(a(0,j),20,MPI_REAL,other,tag,MPI_COMM_WORLD,status,ierror)
		endif
	enddo
	call MPI_BARRIER(MPI_COMM_WORLD, ierror)

	if (mypid.eq.0) then
		write(10+mypid,'(21f7.4)') ((a(i,j),i=0,20),j=0,20)
	endif

	call MPI_FINALIZE(ierror)

	end program AS73f




