	program AS48f

	implicit none

	real*8, dimension(0:20,0:20) :: a, b, c
	real*8, dimension(19) :: aed, aec, d, v, gm, bt
	real*8 rms,dad,dac/0./
	real*8 a1,b1,c1,d1,R,pi
	integer i, j, k
	logical mask(0:20,0:20)

	aed = (/.4972559,.4890261,.4753282,.4562242,.4318678,.4025638,.3688278,.3314301,.2914033,.2499993,.2085953,.1685686,.1311710, & 
		.0974351,.0681313,.0437752,.0246714,.0109737,.0027441/)
	aec = (/.4972645,.4890595,.4753997,.4563411,.4320283,.4027535,.3690196,.3315883,.2914935,.2500000,.2085065,.1684117,.1309804, & 
		.0972465,.0679717,.0436589,.0246003,.0109505,.00273552/)
	pi = 4.d0*atan(1.d0)

	write(*,'(a)')' AS48/P1 33-241, F08 GREER'

	mask=.false. ; do i=1,19 ; mask(i,i)=.true. ; enddo
	v=0.d0
	d=0.d0
	b=0.d0 ; b(1:19,0)=1.d0
	a=0.d0 ; a(1:19,0)=1.d0 ! start clever init
	do i=1,19
		do j=1,19
			a1 = (20.d0-i)/20.d0
			b1 = (i)/20.d0
			c1 = (j)/20.d0
			d1 = (20.d0-j)/20.d0
			a(i,j) = (a(20,j)/a1 + a(0,j)/b1 + a(i,0)/c1 + a(i,20)/d1)/(1.d0/(a1*b1) + 1.d0/(c1*d1))
		enddo
	enddo
	write(*,'(x21i4)') int(1000*a+.5d0)

	dad = sqrt(sum((pack(a,mask)-aed)**2)/19.d0)
	dac = sqrt(sum((pack(a,mask)-aec)**2)/19.d0)

	write(*,'(a,2e15.5)'),' AS48/P2 dad,dac=',dad,dac ! everything good to here, i think

	do k=1,10
		R = 4.d0*(sin((pi*k)/(2.d0*20.d0)))**2.d0
        	c=a
! do the columns going from a to b ...
		do i=1,19 ! columns.
			!{... compute of constant d using a, mostly array notation ...}

	! unsure of how to calculate d() array

			d = -a(i+1,1:19) - (R-2.d0)*a(i,1:19) - a(i-1,1:19)
			d(1) = d(1) - a(i,0)
			d(19) = d(19) - a(i,20)


			call tridi(-2.d0-R)
			! {. . . place new column in b with array notation . . .}
			b(i,1:19) = v
		enddo
! do the rows going from b to a ...
		do j=1,19 ! rows.
			!{... compute of constants d using b, mostly array notation ...}


			d = -b(1:19,j+1) - (R-2.d0)*b(1:19,j) - b(1:19,j-1)
			d(1) = d(1) - b(0,j)
			d(19) = d(19) - b(20,j)


			call tridi(-2.d0-R)
			!{. . . place new row in a with array notation . . .}
			a(1:19,j) = v
		enddo
		dad = sqrt(sum((pack(a,mask)-aed)**2)/19.d0)
		dac = sqrt(sum((pack(a,mask)-aec)**2)/19.d0)
		rms = sqrt(sum((a-c)**2))/19.d0
		write(*,'(a,i4,3e15.5)'),' AS48/P3 k,rms,dad,dac=',k,rms,dad,dac
	enddo

	write(*,'(x21i4)') int(1000*a+.5d0)

	contains
		subroutine tridi(b_val) ! fairly certain all of this is correct
		! local declarations
		real*8 b_val
		integer q

		bt=0.d0
		gm=0.d0
		v=0.d0

		bt(1) = b_val
		gm(1)=d(1)/bt(1)

		do q=2,19
			! forward recursion relation using d() to get dm() & bt()
			bt(q)=b_val-1.d0/bt(q-1)
			gm(q)=(d(q)-gm(q-1))/bt(q)
		enddo

		v(19) = gm(19)

		do q=18,1,-1
			! backward recursion relation using gm() & bt() to get v()
			v(q)=gm(q)-v(q+1)/bt(q)
		enddo

		return
		end subroutine tridi

	end program AS48f




