	program isingtest3

	! COLLECTIVE PHENOMENA PROJECT - 2D ISING MODEL
	! Patrick Sánchez Galea
	! MC algorithm using the Metropolis and Wolff methods for the 2D Ising model.


	implicit none

	integer :: i, j, m, n, d, k					! dummy integers
	integer, allocatable :: A(:,:)				! spins' matrix
	integer :: L							! length of the lattice, # of spins
	integer :: newspin						! values of changed spin
	integer :: nscan, iscan					! number/current scans at different T
	integer :: nstep, istep					! MC-time and current pass number
	integer :: ntherm						! number of thermalization steps
	integer :: nvalues						! # of points after thermalization
	real :: T, beta						! T and 1/KT	
	real :: minT, maxT, dT					! limit temperatures to scan, and step
	real :: dU							! change of energy between configs
	real :: pa							! probability for the Wolff algorithm
	double precision :: mag, mag_avg, mag2_avg, mag4_avg			! magnetization, cumulative avg m/m2
	double precision :: energy, energy_avg, energy2_avg, energy4_avg	! same than mag for energy
	double precision :: chi, err_chi, sh, err_sh				! susceptibility and heat capacity
	double precision u, dran_u, logran, corr(127)					! random number generator
	
	

	! read the input parameters from "parameters"
	open(unit=11, file='parameters', status='old', action='read')
	read(11,*); read(11,*) L
	read(11,*); read(11,*) nstep
	read(11,*); read(11,*) ntherm
	read(11,*); read(11,*) maxT
	read(11,*); read(11,*) minT
	read(11,*); read(11,*) dT
	close(11)
	nstep = nstep*L**2 + ntherm


	! open the output files
	open(unit=21,file='correlation',status='replace',action='write')

	
	! we use periodic bc, so we have the same first and last rows and columns
	allocate(A(L+2,L+2))
	nscan = int((maxT-minT)/dT) + 1
	call dran_ini(481992)


	! Initial configuration of the spins. Chosen to be random
	do i = 2, L +1 
		do j =2, L +1
			if (dran_u() .le. 0.5d0) then 
				A(i,j) = int(-1)
			else 
				A(i,j) = int(1)
			endif
		enddo
	enddo
	do i = 1, L+2		! begin boundary conditions
		A(i,1) = A(i,L+1)
		A(i,L+2) = A(i,2) 
	enddo
	do j = 1, L+2
		A(1,j) = A(L+1,j)
		A(L+2,j) = A(2,j)
	enddo			! end boundary conditions



	! Scan for the different temperatures
	temp_loop: do iscan = 1, nscan
		T = maxT-dT*(iscan-1)
		print*, "Program running for T = ", T
		
		! Initialization of the variables
		beta = 1.d0/T
		nvalues = 0
		mag_avg = 0.d0
		pa = 1-dexp(-2.d0*beta)

		! Core of the MC algorithm
		if (T .le. 2.31d0 .and. T .ge. 2.29d0) then		! Wolff near Tc
			MC_loop1: do istep = 0, 60000
				m = nint((L-1)*dran_u() + 2)
				n = nint((L-1)*dran_u() + 2)
				call clustering(A, m, n, L, pa)

				if (istep > 10000) then
				  if (mod(istep,100) == 0) then  

				    nvalues = nvalues + 1	! counter of steps to average
				    mag = sum(A(2:L+1,2:L+1))/(1.d0*L**2)
				    mag_avg = mag_avg + abs(mag)
				    
				    do d=1, L/2-1
				     
				    enddo

				  endif
				endif

			enddo MC_loop1
			do d=1, L/2-1
			 corr(d)=corr(d)/(4*L**2*nvalues)
			 corr(d)=corr(d)-mag_avg/nvalues
			 write(21,*) corr(d), d
			enddo
			print*, mag_avg/nvalues
		endif
	enddo temp_loop
	
	close(21)
	close(22)
	close(23)
	print*, "Program Complete"


	contains
	! Subroutine to form clusters : Wolff algorithm
	recursive subroutine clustering(A, m, n, L, pa)
	integer L, m, n, k, m2, n2
	real :: pa
	integer :: A(L+2,L+2)

	A(m,n) = -A(m,n)	! flip the spin
	if (m==1 .or. m == L+2 .or. n == 1 .or. n == L+2) then	! boundary
		if (m == 1) m = L+1
		if (m == L+2) m = 2
		if (n == 1) n = L+1
		if (n == L+2) n = 2
		A(m,n) = -A(m,n)	! flip the symmetric spin
	endif
	do i = 1, L+2		! begin boundary conditions
		A(i,1) = A(i,L+1)
		A(i,L+2) = A(i,2) 
	enddo
	do j = 1, L+2
		A(1,j) = A(L+1,j)
		A(L+2,j) = A(2,j)
	enddo			! end boundary conditions

	if (A(m-1,n) == -A(m,n)) then
		if (dran_u() .le. pa) call clustering(A,m-1,n,L,pa)
	endif
	if (A(m+1,n) == -A(m,n)) then
		if (dran_u() .le. pa) call clustering(A,m+1,n,L,pa)
	endif
	if (A(m,n-1) == -A(m,n)) then
		if (dran_u() .le. pa) call clustering(A,m,n-1,L,pa)
	endif
	if (A(m,n+1) == -A(m,n)) then
		if (dran_u() .le. pa) call clustering(A,m,n+1,L,pa)
	endif
	end subroutine clustering


	end program isingtest3
