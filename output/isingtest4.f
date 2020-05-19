	program isingtest3

	! COLLECTIVE PHENOMENA PROJECT - 2D ISING MODEL
	! Patrick Sánchez Galea
	! MC algorithm using the Metropolis and Wolff methods for the 2D Ising model.


	implicit none

	integer :: i, j, m, n					! dummy integers
	integer, allocatable :: A(:,:)				! spins' matrix
	integer :: L							! length of the lattice, # of spins
	integer :: newspin						! values of changed spin
	integer :: nscan, iscan					! number/current scans at different T
	integer :: nstep, istep					! MC-time and current pass number
	integer :: ntherm						! number of thermalization steps
	integer :: nvalues						! # of points after thermalization
	real :: T, beta						! T and 1/KT	
	real :: minT, maxT, dT					! limit temperatures to scan, and step
	real :: dU							! change of energy between configs.
	real :: mag, mag_avg, mag2_avg, mag4_avg			! magnetization, cumulative avg m/m2
	real :: energy, energy_avg, energy2_avg, energy4_avg	! same than mag for energy
	real :: pa							! probability for the Wolff algorithm
	double precision u, dran_u, logran				! random number generator
	
	

	! read the input parameters from "parameters"
	open(unit=11, file='parameters', status='old', action='read')
	read(11,*); read(11,*) L
	read(11,*); read(11,*) nstep
	read(11,*); read(11,*) ntherm
	read(11,*); read(11,*) maxT
	read(11,*); read(11,*) minT
	read(11,*); read(11,*) dT
	close(11)

	ntherm = 10000
	nstep = 10000*L**2 + ntherm


	! open the output files
	open(unit=21,file='magnetization',status='replace',action='write')
	write(21,*) "temp      avg_mag      avg_mag2      susceptibility"
	open(unit=22, file='energy', status='replace', action='write')
	write(22,*) "temp      avg_energy      avg_energy2      Cv"

	
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
		mag2_avg = 0.d0
		mag4_avg = 0.d0
		energy_avg = 0.d0
		energy2_avg = 0.d0
		pa = 1-dexp(-2.d0*beta)


		! Core of the MC algorithm

		if (T .le. 2.6d0 .and. T .ge. 2.1d0) then		! Wolff near Tc
			MC_loop1: do istep = 0, nstep/100

				!Choose a random spin to swap and clusterize
				m = nint((L-1)*dran_u() + 2)
				n = nint((L-1)*dran_u() + 2)
				call clustering(A, m, n, L, pa)

				! When the thermalization time is reached, start measuring
				! 10 cluster-flips to measure, 1000 values
				if (istep > nstep/100-10000) then
				  if (mod(istep,10) == 0) then  
				    nvalues = nvalues + 1	! counter of steps to average
				    mag = sum(A(2:L+1,2:L+1))/(1.d0*L**2)
				    mag_avg = mag_avg + abs(mag)
				    mag2_avg = mag2_avg + mag**2
				    mag4_avg = mag2_avg + mag**4
				    energy = 0.d0
				    do i = 2, L+1
				      do j = 2, L+1
					 energy = energy - A(m,n)*(A(m-1,n)+A(m+1,n)+ 
     * 			  		     A(m,n-1)+A(m,n+1))
				      enddo
				    enddo
				    energy = energy/(2*L**2)
				    energy_avg = energy_avg + energy
				    energy2_avg = energy2_avg + energy**2
				    energy4_avg = energy4_avg + energy**4
				  endif
				endif
			enddo MC_loop1
		else							! Metropolis if T is far from Tc
			! Core of the MC algorithm
			MC_loop2: do istep = 0, nstep

				!Choose a random spin to swap
				m = nint((L-1)*dran_u() + 2)
				n = nint((L-1)*dran_u() + 2)

				! Metropolis 
				! If dU<0, swap
				! If dU>0, swap when exp(-beta*dU) > random(0,1)
				newspin = -A(m,n)
				dU=-2*newspin*(A(m-1,n)+A(m+1,n)+A(m,n-1)+A(m,n+1))
				if (dU < 0.d0) then
				  A(m,n) = newspin
				  if (m == 2) A(L+2,n) = newspin
				  if (m == L+1) A(1,n) = newspin
				  if (n == 2) A(m,L+2) = newspin
				  if (n == L+1) A(m,1) = newspin
				else
				  logran = dlog(dran_u())	
				  if (-beta*dU > logran) then
				    A(m,n) = newspin
				    if (m == 2) A(L+2,n) = newspin
				    if (m == L+1) A(1,n) = newspin
				    if (n == 2) A(m,L+2) = newspin
				    if (n == L+1) A(m,1) = newspin
				  endif
				endif

				! When the thermalization time is reached, start measuring
				! MC iteration when each spin is flipped once in average
				if (istep > ntherm) then
				  if (mod(istep,L**2) == 0) then  
				    nvalues = nvalues + 1	! counter of steps to average
				    mag = sum(A(2:L+1,2:L+1))/(1.d0*L**2)
				    mag_avg = mag_avg + abs(mag)
				    mag2_avg = mag2_avg + mag**2
				    mag4_avg = mag2_avg + mag**4
				    energy = 0.d0
				    do i = 2, L+1
				      do j = 2, L+1
					 energy = energy - A(m,n)*(A(m-1,n)+A(m+1,n)+ 
     * 			  		     A(m,n-1)+A(m,n+1))
				      enddo
				    enddo
				    energy = energy/(2*L**2)
				    energy_avg = energy_avg + energy
				    energy2_avg = energy2_avg + energy**2
				    energy4_avg = energy4_avg + energy**4
				  endif
				endif
			enddo MC_loop2
		endif 
		! write the output files
		write(21,*) T,',',mag_avg/nvalues,',',
     *        mag2_avg/nvalues, 
     *		',',beta*(mag2_avg/nvalues - (mag_avg/nvalues)**2),
     *		',',mag4_avg/nvalues
		write(22,*) T,',',energy_avg/nvalues,',',
     *        energy2_avg/nvalues,',' , (beta**2)*(energy2_avg/nvalues
     *        -(energy_avg/nvalues)**2), ',', energy4_avg/nvalues
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
		A(m,n) = -A(m,n)					! flip the symmetric spin
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
