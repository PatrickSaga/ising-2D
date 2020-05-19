	program ising2d

	! COLLECTIVE PHENOMENA PROJECT - 2D ISING MODEL
	! Patrick Sánchez Galea
	! MC algorithm using the Metropolis method for the 2D Ising model.


	implicit none

	integer :: i, j, m, n		! dummy integers
	integer, allocatable :: A(:,:)	! spins' matrix
	integer :: L				! length of the lattice, # of spins
	integer :: newspin			! values of changed spin
	integer :: nscan, iscan		! number/current scans at different T
	integer :: nstep, istep		! MC-time and current pass number
	integer :: ntherm			! number of thermalization steps
	integer :: nvalues			! # of points after thermalization
	real :: T, beta			! T and 1/KT	
	real :: minT, maxT, dT		! limit temperatures to scan, and step
	real :: dU				! change of energy between configs.
	real :: mag, mag_avg, mag2_avg	! magnetization, cumulative avg m/m2
	real :: energy, energy_avg, energy2_avg	! same than mag for energy
	double precision u, dran_u, logran	! random number generator
	

	! read the input parameters from "parameters"
	open(unit=11, file='parameters', status='old', action='read')
	read(11,*); read(11,*) L
	read(11,*); read(11,*) nstep
	read(11,*); read(11,*) ntherm
	read(11,*); read(11,*) maxT
	read(11,*); read(11,*) minT
	read(11,*); read(11,*) dT
	close(11)


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
	do i = 1, L +2	! begin boundary conditions
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
		energy_avg = 0.d0
		energy2_avg = 0.d0
		

		! Core of the MC algorithm
		MC_loop: do istep = 0, nstep

			!Choose a random spin to swap and compute dU
			m = nint((L-1)*dran_u() + 2)
			n = nint((L-1)*dran_u() + 2)
			newspin = -A(m,n)
			dU=-2*newspin*(A(m-1,n)+A(m+1,n)+A(m,n-1)+A(m,n+1))

			! Metropolis step 
			! If dU<0, swap
			! If dU>0, swap when exp(-beta*dU) > random(0,1)
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
			    mag_avg = mag_avg + mag
			    mag2_avg = mag2_avg + mag**2
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
			  endif
			endif
			
		enddo MC_loop
		
		! write the output files
		write(21,*) T,',',abs(mag_avg/nvalues),',',
     *        mag2_avg/nvalues, 
     *		',',beta*(mag2_avg/nvalues - (mag_avg/nvalues)**2)
		write(22,*) T,',',energy_avg/nvalues,',',
     *        energy2_avg/nvalues,',',(beta**2)*(energy2_avg/nvalues
     *        -(energy_avg/nvalues)**2)

	enddo temp_loop
	
	close(21)
	close(22)
	close(23)
	print*, "Program Complete"

	end program ising2d