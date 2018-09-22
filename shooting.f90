 Program shooting
!---------------------------------------------------------------------
! Solver for 1D boundary-value problem second-order ODE 
! Method: unilizes the shooting method based on the method of secants
! (calls 4th-order Runge-Kutta to solve the initial value problem)
! Version for d2y/dx2 = f(x,y,dy/dx) equation
! COMPILE: gfortran -o shooting shooting.f90
! RUN: ./shooting
! GNUPLOT SCRIPT FOR PLOTTING: gnuplot 'plot.txt'
!----------------------------------------------------------------------

implicit none

	integer		 			:: i,m,n,t,s  ! i=counter,n=number of base points
	!m,t,s used for file operations
	real*8, external 			:: f	      ! external function
	real*8, allocatable, dimension(:)	:: x, y, dy   ! x, y, y'
	write(*,*)"Enter the value of number of base points:"
	read(*,*)n
	allocate(x(n),y(n),dy(n))

	OPEN(20,ACCESS='SEQUENTIAL',STATUS='REPLACE',FILE='shooting.dat',IOSTAT=m)

	! boundary values
	x(1) =  1.0
	y(1) =  1.0
	x(n) =  1.5
	y(n) =  1.2

	! assumptions for y'(1) - use dy(1) and dy(2) here only as a storage
	dy(1) =  2.0
	dy(2) =  2.5

	call shoot(f,x,y,dy,n)

	write(20,100)
	do i=1,n
  		write (20,101) x(i), y(i), dy(i)
	end do

100     format(5x,'x',11x,'y',11x,'dy')
101	format(3(1pe12.4))
! check file integrity
t=0;
do 
	read(20,*,iostat=m)
	if(m<0) exit
	t=t+1
end do

	rewind(20)
	close (20) !Close file

!Open gnuplot command file
OPEN(30,ACCESS='SEQUENTIAL',STATUS='REPLACE',FILE='plot.txt',IOSTAT=s,position='rewind')
write(30,*)"set terminal png"
write(30,*)"set output 'solution.png' "
write(30,*)"set title 'SOLUTION OF SECOND ORDER DE USING SHOOTING METHOD' "
write(30,*)"set xlabel 'xval' "
write(30,*)"set ylabel 'yval'  "
write(30,*)"plot  'shooting.dat' using 1:2 w lp  "

	rewind(30)
	close(30) !Close file

deallocate(x,y,dy)

	stop

end program shooting

  subroutine shoot(f,x,y,dy,n)
!----------------------------------------------------------------------
! Solution of the boundary-value second-order 1D ODE
! d2y/dx2 = f(x,y,dy/dx) with Dirichlet boundary conditions
! y(xmin) = ..., and y(xmax) = ...
! Method: unilizes the shooting method based on the method of secants
! (calls 4th-order Runge-Kutta to solve the initial value problem)
!----------------------------------------------------------------------
! input ...
!  f(x,y,dy)  - function d2y/dx2 (supplied by a user)
!  x(1), x(n) - boundary points
!  y(1), y(n) - boundary values (Dirichlet boundary conditions)
! dy(1),dy(2) - two guesses for y'(x(1))
!
! output ...
!  y(i) and dy(i) solutions at points x(i) (i=1,...,n)
! note: dy corresponds to y' (the first derivative)
!----------------------------------------------------------------------
implicit none
	real*8, parameter  			:: eps =1e-12                ! target tolerance
	integer,intent(in)			:: n
	real*8, external       			:: f 
	real*8, dimension(1:n),intent(inout)	:: x, y, dy
	integer 				:: i,j
	real*8 		       			:: dx, yn
	real*8, allocatable, dimension(:)	:: g,c
	allocate(g(n),c(n))

! g(n) - array of dy(n) values [iterative improvement]
! c(n) - array of y(n) values [iterative solutions for g(it)]
! first guesses for g(n)
	g(1) = dy(1)
	g(2) = dy(2)

! remember the second boundary condition 
! since y(n) recalculated for each new g(i)
	yn = y(n)

! generate base points x(i) from x(1), x(n) and n
	dx = (x(n)-x(1))/float(n-1)
	do i=2,n
  		x(i) = x(i-1)+dx
	end do

! shooting iterations (for the first two - we use assumed values of dy(1))
	do j=1,n
  		dy(1) = g(j)
  		call rk4_2d(f,x,y,dy,n) !Solves initial value ODE on n-base points
  		c(j) = y(n)
  		if(abs(yn-c(j)) <= eps) exit	!Check convergence
  		if(j.ge.2) g(j+1)=g(j)-(c(j)-yn)*(g(j)-g(j-1))/(c(j)-c(j-1)) ! secant method
	end do
end subroutine shoot

  function f(x,y,dy)
!------------------------------------------
! the second derivative - use original ODE
! d2y/dx2 = f(x,y,dy)
!------------------------------------------
  implicit none
  	real*8 :: f
	real*8,intent(inout):: x, y, dy
    	f = exp(-x)**3
  end function f

  subroutine rk4_2d(f,x,y,dy,n)
!-----------------------------------------------------------------
! Solution of the second-order 1D ODE IVP
! method:  Runge-Kutta 4th-order
!-----------------------------------------------------------------
! input ...
!  f(x,y,dy)- function from d2y/dx2=f(x,y,dy) (supplied by a user)
!  x(i)  - array of x base point (n points)
!  y(1)  - initial value for y(1)
! dy(1) - initial value for y'(1)
!
! output ...
!  y(i)  - solutions for y  in n points
! dy(i)  - solutions for y' in n points 
!-----------------------------------------------------------------
	implicit none
	real*8,external				:: f
	integer,intent(in)			:: n
	real*8, dimension(1:n),intent(inout) 	:: x, y, dy
	integer 				:: i
	real*8 					:: h,k11,k12,k21,k22,k31,k32,k41,k42

	do i=2,n
	   h   = x(i)-x(i-1)
	   k11 = h*dy(i-1)
	   k12 = h*f(x(i-1),y(i-1),dy(i-1))
	   k21 = h*(dy(i-1)+k12/2.0)
	   k22 = h*f(x(i-1)+h/2.0, y(i-1)+k11/2.0, dy(i-1)+k12/2.0)
	   k31 = h*(dy(i-1)+k22/2.0)
	   k32 = h*f(x(i-1)+h/2.0, y(i-1)+k21/2.0, dy(i-1)+k22/2.0)
	   k41 = h*(dy(i-1)+k32)
	   k42 = h*f(x(i-1)+h,y(i-1)+k31,dy(i-1)+k32)

	   y(i) = y(i-1) + (k11 + 2.0*(k21+k31) + k41)/6.0
	  dy(i) = dy(i-1)+ (k12 + 2.0*(k22+k32) + k42)/6.0
	end do 
end subroutine rk4_2d
