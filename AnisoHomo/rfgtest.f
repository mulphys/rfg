*     This is a test of RFG routine
*     Author: A.Smirnov, www.cemr.wvu.edu/~andrei
*
*     run the compiled file as
*
*     ./rfgtest > rfgtest.dat
*
*     The produced rfgtest.dat file contains the 
*     random velocity data in the format defined
*     by the write statement below.
*
*     To visualize the data you may split the file
*     into three (on Unix):
*
*	tr -s ' ' '\t' < rfgtest.dat | cut -f 2-3,4 > U.dat
*	tr -s ' ' '\t' < rfgtest.dat | cut -f 2-3,5 > V.dat
*	tr -s ' ' '\t' < rfgtest.dat | cut -f 2-3,6 > W.dat
*
*     corresponding to each velocity component.
*     To display each component as a function of 
*     time and space run the gnuplot splot:
*
*	gnuplot> splot 'U.dat' w l
*
*     you may use set surf and set hidden3d to 
*     get a better view
*
      program RFGTEST
      DOUBLE PRECISION 
     & x(3),v(3) ! space coordinates and velocity vectors
     & ,turbtime,turblength(3),t,dt,dx
      nspec=1000     ! spectral sample size
      turbtime=1.0   ! turbulence time-scale
      turblength(1)=1.0 ! turbulence length-scale
      turblength(2)=1.0 ! turbulence length-scale
      turblength(3)=1.0 ! turbulence length-scale
*
*     Initialize the spectrum
*
      call genspec(nspec,turbtime,turblength)
*
*     the call above will allocate memory for
*     10*nspec double presision words and 
*     set up the turbulence spectrum to be
*     used in subsequent calls to genvec-function
*
*     Run the time-space loop
*
      nt=10  ! number of time steps
      dt=1.  ! time-step
      nx=100 ! number space-points
      dx=1.  ! their separation
      do it=0,nt-1
      do ix=0,nx-1
*        Only one space dimension is used 
*        in this example
         t=it*dt
         x(1)=0.0
         x(2)=0.0
         x(3)=ix*dx
         call genvec(t,x,v)
*
*        the call above returns a random vector v(1:3)
*        for every time t and space location x(1:3)
*        These vectors form an isotropic homogeneous 
*        divergence-free random flow field with the 
*        avarage <v(i)>=0.0 and the cross-correlations 
*        <v(i)*v(i)> -> 1.0 
*        <v(i)*v(j)> -> 0.0, i.ne.j
*        where <.> is an average over a time interval
*        T>>dt, or/and over a space domain of size
*        D>>sqrt(dx**2+dy**2+dz**2)
*
         write (*,*) t,x(3),v(1),v(2),v(3)
      enddo
      write(*,*)
      enddo
      call delspec
      end
