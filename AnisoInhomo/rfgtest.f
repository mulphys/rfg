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
     & ,turbtime,turblength(3),uu(6)
     & ,t,dt,dx
      nspec=1000     ! spectral sample size
*
*     Initialize the spectrum
*
      call genspec(nspec)
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
*        Set time:
         t=it*dt
*        Set space coordinates:
         x(1)=0.0
         x(2)=0.0
         x(3)=ix*dx
*        Turbulent time scale:
         turbtime=1.0   ! turbulence time-scale
*        Components of Reynolds stresses:
         uu(1)=1.    !UU: <1,1>
         uu(2)=0.    !UV: <1,2>
         uu(3)=1.    !VV: <2,2>
         uu(4)=0.    !UW: <1,3>
         uu(5)=0.    !VW: <2,3>
         uu(6)=1.    !WW: <3,3>
*        Call the routine to generate the random vector
         call genvec(t,x,turbtime,uu,v)
*        The call above returns a random vector v(1:3)
*        for every time t and space location x(1:3),
*        and for given turbulent time and length-scales,
*        and Reynolds stress tensor components.
*        For a homogeneous case these vectors form a
*        divergence-free random flow field with the 
*        avarage <v(i)>=0.0 and the cross-correlations 
*        <v(i)*v(j)> -> uu(k)
*        where <.> is an average over a time interval
*        T>>dt, or/and over a space domain of size
*        D>>sqrt(dx**2+dy**2+dz**2),
*        and k=k(i,j) as given by the correspondence
*        between uu(k) and <i,j> above.
*        Output the returned vectors:
         write (*,*) t,x(3),v(1),v(2),v(3)
      enddo
      write(*,*)
      enddo
      call delspec
      end
