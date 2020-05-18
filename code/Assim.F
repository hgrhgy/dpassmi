*     ==================================================================
*     program DPenAssim
*     Program to do (1+1)dVar with the double pendulum system
*     Modified double pendulum (system 2)
*     Ross Bannister
 
*     History
*     -------
*     12.09.01 Original code.
*     15.11.01 Include output of init conditions (each cycle) with errs,
*              J, Jb, Jo for each iter.
*     18.12.01 Create new section for simple insertion of observations.

implicit none
*     Constants
*     ---------
      integer maxobs,winsize
      parameter (maxobs=5000, winsize=500)
*     maxobs: maximum # of observations
*     winsize: maximum # of timesteps in an assim cycle

*     Integer variables
*     -----------------
      integer Nt,Nout,Nobs,Nwindow,ObTypes(maxobs),Cycles,AssCy,Nsteps
      integer ObFirst,ObLast,Index(maxobs),remaining,loop,safety,t
*     Nt: total # time steps
*     Nout: Frequency of obs output
*     Nobs: # observations in total
*     Nwindow: # time steps in typical assimilation cycle
*     ObTypes: 'degree-of-freedom' index for ob
*     Cycles: # of cycles for whole run
*     Nsteps: # time steps in current assimilation cycle
*     ObFirst: index to first ob in current cycle
*     ObLast: index to last ob in current cycle
*     Index: links observations to model states (like pointer)
*     safety: max number of Var. iterations
*     t: time loop (used when 'inserting' observations)

*     Real variables
*     --------------
      real x0(4),g,dt,pi,l1,l2,l3,m1,m2,m3,noisepos,noisevel,t1,t2
      real ActualObs(maxobs),ObTimes(maxobs),ObErrs(maxobs),Binv(4,4)
      real x(winsize,4),xB(4),xANp(4),xLS(4),xBp(4),GradJOb(maxobs)
      real U1(4,4),U2(4,4),GradJTrans(4),ms,crit,H(maxobs,4),Hess(4,4)
      real xfore(winsize,4),x0f(4),J,Jb,Jo,HessInv(4,4)
*     x0: Initial conditions for whole run
*     g,l1,l2,l3,m1,m2,m3: model parameters
*     noisepos, noisevel: random noise for angles and angular velocities
*     ActualObs: the observations
*     ObTimes: the times of the observations (real time)
*     ObErrs: the uncertainties in the observtions
*     x: for model integration
*     xb: background state of current cycle
*     xANp: current version of the analysis increment xANp = xAN - xLS
*     xLS: current linearization state
*     xBp: background increment (xBp=xB-xLS)
*     GradJOb: gradient of J wrt model obs
*     U2: eigenvecs of transformation (rows of trans(U2) are eigenvecs)
*     U1: eigenvalues of transformation (diagonal elements)
*     GradJTrans: gradient of J wrt transformed state vector
*     crit: cut-off for convergence of Var.
*     H: d(modelob(t))/dx(0)
*     Hess: Hessian matrix (model space)
*     HessInv: Inverse of Hessian
*     Binv: Inverse of background error co-variance matrix
*     t1,t2: Times bounding the assimilation cycle
*     xfore: Forecast using originial initial conditions (for info)
*     x0f: Initial conditions for current cycle from originial initial
*          conditions (with no assimilation)
*     J, Jb, Jo: cost function and components

*     Characters
*     ----------
      character*7 Status
      character*8 fname

*     Boolean
*     -------
      logical finish

*     Workspace arrays
*     ----------------
      real Work1(maxobs),Work2(maxobs),Work3(winsize,4)

*     Functions
*     ---------
      real modsq,mymod

*     ==================================================================
      call Init (pi,crit,safety)

*     Read-in the input data
      call ReadConfig (Status,x0,g,dt,Nt,Nout,l1,l2,l3,m1,m2,m3,pi,
     :                 noisepos,noisevel,ActualObs,ObTypes,ObTimes,
     :                 ObErrs,Nobs,Nwindow,maxobs)

      if (Status.eq.'MakeObs') then

*       ----------------------------------------------------------------
*       || *** SECTION I: MAKE OBSERVATIONS ***                       ||
*       ----------------------------------------------------------------
*       This run is to make-up some observations
*       ----------------------------------------
        print*,'MAKING OBSERVATIONS'
        print*,'-------------------'
        call MakeObs (x0,g,dt,pi,Nt,Nout,l1,l2,l3,m1,m2,m3,noisepos,
     :                noisevel)


      elseif (Status.eq.'Analyse') then


*       ----------------------------------------------------------------
*       || *** SECTION II - ANALYSE OBSERVATIONS ***                  ||
*       ----------------------------------------------------------------
*       This run is to analyse some observations
*       ----------------------------------------
        print*,'ANALYSE OBSERVATIONS'
        print*,'--------------------'

*       Open the analysis integration output file
        open (13,file='Assim.anal')
*       Open files containing observations used
        open (14,file='Assim.obs1')
        open (15,file='Assim.obs2')
        open (16,file='Assim.obs3')
        open (17,file='Assim.obs4')
*       Open file containing cycle boundary times
        open (18,file='Assim.bnds')
*       Open file for complete forecast from originial initial conds
        open (19,file='Assim.fore')
*       Unit 20 is to output the test integration after each iteration
*       Open file for the analysis initial conditions (and errors)
        open (21,file='Assim.init')


*       Calculate the number of assimilation cycles
        Cycles=int(real(Nt)/real(Nwindow))
        if (Cycles*Nwindow.lt.Nt) Cycles=Cycles+1
*       The number of time steps in an assimilation cycle is usually
        Nsteps=Nwindow
*       Point to end of observations (initial value)
        ObLast=0

        do AssCy=1,Cycles
          print*,'-----------------------------------------------------'
          print*,'Assimilation Cycle No. ',AssCy
          print*,'-----------------------------------------------------'

          ObFirst=ObLast+1
*         Check to see if there are < Nwindow left in the time domain
          remaining=Nt-(AssCy-1)*Nwindow
          if (remaining.lt.Nwindow) Nsteps=remaining

*         Work out times bounding this assimilation cycle
          t1=real((AssCy-1)*Nwindow)*dt
          t2=real((AssCy-1)*Nwindow+Nsteps)*dt

*         Output the boundary time
          call BoundTime (t1,18)

*         Find the last observation pertaining to current cycle
*         (and set-up the observation indices)
          call FindLast (ActualObs,ObTimes,ObErrs,ObTypes,ObFirst,
     :                   ObLast,Index,Nobs,t1,t2,dt,14,maxobs)
          print*,'Start and end times : ',t1,t2
          print*,'First ob index,time : ',ObFirst,ObTimes(ObFirst)
          print*,'Last ob index,time  : ',ObLast,ObTimes(ObLast)

*         Store the background for this cycle from either given initial
*         conditions, or from previous forecast
          call Background (x0,x,AssCy,xB,Nwindow,winsize)

*         Set the linearization state (l.s.) to the background
          call Equal (xLS,xB,4,4)
          print*,'LS ',xLS(1),xLS(2),xLS(3),xLS(4)

          if (ObLast.ge.ObFirst) then
*           There are obs in this window - do data assimilation
*           ---------------------------------------------------
            print*,'Analysing ',ObLast-ObFirst+1,' observations'

            loop=0
10          loop=loop+1
            print*,'  -----------------------------------'
            print*,'  Loop number ',loop

*           Find the background as an increment to the l.s. (xBp=xB-xLS)
            call PlusMinus (xBp,xB,xLS,-1.,4,4)

*           Do a forward integration, starting from xLS
            call Advance (x,xLS,Nsteps,1,0.,m1,m2,m3,l1,l2,l3,g,dt,pi,
     :                    winsize)

*           Output the current guess at state
            call SetUpName (AssCy,loop,fname)
            open (20,file=fname)
            call OutputIntegration (x,Nsteps,Nwindow,AssCy,20,dt,
     :                              winsize)
            close (20)

            if (loop.eq.1) then
*             If first run, calculate Binv based on this integration
              if (AssCy.eq.1) call CalcBinv (x,Binv,Nsteps,pi,winsize)
*             Start/continue forecast (from original initial conds)
            endif

*           Calculate the the gradient of J w.r.t. model observations
*           (this is E^(-1) (Obs - Modelobs)
*           Also convenient to calculate Jo also (for diagnosis)
            call GradientJOb (x,ActualObs,ObTypes,ObErrs,Index,GradJOb,
     :                        Jo,ObFirst,ObLast,pi,winsize,maxobs)

*           Calculate Jb
            call CalcJb (xB,x,Jb,Binv,winsize)
*           The total cost function
            J=Jb+Jo
            print 113,'  J,Jb,Jo ',J,Jb,Jo

*           Compute the forward model matrix, H
            call ForwardModel (xLS,Nsteps,m1,m2,m3,l1,l2,l3,g,dt,pi,
     :                         ObFirst,ObLast,ObTypes,Index,Work1,Work2,
     :                         Work3,H,winsize,maxobs)

*           Compute the Hessian in model space
            call Hessian (Binv,ObErrs,H,Hess,Obfirst,Oblast,maxobs)

*           Diagonalize the Hessian
            call Jacobi (Hess,U1,U2,4,4,0.00001)
*           Check eigenvectors
*            call CheckDiag (Hess,U2)

*           Check signs of eigenvalues
            call CheckPos (U1,4)

*           Currently: Diagonal elements of U1 are eigenvalues
*                      Rows of U2 are eigenvectors
*           Want     : Diagonal elements of U1 to be 1/sqrt(e-values)
*                      Columns of U2 to be e-vectors
            call InvSqrt (U1,4)
            call Transpose (U2,4)
*           Check that the Hessian in transformed space = I
*            call CheckHess (Hess,U1,U2)

*           Calculate the gradient of J wrt the transformed space
*           (note U = U2 U1)
            call GradientJTrans (GradJTrans,GradJOb,ObFirst,ObLast,Binv,
     :                           U1,U2,xBp,H,loop,maxobs)

*           The analysis increment in v-space is -GradJTrans = -dJ/dv
*           Transform dJ/dv to x-space
            call Utransform (xANp,GradJTrans,U1,U2,4)

*           Find new linearization state (ie the analysis) xLS=xLS-xANp
            call PlusMinus (xLS,xLS,xANp,-1.,4,4)

*           Put angles into range (-pi to pi)
            xLS(1)=mymod(xLS(1),pi)
            xLS(2)=mymod(xLS(2),pi)

*           Calculate the modsquared of xANp
            ms=modsq(xANp,4)
            print*,'  Outer loop ',loop,', increment len. ',ms
            if ((ms.gt.crit).and.(loop.lt.safety)) goto 10
*           ------------------------------------------------------------
            else
            print*,'NO OBESERVATIONS'
          endif


*         Output analysed initial conditions (with errors)
          write (21,111) 'Anal. init. conds. Time = ',
     :                   real((AssCy-1)*Nwindow)*dt,' with error'
          write (21,112) 180.*xLS(1)/pi, 180.*xLS(2)/pi,
     :                   180.*xLS(3)/pi,180.*xLS(4)/pi
*         Calculate and output the analysis errors (from the Hessian)
*         First calculate the forward model (required to compute Hess.)
          call ForwardModel (xLS,Nsteps,m1,m2,m3,l1,l2,l3,g,dt,pi,
     :                       ObFirst,ObLast,ObTypes,Index,Work1,Work2,
     :                       Work3,H,winsize,maxobs)
          call Hessian (Binv,ObErrs,H,Hess,Obfirst,Oblast,maxobs)
*         Invert the Hessian (use U1 and U2 as workspace)
          call Inverse (Hess,Hessinv,U1,U2,4)
          write (21,112) 180.*sqrt(Hessinv(1,1))/pi,
     :                   180.*sqrt(Hessinv(2,2))/pi,
     :                   180.*sqrt(Hessinv(3,3))/pi,
     :                   180.*sqrt(Hessinv(4,4))/pi


*         Make a forcast based on xLS for next background
          remaining=Nt-(AssCy-1)*Nwindow
          if (remaining.lt.Nwindow) Nsteps=remaining
          call Advance (x,xLS,Nsteps,1,0.,m1,m2,m3,l1,l2,l3,g,dt,pi,
     :                  winsize)
*         Output analysis
          call OutputIntegration (x,Nsteps,Nwindow,AssCy,13,dt,winsize)

*         Calculate B and Binv based on this integration
          call CalcBinv (x,Binv,Nsteps,pi,winsize)

*         Continue forecast with no assimilation
          call Background (x0,xfore,AssCy,x0f,Nwindow,winsize)
          call Advance (xfore,x0f,Nsteps,1,0.,m1,m2,m3,l1,l2,l3,g,dt,pi,
     :                  winsize)
          call OutputIntegration (xfore,Nsteps,Nwindow,AssCy,19,dt,
     :                            winsize)


        enddo

*       Close the output files
        close (13)
        close (14)
        close (15)
        close (16)
        close (17)
        close (18)
        close (19)
        close (21)


      else


*       ----------------------------------------------------------------
*       || *** SECTION III - INSERT OBSERVATIONS ***                  ||
*       ----------------------------------------------------------------
*       This run is to insert observations
*       ----------------------------------
        print*,'INSERT OBSERVATIONS'
        print*,'-------------------'

*       Open the integration output file for the 'insertion' run
        open (13,file='Assim.ins')
*       Open files containing observations used
        open (14,file='Assim.obs1')
        open (15,file='Assim.obs2')
        open (16,file='Assim.obs3')
        open (17,file='Assim.obs4')

*       Calculate the number of assimilation cycles
        Cycles=int(real(Nt)/real(Nwindow))
        if (Cycles*Nwindow.lt.Nt) Cycles=Cycles+1
*       The number of time steps in an assimilation cycle is usually
        Nsteps=Nwindow
*       Point to end of observations (initial value)
        ObLast=0

        do AssCy=1,Cycles
          print*,'-----------------------------------------------------'
          print*,'Assimilation Cycle No. ',AssCy
          print*,'-----------------------------------------------------'

          ObFirst=ObLast+1
*         Check to see if there are < Nwindow left in the time domain
          remaining=Nt-(AssCy-1)*Nwindow
          if (remaining.lt.Nwindow) Nsteps=remaining

*         Work out times bounding this assimilation cycle
          t1=real((AssCy-1)*Nwindow)*dt
          t2=real((AssCy-1)*Nwindow+Nsteps)*dt

*         Find the last observation pertaining to current cycle
*         (and set-up the observation indices)
          call FindLast (ActualObs,ObTimes,ObErrs,ObTypes,ObFirst,
     :                   ObLast,Index,Nobs,t1,t2,dt,14,maxobs)
          print*,'Start and end times : ',t1,t2
          print*,'First ob index,time : ',ObFirst,ObTimes(ObFirst)
          print*,'Last ob index,time  : ',ObLast,ObTimes(ObLast)

          print*,'Inserting ',ObLast-ObFirst+1,' observations'

*         Set up the initial conditions in x array
          call InitialCondsInsert (x0,x,AssCy,Nsteps,winsize)

*         Integrate over time and insert observations as and when
*         they appear.
          loop=ObFirst
          do t=1,Nsteps
*           If there are observations at this time, insert them.
            finish=.false.
20          if (loop.le.ObLast) then
              if (Index(loop).eq.t) then
*               This is observation is for time t - insert it.
                x(t,ObTypes(loop))=ActualObs(loop)
                loop=loop+1
              else
                finish=.true.
              endif
            else
              finish=.true.
            endif
            if (.not.finish) goto 20

*           Integrate one time step forward
            if (t.lt.Nsteps)
     :        call Adv1Step (x,t,m1,m2,m3,l1,l2,l3,g,dt,pi,winsize)

          enddo

*         Output 'analysis'
          call OutputIntegration (x,Nsteps,Nwindow,AssCy,13,dt,winsize)

        enddo

*       Close the output files
        close (13)
        close (14)
        close (15)
        close (16)
        close (17)


      endif
      print*,'Program Completed'
      print*,'========================================================='

110   format ((A),4F18.3)
111   format ((A),F10.4,(A))
112   format (4F18.3)
113   format ((A),3F18.3)

      end

*     ==================================================================
*     ==================================================================

      subroutine Init (pi,crit,safety)
*     Welcome message, and initialize variables
      implicit none
      integer safety
      real pi,crit
      print*,'========================================================='
      print*,'(1+1)D VAR. SYSTEM FOR THE MODIFIED DOUBLE PENDULUM'
      print*,'---------------------------------------------------------'
      print*,'Ross Bannister, Data Assimilation Research Centre'
      print*,'                University of Reading, UK'
      print*,'---------------------------------------------------------'
      pi=4.*atan(1.)
*     Var. convergence conditions
      crit=0.001
*     Limit for number of iterations/cycle
      safety=20
      return
      end
*     ------------------------------------------------------------------

      subroutine ReadConfig (Status,x0,g,dt,Nt,Nout,l1,l2,l3,m1,m2,m3,
     :                       pi,noisepos,noisevel,ActualObs,ObTypes,
     :                       ObTimes,ObErrs,Nobs,Nwindow,maxobs)
*     To read-in how the code is to be configured
*     To read-in the configuration (and observations if required)
      implicit none
      integer Nt,Nout,Nobs,maxobs,ObTypes(maxobs),Nwindow,lp
      real x0(4),g,dt,l1,l2,l3,m1,m2,m3,noisepos,noisevel,pi
      real ActualObs(maxobs),ObTimes(maxobs),ObErrs(maxobs)
      character*7 Status
      character*10 Temp

*     Nt     : length of integration in time steps
*     Nout   : Frequency of observation output (when creating obs)
*     Nobs   : the number of observations
*     Nwindow: the size of the assimilation 'windows' in time steps
      print 100,'Double Pendulum 4d-Var Assimilation'
      print 100,'Reading in Configuration'
      open (12,file='Assim.conf')
*     Read-in general parameters
      read (12,100) Temp
      read (12,100) Temp
      read (12,100) Temp
      read (12,100) Temp
      read (12,100) Temp
      read (12,100) Status
      print 100,'Status           :  '//Status
      read (12,100) Temp
      read (12,101) l1,l2,l3
      print 111, 'Lengths of arms  : ',l1,l2,l3
      read (12,100) Temp
      read (12,101) m1,m2,m3
      print 111, 'Masses of nodes  : ',m1,m2,m3
      read (12,100) Temp
      read (12,102) g
      print 111, 'g                : ',g
      read (12,100) Temp
      read (12,102) dt
      print 112,'Time step        : ',dt
      read (12,100) Temp
      read (12,103) Nt
      print 113,'Total No of steps: ',Nt
      read (12,100) Temp
      read (12,100) Temp
      read (12,104) x0(1),x0(2),x0(3),x0(4)
      print 114,'Initial state    : ',x0(1),x0(2),x0(3),x0(4)
*     Convert to radians and radians per second
      do lp=1,4
        x0(lp)=pi*x0(lp)/180.
      enddo
*     Read-in parameters to do with making observations
      read (12,100) Temp
      read (12,100) Temp
      read (12,100) Temp
      read (12,103) Nout
      if (Status.eq.'MakeObs') print 113,'Ob create frq: ',Nout
      read (12,100) Temp
      read (12,102) noisepos
      if (Status.eq.'MakeObs') print 112,'Pos. noise   : ',noisepos
*     Convert to radians
      noisepos=pi*noisepos/180.
      read (12,100) Temp
      read (12,102) noisevel
      if (Status.eq.'MakeObs') print 112,'Noise on roc : ',noisevel
*     Convert to radians per second
      noisevel=pi*noisevel/180.
      if (Status.ne.'MakeObs') then
*       Read-in parameters to do with the observations
        read (12,100) Temp
        read (12,100) Temp
        read (12,100) Temp
        read (12,103) Nwindow
        print 113,'Steps/cycle      : ',Nwindow
        read (12,100) Temp
        read (12,100) Temp
        read (12,100) Temp
        print 100,'OBSERVATIONS     :'
*       Loop round and read-in the observations
        Nobs=0
10      Nobs=Nobs+1
          if (Nobs.gt.maxobs) then
            print*,'Insufficient space to hold observations.'
            print*,'Increase value of maxobs in code.'
            stop
          endif
          read (12,105) ActualObs(Nobs),ObTypes(Nobs),ObErrs(Nobs),
     :                  ObTimes(Nobs)
          print 105,ActualObs(Nobs),ObTypes(Nobs),ObErrs(Nobs),
     :              ObTimes(Nobs)
          if ((ActualObs(Nobs).eq.0.).or.(ObTypes(Nobs).eq.0)) goto 11
        goto 10
11      Nobs=Nobs-1
      endif
      close (12)

100   format ((A))
101   format (3F)
102   format (F)
103   format (I)
104   format (4F)
105   format (F9.4,I4,F9.4,F12.6)
111   format ((A),3F)
112   format ((A),F)
113   format ((A),I)
114   format ((A),4F)
      return
      end
*     ------------------------------------------------------------------

      subroutine Equal (x,y,N,size)
*     Set the vector x=y
      implicit none
      integer size,N,i
      real x(size),y(size)
      do i=1,N
        x(i)=y(i)
      enddo
      return
      end
*     ------------------------------------------------------------------

      subroutine PlusMinus (x,y,z,s,N,size)
*     x=y+sz (x,y,z, vectors)
      implicit none
      integer size,N,i
      real x(size),y(size),z(size),s
      do i=1,N
        x(i)=y(i)+s*z(i)
      enddo
      return
      end
*     ------------------------------------------------------------------

      subroutine SetZero (x,N,size)
*     Set the vector x=0
      implicit none
      integer size,N,i
      real x(size)
      do i=1,N
        x(i)=0.
      enddo
      return
      end
*     ------------------------------------------------------------------

      subroutine MakeObs (x0,g,dt,pi,Nt,Nout,l1,l2,l3,m1,m2,m3,noisepos,
     :                    noisevel)
*     To run the double pendulum, outputting a set of simulated obs
*     x0 is the initial conditions
*     Nout is the frequency that the modified state is output
*     noise is the range of the noise distribution

      implicit none
      integer Nt,Nout,t,lp,seed,lim
      real x0(4),x(4),g,dt,pi,l1,l2,l3,m1,m2,m3,noisepos,noisevel,err(4)
      real ob,ran,ran0

*     Seed the random number generator
      seed=549872180

*     Transfer some data around
      err(1)=noisepos
      err(2)=noisepos
      err(3)=noisevel
      err(4)=noisevel

*     Output angular positions, and roc of angles or positions only?
      if (abs(noisevel).lt.(0.0000001)) then
        lim=2
        else
        lim=4
      endif

*     Put initial conditions into state vector
      do lp=1,4
        x(lp)=x0(lp)
      enddo

      open (12,file='Observations')
      open (13,file='Assim.ObsBase')
      write (12,100) '-----------------------------------------'
      do t=1,Nt
        call RungeKutt4 (x,m1,m2,m3,l1,l2,l3,g,dt,pi)
        write (13,102) real(t-1)*dt,x(1),x(2),x(3),x(4)
        if (mod(t,Nout).eq.0) then
*         Make up some 'observations', lim is either 2 or 4 (see above)
          do lp=1,lim
*           Throw a random number to add noise to model value
            ran=ran0(seed)
            ob=2.*(ran-0.5)*err(lp)+x(lp)
            write (12,101) ob,lp,err(lp),real(t-1)*dt
          enddo
        endif
      enddo
      write (12,101) 0.,0,0.,0.
      write (12,100) '-----------------------------------------'
      close (12)
      close (13)

100   format ((A))
101   format (F9.4,I4,F9.4,F12.6)
102   format (5F12.6)
      return
      end
*     ------------------------------------------------------------------

      subroutine FindLast (ActualObs,ObTimes,ObErrs,ObTypes,First,Last,
     :                     Index,Nobs,tfirst,tlast,dt,unit,maxobs)
*     To find the index of the last observation in the current cycle
*     Also to mark the presence of each observation at the appropriate
*     time in the array 'index'
*     Also to output the observations to special files
*     types 1 => unit,   types 2 => unit+1,
*     types 3 => unit+2, types 4 => unit+3  (unit is a file)
      implicit none
      integer maxobs,First,Last,Index(maxobs),Nobs,i,ObTypes(maxobs)
      integer unit
      real ActualObs(maxobs),ObTimes(maxobs),ObErrs(maxobs)
      real tfirst,tlast,dt

*     Check that the next observation is actually in this window
      if ((ObTimes(First).ge.tfirst).and.(ObTimes(First).lt.tlast)) then

*       Find extent of observations
        Last=First-1
10      if ((ObTimes(Last+1).ge.tlast).or.(Last.eq.Nobs)) goto 11
        Last=Last+1
        goto 10
11      continue

*       Fill-out the Index array
        if (Last.ge.First) then
          do i=First,Last
            Index(i)=INT((ObTimes(i)-tfirst+dt/2.)/dt)
            write (unit+ObTypes(i)-1,101) ObTimes(i),ActualObs(i),
     :                                    ObErrs(i)
          enddo
        endif
      endif

101   format (3F12.6)
      return
      end
*     ------------------------------------------------------------------

      subroutine Advance (yt,y,Nt,i,dyi,m1,m2,m3,l1,l2,l3,g,dt,pi,
     :                    winsize)
*     To advance the system through Nt time steps
*     Element i of the initial state is incremented by dyi
      implicit none
      integer Nt,i,winsize,j,t
      real yt(winsize,4),dyi,m1,m2,m3,l1,l2,l3,g,dt,pi,ytemp(4),y(4)

*     Setup the initial state
      do j=1,4
        yt(1,j)=y(j)
        ytemp(j)=y(j)
      enddo
      yt(1,i)=y(i)+dyi
      ytemp(i)=y(i)+dyi

*     Integrate over time
      do t=2,Nt
*       Advance in time
        call RungeKutt4 (ytemp,m1,m2,m3,l1,l2,l3,g,dt,pi)
*       Copy resulting state into 2d structure
        do j=1,4
          yt(t,j)=ytemp(j)
        enddo
      enddo

      return
      end
*     ------------------------------------------------------------------

      subroutine Adv1Step (x,t,m1,m2,m3,l1,l2,l3,g,dt,pi,winsize)
*     To advance the state by one time step, from t to t+1
      implicit none
      integer t,winsize,i
      real x(winsize,4),m1,m2,m3,l1,l2,l3,g,dt,pi,xtemp(4)

      do i=1,4
        xtemp(i)=x(t,i)
      enddo
      call RungeKutt4 (xtemp,m1,m2,m3,l1,l2,l3,g,dt,pi)
      do i=1,4
        x(t+1,i)=xtemp(i)
      enddo
      return
      end
*     ------------------------------------------------------------------

      subroutine Background (x0,x,AssCy,xb,Nwindow,winsize)
*     To set up the background state.
*     The background state is taken from the 'last' element of x
*     if AssCy>1.  If AssCy=1 then take the background from x0.
      implicit none
      integer AssCy,Nwindow,winsize,i
      real x0(4),x(winsize,4),xb(4)
      if (AssCy.eq.1) then
*       Use x0 as background
        do i=1,4
          xb(i)=x0(i)
        enddo
        else
*       Use the 'last' element of x
        do i=1,4
          xb(i)=x(Nwindow,i)
        enddo
      endif
      return
      end
*     ------------------------------------------------------------------

      subroutine InitialCondsInsert (x0,x,AssCy,Nwindow,winsize)
*     To set up the initial conditions for the 'insert' mode of the code
*     The initial state is taken from the 'last' element of x
*     if AssCy>1.  If AssCy=1 then take the values from x0.
      implicit none
      integer AssCy,Nwindow,winsize,i
      real x0(4),x(winsize,4)
      if (AssCy.eq.1) then
*       Use x0 as initial conditions
        do i=1,4
          x(1,i)=x0(i)
        enddo
        else
*       Use the 'last' element of x
        do i=1,4
          x(1,i)=x(Nwindow,i)
        enddo
      endif
      return
      end
*     ------------------------------------------------------------------

      subroutine GradientJOb (x,ActualObs,ObType,ObErrs,Index,GradJOb,
     :                        Jo,First,Last,pi,winsize,maxobs)
*     To calculate E^(-1) (obs-modelobs) and Jo
      implicit none
      integer maxobs,winsize,Index(maxobs),ObType(maxobs),First,Last,i
      real x(winsize,4),ActualObs(maxobs),ObErrs(maxobs),GradJOb(maxobs)
      real diff,pi,Jo

      Jo=0.
      do i=First,Last
*       Work out the difference
        diff=ActualObs(i)-x(Index(i),ObType(i))
*       if the observation is an angle, need to do something special
        if (ObType(i).le.2) then
*         The observation - model prediction may need adjusting
          call AdjustDiff (diff,pi)
        endif
        GradJob(i-First+1)=-1.*diff/(ObErrs(i)*ObErrs(i))
        Jo=Jo-diff*GradJob(i-First+1)
      enddo
      Jo=Jo/2.

      return
      end
*     ------------------------------------------------------------------

      subroutine CalcJb (xB,x,Jb,Binv,winsize)
*     Calculation of Jb (Jo is found within subroutine GradientJOb
      implicit none
      integer winsize,i,j
      real xB(4),x(winsize,4),Jb,Binv(4,4),sum

*     Jb = tpose(xB-x(1)) Binv (xB-x(1))
      Jb=0.
      do i=1,4
        sum=0.
        do j=1,4
          sum=sum+Binv(i,j)*(xB(j)-x(1,j))
        enddo
        Jb=Jb+(xB(i)-x(1,i))*sum
      enddo
      Jb=Jb/2.

      return
      end
*     ------------------------------------------------------------------

      subroutine ForwardModel (xLS,Nsteps,m1,m2,m3,l1,l2,l3,g,dt,pi,
     :                         First,Last,ObTypes,Index,MobsMinus,
     :                         MobsPlus,xt,H,winsize,maxobs)
*     Compute the forward model matrix dyob(j)/dxi
      implicit none
      integer i,j,Nsteps,winsize,maxobs,ObTypes(maxobs),Index(maxobs)
      integer First,Last
      real xLS(4),xt(winsize,4),delta,m1,m2,m3,l1,l2,l3,g,dt,pi
      real MobsMinus(maxobs),MobsPlus(maxobs),H(maxobs,4)

      delta=0.0001
      do i=1,4
*       Do two integrations by adjusting state i

*       Decrement element i by -delta
*       算出窗口每个时刻的状态
        call Advance (xt,xLS,Nsteps,i,-1.*delta,m1,m2,m3,l1,l2,l3,g,dt,
     :                pi,winsize)
*       Contruct a vector of observations from this forecast
        call FindModelObs (MobsMinus,ObTypes,Index,xt,First,Last,maxobs,
     :                     winsize)
*       Increment element i by delta
        call Advance (xt,xLS,Nsteps,i,delta,m1,m2,m3,l1,l2,l3,g,dt,pi,
     :                winsize)
*       Contruct a vector of observations from this forecast
        call FindModelObs (MobsPlus,ObTypes,Index,xt,First,Last,maxobs,
     :                     winsize)
*       Now fill out the forward model matrix
        do j=1,Last-First+1
          H(j,i)=(MobsPlus(j)-MobsMinus(j))
          if (ObTypes(j+First-1).le.2) then
*           The model observation difference may need adjusting
            call AdjustDiff (H(j,i),pi)
          endif
          H(j,i)=H(j,i)/(2.*delta)
        enddo
      enddo

      return
      end
*     ------------------------------------------------------------------

      subroutine FindModelObs (ModelObs,ObTypes,Index,x,first,last,
     :                         maxobs,winsize)
*     Construct a vector of model observations
      implicit none
      integer maxobs,winsize,Index(maxobs),ObTypes(maxobs),i,first,last
      real ModelObs(maxobs),x(winsize,4)
      do i=first,last
        ModelObs(i-first+1)=x(Index(i),ObTypes(i))
      enddo
      return
      end
*     ------------------------------------------------------------------

      subroutine Hessian (Binv,ObErrs,H,Hess,first,last,maxobs)
*     To compute the Hessian in model (physical) space
      implicit none
      integer first,last,maxobs,i,j,k,l
      real Binv(4,4),ObErrs(maxobs),H(maxobs,4),Hess(4,4),sum

      do i=1,4
        do j=1,4
          sum=0.
          do k=first,last
            l=k-first+1
            sum=sum+H(l,i)*H(l,j)/(ObErrs(k)*ObErrs(k))
          enddo
          Hess(i,j)=Binv(i,j)+sum
        enddo
      enddo
      return
      end
*     ------------------------------------------------------------------

      subroutine Transpose (Matrix,size)
*     Transpose a matrix
      implicit none
      integer size,i,j
      real Matrix(size,size),element
      do i=1,size-1
        do j=i+1,size
          element=Matrix(i,j)
          Matrix(i,j)=Matrix(j,i)
          Matrix(j,i)=element
        enddo
      enddo
      return
      end
*     ------------------------------------------------------------------

      subroutine GradientJTrans (GradJTrans,GradJOb,first,last,Binv,U1,
     :                           U2,xBp,H,iteration,maxobs)
*     To calculate the gradient of J wrt the transformed space vector
      implicit none
      integer first,last,NumObCycle,i,k,l,iteration,maxobs
      real GradJTrans(4),GradJOb(maxobs),Binv(4,4),U1(4,4),U2(4,4)
      real H(maxobs,4),xBp(4),sumb,sumo,sum,temp(4)

*     The number of observations in this cycle
      NumObCycle=last-first+1

*     Stage 1: Compute -(B^-1)xBp + tpose(H)GradJOb
      do k=1,4
        sumb=0.
*       If iteration=1 then know that sumb=0 as xBp=0
        if (iteration.gt.1) then
          do l=1,4
            sumb=sumb+Binv(k,l)*xBp(l)
          enddo
        endif
        sumo=0.
        do l=1,NumObCycle
          sumo=sumo+H(l,k)*GradJob(l)
        enddo
        temp(k)=sumo-sumb
      enddo

*     Stage 2: Compute GradJTrans (apply tpose(U) transform)
      do i=1,4
        sum=0.
        do k=1,4
          sum=sum+U2(k,i)*temp(k)
        enddo
        GradJTrans(i)=U1(i,i)*sum
      enddo

      return
      end
*     ------------------------------------------------------------------

      subroutine Utransform (dx,dv,U1,U2,size)
*     Performs U-transform, dx=Udv,  U = U2 U1
      implicit none
      integer i,j,size
      real dx(size),dv(size),U1(size,size),U2(size,size),sum
      do i=1,size
        sum=0.
        do j=1,size
          sum=sum+U2(i,j)*U1(j,j)*dv(j)
        enddo
        dx(i)=sum
      enddo
      return
      end
*     ------------------------------------------------------------------

      real function modsq (x,size)
*     Calculate transpose(x)(x)
      implicit none
      integer size,i
      real sum,x(size)
      sum=0.
      do i=1,size
        sum=sum+x(i)*x(i)
      enddo
      modsq=sum
      return
      end
*     ------------------------------------------------------------------

      subroutine CalcBinv (x,Binv,Nsteps,pi,winsize)
*     To calculate B (time average correlation), returning its inverse
      implicit none
      integer Nsteps,winsize,i,j,k
      real x(winsize,4),B(4,4),Binv(4,4),Mean(4),sum,diff1,diff2,pi
      real Temp1(4,4),Temp2(4,4),sgn,inc

*     Calculate the mean value of each degree of freedom
      do i=1,4
        sum=0.
        sgn=sign(1.,x(1,i))
        do k=1,Nsteps
          inc=x(k,i)
*         For angles, if the sign has flipped - need to adjust inc
          if (i.le.2) then
            if (abs(sign(1.,inc)-sgn).gt.(0.1)) inc=inc+2.*pi*sgn
          endif
          sum=sum+inc
        enddo
        Mean(i)=sum/real(Nsteps)
      enddo

*     Calculate the mean correlations (B)
      do i=1,4
        do j=1,4
          sum=0.
          do k=1,Nsteps
            diff1=x(k,i)-Mean(i)
            if (i.le.2) call AdjustDiff (diff1,pi)
            diff2=x(k,j)-Mean(j)
            if (j.le.2) call AdjustDiff (diff2,pi)
            sum=sum+diff1*diff2
          enddo
          B(i,j)=sum/real(Nsteps)
        enddo
      enddo

*     Invert B for Binv
      call Inverse (B,Binv,Temp1,Temp2,4)
      return
      end
*     ------------------------------------------------------------------

      subroutine InvSqrt (A,size)
*     Invert and square root the diagonal matrix A
      implicit none
      integer i,size
      real A(size,size)
      do i=1,size
        A(i,i)=1./sqrt(A(i,i))
      enddo
      return
      end
*     ------------------------------------------------------------------

      subroutine OutputIntegration (xt,Nsteps,Nwindow,AssCycle,unit,
     :                              dt,winsize)
*     To output the result of an integration
      implicit none
      integer Nsteps,Nwindow,AssCycle,unit,winsize,t
      real xt(winsize,4),time,dt
      do t=1,Nsteps
        time=real((AssCycle-1)*Nwindow+t-1)*dt
        write (unit,101) time,xt(t,1),xt(t,2),xt(t,3),xt(t,4)
      enddo
101   format (5F12.6)
      return
      end
*     ------------------------------------------------------------------

      subroutine CheckPos (U1,size)
*     To check sign of eigenvalues (diagonal elements of U1)
      implicit none
      integer size,i
      real U1(size,size)
      do i=1,size
        if (U1(i,i).le.(0.)) then
          print*,'  Negative eigenvalue ',U1(i,i)
          print*,'  ABORTING CODE'
          stop
        endif
      enddo
      return
      end
*     ------------------------------------------------------------------

      subroutine BoundTime (time,unit)
*     To output the time at the start of the boundary
      implicit none
      integer unit
      real time
      write (unit,101) time,0.
101   format (F12.6,F4.1)
      return
      end
*     ------------------------------------------------------------------

      subroutine AdjustDiff (diff,range)
*     to check, adn adjust if necessary, angular ranges
      implicit none
      real diff,range
      if (diff.ge.range) then
        diff=diff-2.*range
      else if (diff.lt.(-1.*range)) then
        diff=diff+2.*range
      endif
      return
      end
*     ------------------------------------------------------------------

      subroutine SetUpName (AssCy,outerloop,fname)
*     To concoct a filename for the output of the current integration
      implicit none
      integer AssCy,outerloop
      character*8 fname
      character*2 part1,part2
      if (AssCy.lt.10) then
        write (part1,'(A1,I1)') '0',AssCy
        else
        write (part1,'(I2)') AssCy
      endif
      if (outerloop.lt.10) then
        write (part2,'(A1,I1)') '0',outerloop
        else
        write (part2,'(I2)') outerloop
      endif
      fname='Cy'//part1//'Lp'//part2
      return
      end

*     ------------------------------------------------------------------

      subroutine CheckHess (Hess,U1,U2)
*     Temporary routine to check that Hessian in trans. space = I
*     Also check that U2 is orthonormal
      implicit none
      integer i,j,l,m
      real Hess(4,4),Htrans(4,4),U1(4,4),U2(4,4),sum1,sum2

      print*,'  Check that U2 is orthonormal'
      do i=1,4
        do j=1,4
          sum1=0.
          do l=1,4
            sum1=sum1+U2(l,i)*U2(l,j)
          enddo
*         Use Htrans space for convenience
          Htrans(i,j)=sum1
        enddo
        print*,'  U2TU2 ',(Htrans(i,j),j=1,4)
      enddo

      print*,'  Check that Hess gives unit matrix'
      do i=1,4
        do j=1,4
          sum1=0.
          do l=1,4
            sum2=0.
            do m=1,4
              sum2=sum2+Hess(l,m)*U2(m,j)
            enddo
            sum1=sum1+U2(l,i)*sum2
          enddo
          Htrans(i,j)=U1(i,i)*U1(j,j)*sum1
        enddo
        print*,'  H-trans2 ',(Htrans(i,j),j=1,4)
      enddo
      return
      end
*     ------------------------------------------------------------------

      subroutine CheckDiag (Hess,U2)
*     Temporary routine to check that U2 is the evector matrix of Hess
      implicit none
      integer i,j,k,l
      real Hess(4,4),Htrans(4,4),U2(4,4),sum1,sum2


      print*,'  Check that eigenvectors'
      do i=1,4
        do j=1,4
          sum1=0.
          do k=1,4
            sum2=0.
            do l=1,4
              sum2=sum2+Hess(k,l)*U2(j,l)
            enddo
            sum1=sum1+U2(i,k)*sum2
          enddo
          Htrans(i,j)=sum1
        enddo
        print*,'  H-trans1 ',(Htrans(i,j),j=1,4)
      enddo
      return
      end
*     ------ FORWARD MODEL ROUTINES ------------------------------------

      subroutine matrixA (y,m1,m2,m3,l1,l2,l3,A)
*     To determine the elements of matrix A
      implicit none
      real y(4),m1,m2,m3,l1,l2,l3,A(2,2)
      real th1,th2,al1,al2,term1,term2,term3
*     Angle 1
      th1=y(1)
*     Angle 2
      th2=y(2)
*     Rate of change angle 1
      al1=y(3)
*     Rate of change angle 2
      al2=y(4)

      term1=l2*l3*m3*sin(th1-th2)
      term2=al2*term1
      term3=al1*term1

      A(1,1)=-1.*term2
      A(1,2)=term2
      A(2,1)=-1.*term3
      A(2,2)=term3

      return
      end
*     ------------------------------------------------------------------

      subroutine matrixBinv (y,m1,m2,m3,l1,l2,l3,Binv)
*     To determine the elements of matrix B-inverse
      implicit none
      real y(4),m1,m2,m3,l1,l2,l3,Binv(2,2)
      real B(2,2),th1,th2,al1,al2,term1,det
*     Angle 1
      th1=y(1)
*     Angle 2
      th2=y(2)
*     Rate of change angle 1
      al1=y(3)
*     Rate of change angle 2
      al2=y(4)

      term1=m3*l2*l3*cos(th1-th2)

*     Construct the actual matrix B ...
      B(1,1)=m1*l1*l1+(m2+m3)*l2*l2
      B(1,2)=term1
      B(2,1)=term1
      B(2,2)=m3*l3*l3

*     ... but we need its inverse
      det=B(1,1)*B(2,2)-B(1,2)*B(2,1)
      Binv(1,1)=B(2,2)/det
      Binv(1,2)=-1.*B(1,2)/det
      Binv(2,1)=-1.*B(2,1)/det
      Binv(2,2)=B(1,1)/det

      return
      end
*     ------------------------------------------------------------------

      subroutine vectorf (y,m1,m2,m3,l1,l2,l3,g,f)
*     To determine the elements of vector f
      implicit none
      real y(4),m1,m2,m3,l1,l2,l3,g,f(2)
      real th1,th2,al1,al2,term1,term2
*     Angle 1
      th1=y(1)
*     Angle 2
      th2=y(2)
*     Rate of change angle 1
      al1=y(3)
*     Rate of change angle 2
      al2=y(4)

      term1=m3*l2*l3*al1*al2*sin(th1-th2)
      term2=sin(th1)

      f(1)=-1.*term1+g*term2*(m1*l1-(m2+m3)*l2)
      f(2)=term1-g*m3*l3*sin(th2)

      return
      end
*     ------------------------------------------------------------------

      subroutine operate (A,Binv,f,dt,y,kk)
*     To perform operation k=Ty+g where T(A,Binv) and g(f)
*     (see notes)
      implicit none
      integer i,j,k
      real A(2,2),Binv(2,2),f(2),dt,y(4),kk(4)
      real Prod(2,2),total,vec(2)

*     Calculate -Binv A
      do i=1,2
        do j=1,2
          total=0.
          do k=1,2
            total=total-Binv(i,k)*A(k,j)
          enddo
          Prod(i,j)=total
        enddo
      enddo

*     Calculate Binv f
      do i=1,2
        total=0.
        do j=1,2
          total=total+Binv(i,j)*f(j)
        enddo
        vec(i)=total
      enddo

*     Now operate
      kk(1)=dt*y(3)
      kk(2)=dt*y(4)
      kk(3)=dt*(Prod(1,1)*y(3)+Prod(1,2)*y(4)+vec(1))
      kk(4)=dt*(Prod(2,1)*y(3)+Prod(2,2)*y(4)+vec(2))

      return
      end
*     ------------------------------------------------------------------

      subroutine sum (a,b,c,size,flag)
*     Vector sum c=a+b   (flag=1)
*     Vector sum c=a+b/2 (flag=2)
      implicit none
      integer size,i,flag
      real a(size),b(size),c(size)
      if (flag.eq.1) then
        do i=1,size
          c(i)=a(i)+b(i)
        enddo
        else
        do i=1,size
          c(i)=a(i)+b(i)/2.
        enddo
      endif
      return
      end
*     ------------------------------------------------------------------

      subroutine RungeKutt4 (y,m1,m2,m3,l1,l2,l3,g,dt,pi)
*     Fourth order Runge-Kutta for double pendulum
      implicit none
      integer i
      real y(4),m1,m2,m3,l1,l2,l3,g,dt,pi
      real k1(4),k2(4),k3(4),k4(4),A(2,2),Binv(2,2),f(2),ytemp(4),mymod

*     Stage 1
*     -------
      call matrixA (y,m1,m2,m3,l1,l2,l3,A)
      call matrixBinv (y,m1,m2,m3,l1,l2,l3,Binv)
      call vectorf (y,m1,m2,m3,l1,l2,l3,g,f)
      call operate (A,Binv,f,dt,y,k1)

*     Stage 2
*     -------
      call sum (y,k1,ytemp,4,2)
      call matrixA (ytemp,m1,m2,m3,l1,l2,l3,A)
      call matrixBinv (ytemp,m1,m2,m3,l1,l2,l3,Binv)
      call vectorf (ytemp,m1,m2,m3,l1,l2,l3,g,f)
      call operate (A,Binv,f,dt,ytemp,k2)

*     Stage 3
*     -------
      call sum (y,k2,ytemp,4,2)
      call matrixA (ytemp,m1,m2,m3,l1,l2,l3,A)
      call matrixBinv (ytemp,m1,m2,m3,l1,l2,l3,Binv)
      call vectorf (ytemp,m1,m2,m3,l1,l2,l3,g,f)
      call operate (A,Binv,f,dt,ytemp,k3)

*     Stage 4
*     -------
      call sum (y,k3,ytemp,4,1)
      call matrixA (ytemp,m1,m2,m3,l1,l2,l3,A)
      call matrixBinv (ytemp,m1,m2,m3,l1,l2,l3,Binv)
      call vectorf (ytemp,m1,m2,m3,l1,l2,l3,g,f)
      call operate (A,Binv,f,dt,ytemp,k4)

*     Add-up for incremented y-vector
      do i=1,4
        y(i)=y(i) + k1(i)/6. + k2(i)/3. + k3(i)/3. + k4(i)/6.
      enddo

*     Put angles back into range (-pi to pi)
      do i=1,2
        y(i)=mymod(y(i),pi)
      enddo

      return
      end
*     ------------------------------------------------------------------

      real function mymod (x,interval)
*     This function puts x into the interval -interval<=x<interval
      implicit none
      real x,interval,sgn,z
      z=x+interval
      sgn=sign(1.0,z)
      z=abs(z)
      z=z-2.*interval*real(int(z/(2.*interval)))
      z=sgn*(z)-interval
      if (z.lt.(-interval)) z=z+2.*interval
      mymod=z

      return
      end
*     ==================================================================

      subroutine Inverse (A,Ainv,WS1,WS2,size)
*     Computes the inverse of A via Jacobi diagonalization
*     WS1 and WS2 are workspace arrays
*     Not efficient for large matrices
      implicit none
      integer size,i,j,k
      real A(size,size),Ainv(size,size),WS1(size,size),WS2(size,size),s
*     Diagonalize A
      call Jacobi (A,WS1,WS2,size,size,0.00001)
      do i=1,size
        do j=1,size
          s=0.
          do k=1,size
            s=s+WS2(k,i)*WS2(k,j)/WS1(k,k)
          enddo
          Ainv(i,j)=s
        enddo
      enddo
      return
      end
*     ==================================================================

      subroutine Jacobi (A,B,T,N,max,lim)
*     To diagonalize a square symmetric matrix by a sequence of Jacobi
*     rotations
*     Ross Bannister, January 2001
*     A (input) matrix in original representation
*     B (output) matrix in diagonal representation
*     T (output) matrix of eigenvectors (rows)
*     N (input) the order of the matrix
*     max (input) the maximum order of the matrix (for array dim)
*     lim (input) convergence criteria, e.g. lim=0.01

      implicit none
      integer N,max,i,j,p,q,k,sweeps
*      integer lp1,lp2
      real A(max,max),B(max,max),T(max,max),lim,off,diag,theta,t1,t2,tt
      real cc,ss,sc,Bip,Bqj,Tpj,Bpp,Bqq,dis,s,c,norm

      if (N.gt.max) then
        print*,'Report from matrix JACOBI.'
        print*,'N =',N
        print*,'max =',max
        print*,'Please increase the value of max.'
        stop
      endif

*     The matrix B is initially equal to A
      do j=1,N
        do i=1,N
          B(i,j)=A(i,j)
        enddo
      enddo

*     The matrix T is initially the identity matrix
      do j=1,N
        do i=1,N
          if (i.eq.j) then
            T(i,j)=1.0
            else
            T(i,j)=0.0
          endif
        enddo
      enddo

*     Calculate the sum of the square of off-diagonal elements of B
      off=0.0
      do j=2,N
        do i=1,j-1
          off=off+B(i,j)*B(i,j)
        enddo
      enddo
      off=2.0*off

*     Calculate the sum of the square of diagonal elements of B
      diag=0.0
      do i=1,N
        diag=diag+B(i,i)*B(i,i)
      enddo

*     The number of sweeps completed
      sweeps=0
10    continue

*     Loop round each upper-right element (a 'sweep')
      do p=2,N
        do q=1,p-1
          if (abs(B(q,p)).gt.(0.0000001)) then
*           In order to find the new representation of the matrix,
*           work out some details regarding the rotation, U to eliminate
*           element q,p of B (p>q).
            theta=(B(p,p)-B(q,q))/(2.0*B(q,p))
            dis=sqrt(theta*theta+1)
*           Two possible values of t
            t1=-theta+dis
            t2=-theta-dis
*           Choose the smaller value
            if (abs(t1).lt.abs(t2)) then
              tt=t1
              else
              tt=t2
            endif
*           Find the values of s and c
            c=1/sqrt(tt*tt+1)
            s=tt*c
*           Modify the operator UBU(trans) (piecewise)
*           Update only upper-right part of matrix
            cc=c*c
            ss=s*s
            sc=s*c
*           Columns p and q:
            do i=1,p-1
              if (i.ne.q) then
                Bip=B(i,p)
                B(i,p)=c*Bip+s*B(i,q)
                if (i.le.q) B(i,q)=c*B(i,q)-s*Bip
              endif
            enddo
*           Rows p and q:
            do j=q+1,N
              if (j.ne.p) then
                Bqj=B(q,j)
                B(q,j)=c*Bqj-s*B(p,j)
                if (j.gt.p) B(p,j)=c*B(p,j)+s*Bqj
              endif
            enddo
*           Mixed elements
*            B(q,p)=(cc-ss)*B(q,p)+sc*(B(q,q)-B(p,p))
            B(q,p)=0.0

*           New diagonal elements
            Bpp=B(p,p)
            Bqq=B(q,q)
            B(p,p)=cc*Bpp+ss*B(q,q)+2.0*sc*B(p,q)
            B(q,q)=cc*B(q,q)+ss*Bpp-2.0*sc*B(p,q)

*           Modified value of sum of square of off-diagonal elements...
            off=off-2.0*B(p,q)*B(p,q)
*           ...and diagonal elements
            diag=diag-Bpp*Bpp-Bqq*Bqq+B(p,p)*B(p,p)+B(q,q)*B(q,q)

*           Update lower-left of matrix
            do k=1,N
              if (k.lt.p) then
                B(p,k)=B(k,p)
                if (k.lt.q) B(q,k)=B(k,q)
              endif
              if (k.gt.q) then
                B(k,q)=B(q,k)
                if (k.gt.p) B(k,p)=B(p,k)
              endif
            enddo

*           Modify the transformation T(new)=U T(old)
*           Only rows p and q are affected
            do j=1,N
              Tpj=T(p,j)
              T(p,j)=c*Tpj+s*T(q,j)
              T(q,j)=c*T(q,j)-s*Tpj
            enddo
          endif
        enddo
      enddo

      sweeps=sweeps+1
      if (((off/diag).gt.lim).or.(sweeps.lt.3)) goto 10

*     Normalize the eigenvectors
      do i=1,N
        norm=0.0
        do j=1,N
          norm=norm+T(i,j)*T(i,j)
        enddo
        norm=sqrt(norm)
        do j=1,N
          T(i,j)=T(i,j)/norm
        enddo
      enddo

      return
      end
*     ==================================================================

      real function ran0(idum)
      implicit none
      integer idum,IA,IM,IQ,IR,MASK
      real AM
      parameter (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     :           MASK=123459876)
*     Minimal random number generator of Park and Miller.  Returns a
*     uniform random deviate between 0.0 and 1.0.  Set or reset idum
*     to any integer value (except the unlikely value MASK) to
*     initialize the sequence; idum must not be altered between calls
*     for successive deviates in a sequence.
*     Numerical Recipies in Fortran, page 270
      integer k
      idum=ieor(idum,MASK)
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      ran0=AM*idum
      idum=ieor(idum,MASK)
      return
      end
*     ==================================================================