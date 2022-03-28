! FOBOS-binaries.f90 takes two epochs of observation of a two body system and fits
    ! the orbital parameters (semi-major axis, eccentricity, inclination, orientation, and
    ! mean anomaly). The code is currently set up to read two separations (in milliarcseconds,
    ! column 2) and two position angles (in degrees, column 4), plus their associated respective
    ! errors (columns 3 and 5) from a text file called "input". The input file will have the format

    ! t1 (yrs)  sep1 (mas)  sep_err1 (mas)  pa1 (deg)   pa_err1 (deg)
    ! t2 (yrs)  sep2 (mas)  sep_err2 (mas)  pa2 (deg)   pa_err2 (deg)

    ! The value for t1 will typically be zero and t2 the number of years between
    ! the first and second epochs of observation.

    ! The target name, which affects the name of the output files, is given in a text file
    ! called "sys-name". The masses and distance to the objects are contained in the file
    ! "mass-dist", which has the format

    ! M1 M1_err  M2  M2_err  d   d_err

    ! All masses will be given in solar masses and distances in parsecs

PROGRAM binaries
    IMPLICIT NONE
    INCLUDE "omp_lib.h"

    INTEGER:: n,i,idum,kdum,match,nmatch
    INTEGER:: nt,my_id,count,loop
    REAL:: M(1:2),M_err(1:2),dtobs,amin,amax,vobs,d,d_err,er
    REAL:: obs(1:3),err(1:3),params(1:5),sep(1:3)
    REAL:: seps(1:2),sep_err(1:2),pa(1:2),pa_err(1:2)
    REAL:: pi,au,year,msun,G,d2r,dmoved
    REAL:: t_run,ran3
    INTEGER,ALLOCATABLE:: rinit(:)
    !Parameters to read in from file
    REAL:: obsdata(1:5,1:2)
    REAL:: properties(1:6)
    CHARACTER(len=1024):: name

    !Date and time and filenames
    CHARACTER(len=1024):: fname
    CHARACTER(*), PARAMETER:: fpath="./"
    CHARACTER(10):: dat,tim
    INTEGER(8) :: dt(1:8)

    !Required to make ran3 parallelisable
    INTEGER :: iff,inext,inextp,ma(55)
    COMMON /ranvar/ iff,inext,inextp,ma
    !$OMP THREADPRIVATE(/ranvar/)

    OPEN(unit=1,file=fpath//"input",status="old")
    OPEN(unit=2,file=fpath//"mass-dist",status="old")
    OPEN(unit=3,file=fpath//"sys-name",status="old")

    !Read the values into their arrays
    READ(1,*) obsdata
    READ(2,*) properties
    READ(3,*) name
    !Close files after reading
    CLOSE(1)
    CLOSE(2)
    CLOSE(3)

    !Simulation parameters
    n     = 1e9                 !Number of iterations of the main do loop before checking for matches
    kdum  = 625                 !Random seed
    nmatch= 5000                !Number of required matches before terminating the simulation
    nt=omp_get_max_threads()    !Number of threads to use (set to run on all)

    !Constants
    pi  = 4.*ATAN(1.)
    au  = 1.5E11
    year= 365.25*24.*3600.
    msun= 1.989E30
    G   = 6.67E-11
    d2r = pi/180.

    !Counters
    count=0
    loop=0

    !Unpacking the input data (separations and position angles)
    seps    = obsdata(2,:)
    sep_err = obsdata(3,:)
    pa      = obsdata(4,:)*d2r
    pa_err  = obsdata(5,:)*d2r

    !Masses and distance
    M(1)    = properties(1)
    M_err(1)= properties(2)
    M(2)    = properties(3)
    M_err(2)= properties(4)
    d       = properties(5)
    d_err   = properties(6)

    !Converting separations into au
    obs(1) = seps(1)*(1.E-3)*d
    obs(2) = seps(2)*(1.E-3)*d

    !Calculaing the change in position angle between observations
    obs(3) = pa(2)-pa(1)

    !Factor to multiply the errors by (er=1. means use 1 sigma errors,
    !er=2. would use 2 sigma errors, etc.)
    er = 1.

    !Errors on separation and position angles (in au and radians)
    err(1) = obs(1) * SQRT((sep_err(1)/seps(1))**2. + (d_err/d)**2.)*er
    err(2) = obs(2) * SQRT((sep_err(2)/seps(2))**2. + (d_err/d)**2.)*er
    err(3) = (pa_err(1)+pa_err(2))*er

    !Time between observations
    dtobs = obsdata(1,2)-obsdata(1,1)

    !Calculating on sky velocity
    CALL cosine(obs(1),obs(2),obs(3),dmoved)
    vobs = (dmoved*au)/(dtobs*year) 

    !Minimum and maximum semi-major axes
    amin = MIN(obs(1),obs(2))*0.5
    amax = 100.*G*(M(1)+M(2))*msun/(vobs*vobs*au)

    !Creating output file
    WRITE(fname,"(A,A100)") "params-",name
    OPEN(unit=1,file=fpath//trim(fname),status="new")

    !File preamble
        CALL DATE_AND_TIME(DATE=dat,TIME=tim,VALUES=dt)
        WRITE(1,'(A)') "This file contains the results of testing the FOBOS code (Houghton & Goodwin, in prep)"
        WRITE(1,'(A)') "using two astrometric data points on the binary/exoplanetary system: ",trim(name)
        WRITE(1,'(A,I2,A,I2,A,I4,A,I2,A,I2,A,I2)') "This simulation ran on ",dt(3),"/",dt(2),"/",dt(1), & 
        " at ",dt(5),":",dt(6),":",dt(7)
        WRITE(1,'(A,F4.2,A,F4.2,A,F4.2)') "Object masses: m1=",M(1),", m2=",M(2)
        WRITE(1,'(A,F4.2,A)') "Time between observations: ",dtobs," years"
        WRITE(1,'(A,E8.1,A,I2,A)') "Tested ",REAL(n),"fake systems per loop using across ",nt," threads."
        WRITE(1,'(A,I3)') "The random seed value was ",kdum
        WRITE(1,'(A,F3.1)') "The error multiplication factor for this run was:",er
    !End file preamble

    !Random seeds for different threads
    CALL omp_set_num_threads(nt)
    ALLOCATE(rinit(nt))
    do i=1,nt
        rinit(i)=int(ran3(kdum)*1000000)
    end do

    !Main do loop
1   CONTINUE
    !$OMP PARALLEL PRIVATE(my_id,idum,params,match,sep)
    if (loop==0) then
        my_id= omp_get_thread_num()+1
        idum = rinit(my_id)
    end if
    !$OMP DO
    do i=1,n
        CALL fakesystem(amin,amax,dtobs,vobs,M,obs,err,idum,params,match,sep)
        if (match==1) then
            !$OMP CRITICAL
            WRITE(1,*) params
            count=count+1
            !$OMP END CRITICAL
        end if
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    loop = loop + 1
    WRITE(6,*) "Loop:",loop,", Matches:",count
    call cpu_time(t_run)

    if (count<nmatch) goto 1
    
    !End of program admin
    DEALLOCATE(rinit)
    WRITE(1,"(A,F8.2,A)") "Run time =", t_run/nt, "s"
    WRITE(1,*) "Number of loops:",loop, " and final number of matches:",count
    
    !Close files
    CLOSE(1)

    !add some words and things 

END PROGRAM binaries

!----------------------------------------------------------------
!Projects the binary randomly on the sky
SUBROUTINE makeobsbinary(a,ecc,inc,phi,mean,sky)
    IMPLICIT NONE

    REAL, INTENT(in) :: a,ecc,inc,phi,mean
    REAL, INTENT(out):: sky(1:2)
    REAL:: pi,Ee,old,theta,rt,r(1:3),x,y

    pi=4.*ATAN(1.)

    !Calculate eccentric anomaly with Newton-Raphson
    Ee=mean
    old=mean
    do
        old=Ee
        Ee=Ee-(Ee-ecc*SIN(Ee)-mean)/(1.-ecc*COS(Ee))
        if (old-Ee<1.E-6) exit
    end do

    !Calculate true anomaly and true distance
    theta=2.*ATAN(SQRT((1.+ecc)/(1.-ecc))*TAN(0.5*Ee))
    rt = a*(1.-ecc*ecc)/(1.+ecc*COS(theta)) 

    !Get x and y positions
    r(1:3)=0.
    r(1)=rt*COS(theta)
    r(2)=rt*SIN(theta)
    !Orientation around z axis
    x=r(1);y=r(2)
    r(1)=x*COS(phi)-y*SIN(phi)
    r(2)=x*SIN(phi)+y*COS(phi)
    !Modify inclination
    x=r(1);y=r(2)
    r(1)=x
    r(2)=y*COS(inc)
    r(3)=y*SIN(inc)
    !Positions on the sky
    sky(1)=r(1)
    sky(2)=r(3)

END SUBROUTINE makeobsbinary

!----------------------------------------------------------------
SUBROUTINE fakesystem(amin,amax,dtobs,vobs,M,obs,err,idum,params,match,sep)
    IMPLICIT NONE

    REAL, INTENT(in) :: amin,amax,dtobs,vobs,M(1:2),obs(1:3),err(1:3)
    REAL, INTENT(out):: params(1:5),sep(1:3)
    INTEGER, INTENT(in) :: idum
    INTEGER, INTENT(out):: match
    INTEGER:: i
    REAL:: a,ecc,inc,phi,mean1,mean2
    REAL:: sky1(1:2),sky2(1:2),rt
    REAL:: pi,G,au,msun,v,T,dM,f(1:2)
    REAL:: ran3
    EXTERNAL:: ran3

    !Constants
    pi=4.*ATAN(1.)
    G=6.67E-11
    msun=1.989E30
    au=1.496E11

    !Set match-0
    match=0

    !Generate system parameters
    a=amin+(ran3(idum)*amax)
    ecc=ran3(idum)*0.99
    inc=ASIN(ran3(idum))
    phi=ran3(idum)*2.*pi
    mean1=ran3(idum)*2.*pi
    params=(/a,ecc,inc,phi,mean1/)

    !Period
    T=SQRT((a*a*a)/(M(1)+M(2)))

    !Fraction of orbit moved
    dM=(dtobs/T)*2.*pi

    !Project on to the sky
    CALL skyproject(a,ecc,inc,phi,mean1,rt,sky1)
    sep(1)=SQRT(sky1(1)*sky1(1)+sky1(2)*sky1(2))

    !Orbital velocity
    v=G*(M(1)+M(2))*msun*((2./(rt*au))-(1./(a*au)))

    !Throw away systems where v<vobs, true distance<separation
    !or the separation doesn't match observation within errors
    if (v<vobs) goto 10
    if (rt<obs(1)) goto 10
    if (ABS(obs(1)-sep(1))>err(1)) goto 10

    rt=0.
    f=(/1.,-1./)

    !Do loop checking orbit moving both ways
    do i=1,2
        mean2=mean1+f(i)*dM
        CALL skyproject(a,ecc,inc,phi,mean2,rt,sky2)
        if (rt<obs(2)) cycle
        v=G*(M(1)+M(2))*msun*((2./(rt*au))-(1./(a*au)))
        if (v<vobs) goto 10
        sep(2)=SQRT(sky2(1)*sky2(1)+sky2(2)*sky2(2))        
        if (ABS(obs(2)-sep(2))>err(2)) cycle
        sep(3)=ATAN2(sky2(2),sky2(1))-ATAN2(sky1(2),sky1(1))
        if (ABS(ABS(obs(3))-ABS(sep(3)))>err(3)) cycle
        match=1
        exit
    end do

    !End sub with match=0 if no match or match=1 if match

10  CONTINUE

CONTAINS
SUBROUTINE skyproject(a,ecc,inc,phi,mean,rt,sky)
    IMPLICIT NONE

    REAL, INTENT(in) :: a,ecc,inc,phi,mean
    REAL, INTENT(out):: rt,sky(1:2)
    REAL:: Ee,old,theta,r(1:3),x,y

    !Calculate eccentric anomaly with Newton-Raphson
    Ee=mean
    old=mean
    do
        old=Ee
        Ee=Ee-(Ee-ecc*SIN(Ee)-mean)/(1.-ecc*COS(Ee))
        if (old-Ee<1.e-6) exit
    end do

    !Calculate true anomaly and true distance
    theta=2.*ATAN(SQRT((1.+ecc)/(1.-ecc))*TAN(0.5*Ee))
    rt = a*(1.-ecc*ecc)/(1.+ecc*COS(theta)) 

    !Get x and y positions
    r(1:3)=0.
    r(1)=rt*COS(theta)
    r(2)=rt*SIN(theta)
    !Orientation around z axis
    x=r(1);y=r(2)
    r(1)=x*COS(phi)-y*SIN(phi)
    r(2)=x*SIN(phi)+y*COS(phi)
    !Modify inclination
    x=r(1);y=r(2)
    r(1)=x
    r(2)=y*COS(inc)
    r(3)=y*SIN(inc)
    !Positions on the sky
    sky(1)=r(1)
    sky(2)=r(3)

END SUBROUTINE skyproject    

END SUBROUTINE fakesystem

!----------------------------------------------------------------
SUBROUTINE cosine(a,b,theta,c)
    IMPLICIT NONE

    REAL, INTENT(in) :: a,b,theta
    REAL, INTENT(out):: c

    c=SQRT(a*a + b*b - 2.*a*b*COS(theta))

END SUBROUTINE cosine

!======================================================
!--------------------RAN3------------------------------
!Random number generator from Numerical Recipes in Fortran
REAL FUNCTION ran3(idum)
    INTEGER :: idum
    INTEGER :: MBIG,MSEED,MZ
    REAL :: FAC
    PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
    INTEGER :: i,iff,ii,inext,inextp,k
    INTEGER :: mj,mk,ma(55)
    !SAVE iff,inext,inextp,ma
    COMMON /ranvar/ iff,inext,inextp,ma
    !$OMP THREADPRIVATE(/ranvar/)
    DATA iff /0/
    IF(idum.LT.0.OR.iff.EQ.0)THEN
      iff=1
      mj=MSEED-iabs(idum)
      mj=MOD(mj,MBIG)
      ma(55)=mj
      mk=1
      DO 11 i=1,54
        ii=MOD(21*i,55)
        ma(ii)=mk
        mk=mj-mk
        IF(mk.LT.MZ)mk=mk+MBIG
        mj=ma(ii)
    11      CONTINUE
      DO 13 k=1,4
        DO 12 i=1,55
          ma(i)=ma(i)-ma(1+MOD(i+30,55))
          IF(ma(i).LT.MZ)ma(i)=ma(i)+MBIG
    12        CONTINUE
    13      CONTINUE
      inext=0
      inextp=31
      idum=1
    ENDIF
    inext=inext+1
    IF(inext.EQ.56)inext=1
    inextp=inextp+1
    IF(inextp.EQ.56)inextp=1
    mj=ma(inext)-ma(inextp)
    IF(mj.LT.MZ)mj=mj+MBIG
    ma(inext)=mj
    ran3=mj*FAC
    RETURN
END FUNCTION ran3

FUNCTION ran_normal(mu,sigma,idum)
    REAL:: ran3,ran_normal
    REAL:: pi
    EXTERNAL:: ran3
    pi=4.*ATAN(1.)
    ran_normal=(1./(sigma*SQRT(2*pi))*EXP(-0.5*((ran3(idum)-mu)/sigma)**2.))
    return
END FUNCTION ran_normal



