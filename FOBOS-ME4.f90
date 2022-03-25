! Modified 14/01/22
    ! FOBOS-ME4.f90 takes four epochs of observation of a two body system and fits
    ! the orbital parameters (semi-major axis, eccentricity, inclination, orientation
    ! and mean anomaly). The code is currently set up to read in three separations (in milliarcseconds,
    ! column 2) and three position angles (in degrees, column 4), plus their associated respective errors, 
    ! (columns 3 and 5) from a text file called "input". The input file should have the format as shown below.

    ! 0.0000     sep1 (mas)   sep1_err (mas)   pa1 (deg)   pa1_err (deg)
    ! t1 (yrs)   sep2 (mas)   sep2_err (mas)   pa2 (deg)   pa2_err (deg)
    ! t2 (yrs)   sep3 (mas)   sep3_err (mas)   pa3 (deg)   pa3_err (deg)

    ! t1 is the time between the first and second epochs of observation (in years) and 
    ! t2 is the time between the second and third epochs of observation (in years)

    ! The target name, which affects the name of the output files, is given in a text file
    ! called "sys-name". The masses and distance to the objects are contained in the file
    ! "mass-dist", which has the format

    ! M1 M1_err  M2  M2_err  d   d_err

    ! All masses will be given in solar masses and distances in parsecs
! End of information

PROGRAM binariesme
    IMPLICIT NONE
    INCLUDE "omp_lib.h"

    INTEGER:: n,i,idum,kdum,match,nmatch
    INTEGER:: nt,my_id,count,loop
    REAL:: M(1:2),M_err(1:2),dtobs1,dtobs2,dtobs3,d,d_err
    REAL:: obs(1:7),prec(1:7),params(1:5),coords(1:6),sep(1:7)
    REAL:: dts(1:4),seps(1:4),sep_err(1:4),pa(1:4),pa_err(1:4)
    REAL:: dmoved1,dmoved2,dmoved3
    REAL:: pi,au,year,msun,G,d2r
    REAL:: amin,amax,vobs,vobs1,vobs2,vobs3
    REAL:: t_run,texit,er
    REAL:: ran3
    INTEGER,ALLOCATABLE:: rinit(:)
    !Params
    REAL:: obsdata(1:5,1:4)
    REAL:: properties(1:6)
    CHARACTER(len=1024):: name

    !Date and time and filenames
    CHARACTER(len=1024):: fname1
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
    n=1e9           !Number of iterations of the main do loop before checking for matches
    kdum=625        !Random seed
    nmatch=5000     !Number of required matches before terminating the simulation
    texit=3600.     !Time after which to exit the simulation if too few matches have been found (seconds)

    nt=omp_get_max_threads()    !Set number of threads to run on as maximum

    !Constants
    pi  = 4.*ATAN(1.)
    au  = 1.5E11
    year= 365.25*24.*3600.
    msun= 1.989E30
    G   = 6.67E-11
    d2r = (0.5*pi)/90

    !Counters
    count=0
    loop =0

    !Unpacking the input data
    dts     = obsdata(1,:)
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

    !Assign separations and angles to points in the obs array
    obs(1) = seps(1)*(1.E-3)*d
    obs(2) = seps(2)*(1.E-3)*d
    obs(3) = pa(2)-pa(1)
    obs(4) = seps(3)*(1.E-3)*d
    obs(5) = pa(3)-pa(2)
    obs(6) = seps(4)*(1.E-3)*d
    obs(7) = pa(4)-pa(3)

    !Factor to multiply the errors by (er=1. means use 1 sigma errors,
    !er=2. would use 2 sigma errors, etc.)
    er=1.

    !Observational errors
    prec(1) = obs(1) * SQRT((sep_err(1)/seps(1))**2. + (d_err/d)**2.)*er
    prec(2) = obs(2) * SQRT((sep_err(2)/seps(2))**2. + (d_err/d)**2.)*er
    prec(3) = (pa_err(1)+pa_err(2))*er
    prec(4) = obs(4) * SQRT((sep_err(3)/seps(3))**2. + (d_err/d)**2.)*er
    prec(5) = (pa_err(2)+pa_err(3))*er
    prec(6) = obs(6) * SQRT((sep_err(4)/seps(4))**2. + (d_err/d)**2.)*er
    prec(7) = (pa_err(3)+pa_err(4))*er

    !Time between observations
    dtobs1 = dts(2)-dts(1)
    dtobs2 = dts(3)-dts(2)
    dtobs3 = dts(4)-dts(3)

    !Calculating the on sky velocities
    dmoved1 = SQRT(obs(1)**2+obs(2)**2-(2.*obs(1)*obs(2)*COS(obs(3))))
    dmoved2 = SQRT(obs(2)**2+obs(4)**2-(2.*obs(2)*obs(4)*COS(obs(5))))
    dmoved3 = SQRT(obs(4)**2+obs(6)**2-(2.*obs(4)*obs(6)*COS(obs(7))))

    !Minimum semi-major axis
    amin = MIN(obs(1),obs(2),obs(4),obs(6))*0.5
    
    !Observed velocities
    vobs1 = (dmoved1*au)/(dtobs1*year)
    vobs2 = (dmoved2*au)/(dtobs2*year)
    vobs3 = (dmoved3*au)/(dtobs3*year)

    !Maximum semi-major axis
    vobs = MIN(vobs1,vobs2,vobs3)
    amax = 20.*G*(M(1)+M(2))*msun/(vobs*vobs*au)

    !Getting the output file set up
    WRITE(fname1,"(A7,A9)") "params-",name
    OPEN(unit=1,file=fpath//trim(fname1),status="new")

    !Params file preamble
        CALL DATE_AND_TIME(DATE=dat,TIME=tim,VALUES=dt)
        WRITE(1,'(A)') "The results contained in this file are for testing FOBOS-ME: Few Observation Binary Orbit Solver"
        WRITE(1,'(A,A9)') "(Multiple Epoch) on the binary/exoplanetary system ",trim(name)
        WRITE(1,'(A,I2,A,I2,A,I4,A,I2,A,I2,A,I2)') "This simulation ran on  ",dt(3),"/",dt(2),"/",dt(1), & 
        "  at  ",dt(5),":",dt(6),":",dt(7)
        WRITE(1,'(A,E8.1,A)') "Tested ", REAL(n),"fake systems per loop"
        WRITE(1,'(A,I3,A,I1,A)') "The random seed value was ",kdum, " and it ran on ",nt," threads."
        WRITE(1,'(A,F4.2,A,F4.2,A,F4.2)') "Masses: m1=",M(1),", m2=",M(2)
    !End preamble

    !Random seeds for different threads
    CALL omp_set_num_threads(nt)
    ALLOCATE(rinit(nt))
    do i=1,nt
        rinit(i)=int(ran3(kdum)*1000000)
    end do

    !Main do loop
1   continue
    !$OMP PARALLEL PRIVATE(my_id,idum,params,match,sep)
    if (loop == 0) then
        my_id= omp_get_thread_num()+1
        idum = rinit(my_id)
    end if
    !$OMP DO
    do i=1,n
        CALL fakesystem(amin,amax,dtobs1,dtobs2,dtobs3,vobs,M,obs,prec,idum,params,match,coords,sep)
        if (match == 1) then
            !$OMP CRITICAL
            count=count+1
            WRITE(1,*) params
            WRITE(2,*) sep
            WRITE(3,*) coords
            !$OMP END CRITICAL
        end if
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    !Write to screen with how many loops of n interations have been found
    !and the number of matches so far 
    loop = loop + 1
    WRITE(6,*) "Looped:",loop,"Number of matches:",count
    call cpu_time(t_run)

    !Exit the simulation with no matches if it takes too long
    if (t_run/nt < texit) then
       if (count<nmatch) goto 1
    end if

    !If simulation stops because of the time write to screen
    if (count == 0) then
        WRITE(6,*) "Timeout"
    end if
    
    !End of program admin
    DEALLOCATE(rinit)
    WRITE(1,*) "Number of loops:",loop," and number of matches:",count
    WRITE(1,"(A,F8.2,A)") "Run time =", t_run/nt, "s"
    WRITE(6,"(A,F8.2,A)") "Run time =", t_run/nt, "s"

    CLOSE(1)

END PROGRAM binariesme

!----------------------------------------------------------------
!Projects the binary randomly on the sky
SUBROUTINE makeobsbinary(a,ecc,inc,phi,mean,sky)
    IMPLICIT NONE

    REAL, INTENT(in) :: a,ecc,inc,phi,mean
    REAL, INTENT(out):: sky(1:2)
    REAL:: pi,Ee,old,theta,rt,r(1:3),x,y

    pi=4.*ATAN(1.)

    !Calculate eccentric anomoly
    Ee=mean
    old=mean
    do
        old=Ee
        Ee=Ee-(Ee-ecc*SIN(Ee)-mean)/(1.-ecc*COS(Ee))
        if (old-Ee<1.E-6) exit
    end do

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

!-----------------------------------------------------------------
SUBROUTINE fakesystem(amin,amax,dtobs1,dtobs2,dtobs3,vobs,M,obs,prec,idum,params,match,coords,sep)
    IMPLICIT NONE

    REAL, INTENT(in) :: amin,amax,dtobs1,dtobs2,dtobs3,vobs,M(1:2),obs(1:7),prec(1:7)
    REAL, INTENT(out):: params(1:5),coords(1:6),sep(1:7)
    INTEGER, INTENT(in) :: idum
    INTEGER, INTENT(out):: match
    INTEGER:: i,itemp,match2,match3,match4
    REAL:: a,ecc,inc,phi
    REAL:: mean1,mean2,mean3,mean4
    REAL:: dM1,dM2,dM3
    REAL:: sky1(1:2),sky2(1:2),sky3(1:2),sky4(1:4)
    REAL:: pi,G,msun,au,T,f(1:2),v,rt
    
    REAL:: ran3
    EXTERNAL:: ran3

    !Constants and zeroing array
    pi=4.*ATAN(1.)
    G=6.67E-11
    msun=1.989E30
    au=1.496E11
    match=0
    match2=0
    match3=0
    match4=0
    itemp=0

    !Generate system parameters
    a = amin+(ran3(idum)*amax)
    ecc = ran3(idum)*0.99
    inc = ASIN(ran3(idum))
    phi = ran3(idum)*2.*pi
    mean1 = ran3(idum)*2.*pi

    !Period
    T  = SQRT((a*a*a)/(M(1)+M(2)))
    dM1 = (dtobs1/T)*2.*pi
    dM2 = (dtobs2/T)*2.*pi
    dM3 = (dtobs3/T)*2.*pi

    params=(/a,ecc,inc,phi,mean1/)

    !Projecting and testing first point
    CALL skyproject(a,ecc,inc,phi,mean1,rt,sky1)
    sep(1)=SQRT(sky1(1)*sky1(1)+sky1(2)*sky1(2))
    if (rt < obs(1)) goto 10
    if (ABS(obs(1)-sep(1))>prec(1)) goto 10

    rt=0.
    f=(/1.,-1./)

    !Checking second point in prograte and retrograde
    do i=1,2
        mean2 = mean1 + f(i)*dM1
        CALL skyproject(a,ecc,inc,phi,mean2,rt,sky2)
        if (rt < obs(2)) cycle
        !Velocity check
        v=G*(M(1)+M(2))*msun*((2./(rt*au))-(1./(a*au)))
        if (v < vobs) cycle
        !Separation
        sep(2)=SQRT(sky2(1)*sky2(1)+sky2(2)*sky2(2))
        if (ABS(obs(2)-sep(2)) > prec(2)) cycle
        !Angle
        sep(3)=ATAN2(sky2(2),sky2(1))-ATAN2(sky1(2),sky1(1))
        !if (ABS(ABS(obs(3))-ABS(sep(3))) > prec(3)) cycle
        if (obs(3)>(sep(3)-ABS(prec(3))).AND.obs(3)<(sep(3)+ABS(prec(3)))) then
            continue
        else 
            cycle
        end if
        match2=1
        itemp=i
        exit
    end do

    if (match2==0) goto 10

    !Adding the third data point (same direction)----------------------
    mean3 = mean2 + f(itemp)*dM2
    CALL skyproject(a,ecc,inc,phi,mean3,rt,sky3)
    if (rt < obs(4)) goto 10
    !Velocity check
    v = G*(M(1)+M(2))*msun*((2./(rt*au))-(1./(a*au)))
    if (v < vobs) goto 10
    !Separation
    sep(4)=SQRT(sky3(1)*sky3(1)+sky3(2)*sky3(2))
    if (ABS(obs(4)-sep(4)) > prec(4)) goto 10
    !Angle
    sep(5)=ATAN2(sky3(2),sky3(1))-ATAN2(sky2(2),sky2(1))
    !if (ABS(ABS(obs(5))-ABS(sep(5))) > prec(5)) goto 10
    if (obs(5)>(sep(5)-ABS(prec(5))).AND.obs(5)<(sep(5)+ABS(prec(5)))) then
        continue
        match3=1
    else 
        goto 10
    end if
    if (match3==0) goto 10

    !Adding the fourth data point (same direction)----------------------
    mean4 = mean3 + f(itemp)*dM3
    CALL skyproject(a,ecc,inc,phi,mean4,rt,sky4)
    if (rt < obs(6)) goto 10
    !Velocity check
    v = G*(M(1)+M(2))*msun*((2./(rt*au))-(1./(a*au)))
    if (v < vobs) goto 10
    !Separation
    sep(6)=SQRT(sky4(1)*sky4(1)+sky4(2)*sky4(2))
    if (ABS(obs(6)-sep(6)) > prec(6)) goto 10
    !Angle
    sep(7)=ATAN2(sky4(2),sky4(1))-ATAN2(sky3(2),sky3(1))
    !if (ABS(ABS(obs(5))-ABS(sep(5))) > prec(5)) goto 10
    if (obs(7)>(sep(7)-ABS(prec(7))).AND.obs(7)<(sep(7)+ABS(prec(7)))) then
        continue
        match4=1
    else 
        goto 10
    end if
    if (match4==0) goto 10

    coords=(/sky1,sky2,sky3/)
    match=1

10  CONTINUE

CONTAINS
SUBROUTINE skyproject(a,ecc,inc,phi,mean,rt,sky)
    IMPLICIT NONE

    REAL, INTENT(in) :: a,ecc,inc,phi,mean
    REAL, INTENT(out):: rt,sky(1:2)
    REAL:: Ee,old,theta,r(1:3),x,y

    Ee=mean
    old=mean
    do
        old=Ee
        Ee=Ee-(Ee-ecc*SIN(Ee)-mean)/(1.-ecc*COS(Ee))
        if (old-Ee<1.e-6) exit
    end do

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

!-----------------
SUBROUTINE cosine(a,b,thetai,thetaf,c)
    IMPLICIT NONE

    REAL, INTENT(in) :: a,b,thetai,thetaf
    REAL, INTENT(out):: c
    REAL:: pi,d2r,dtheta

    !Convert to radians
    pi=4.*ATAN(1.)
    dtheta=ABS(thetai-thetaf)
    d2r=(2.*pi)/360.
    c=SQRT(a*a + b*b - 2.*a*b*COS(dtheta*d2r))

END SUBROUTINE cosine


!======================================================
!--------------------RAN3------------------------------
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