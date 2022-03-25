! FOBOS-triples.f90 takes two epochs of observation of a three body system and fits
! the orbital parameters (semi-major axis, eccentricity, inclination, orientation
! and mean anomaly) for both companions. The code is currently set up to read from a file called "input", which contains two
! separations and their associated errors (in milliarcseconds) for the secondary (columns 2 and 3) 
! and the tertiary (columns 4 and 5). Each observation also contains two position angles for each
! star with errors (in degrees, columns 5 and 6 for the secondary and 7 and 8 for the tertiary). 
! The input file should have the format as shown below.

! t1 (yrs)  sep1_in (mas)  sep_err1_in (mas)  sep1_out (mas)  sep_err1_out (mas) pa1_in (deg)   pa_err1_in (deg)   pa1_out (deg)   pa_err1_out (deg)
! t2 (yrs)  sep2_in (mas)  sep_err2_in (mas)  sep2_out (mas)  sep_err2_out (mas) pa2_in (deg)   pa_err2_in (deg)   pa2_out (deg)   pa_err2_out (deg)

! The target name, which affects the name of the output files, is given in a text file
! called "sys-name". The masses and distance to the objects are contained in the file
! "mass-dist", which has the format

!M1 M1_err  M2  M2_err  M3  M3_err   d   d_err

! All masses will be given in solar masses and distances in parsecs.

PROGRAM triples
    IMPLICIT NONE
    INCLUDE "omp_lib.h"

    INTEGER:: n,i,kdum,idum,match,nmatch,nt,my_id
    INTEGER:: loop,loopten,count
    REAL:: pi,au,year,msun,G,d2r
    REAL:: M(1:3),M_err(1:3),d,d_err,er
    REAL:: dtobs,seps(1:4),seps_err(1:4),pa(1:4),pa_err(1:4)
    REAL:: vobs(1:2),amin(1:2),amax(1:2),dmoved(1:2)
    REAL:: obs(1:7),err(1:7)
    REAL:: params(1:10)
    REAL:: texit,t_run
    REAL:: ran3
    INTEGER:: pro,retro,pr
    INTEGER,ALLOCATABLE:: rinit(:)
    !Parameters read in from files
    REAL:: obsdata(1:9,1:2)
    REAL:: properties(1:8)
    CHARACTER(len=1024)::name

    !Date and time and filenames
    CHARACTER(len=1024):: fname1
    CHARACTER(*), PARAMETER:: fpath="./"
    CHARACTER(10):: dat,tim
    INTEGER(8) :: dt(1:8)

    !Required to make ran3 parallelisable
    INTEGER :: iff,inext,inextp,ma(55)
    COMMON /ranvar/ iff,inext,inextp,ma
    !$OMP THREADPRIVATE(/ranvar/)

    !Open initial conditions files
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
    n      = 1e8 
    kdum   = 625
    nmatch = 10000
    nt     = omp_get_max_threads()
    texit  = 3600.

    !Constants
    pi  = 4.*ATAN(1.)
    au  = 1.5E11
    year= 3.1558E7
    msun= 1.989E30
    G   = 6.67E-11
    d2r = pi/180.

    !Might not need a bunch of these
    count  = 0
    pro    = 0
    retro  = 0
    loop   = 0
    loopten= 0

    !Unpacking the input data
    !These arrays do not exist yet
    seps(1:2)     = obsdata(2,:)
    seps_err(1:2) = obsdata(3,:)
    seps(3:4)     = obsdata(4,:)
    seps_err(3:4) = obsdata(5,:)
    pa(1:2)       = obsdata(6,:)*d2r
    pa_err(1:2)   = obsdata(7,:)*d2r
    pa(3:4)       = obsdata(8,:)*d2r
    pa_err(3:4)   = obsdata(9,:)*d2r

    !Masses and distance
    M(1)    = properties(1)
    M_err(1)= properties(2)
    M(2)    = properties(3)
    M_err(2)= properties(4)
    M(3)    = properties(5)
    M_err(3)= properties(6)
    d       = properties(7)
    d_err   = properties(8)

    !Entering data
    obs(1) = seps(1)*(1.E-3)*d
    obs(2) = seps(2)*(1.E-3)*d
    obs(3) = pa(2)-pa(1)
    obs(4) = seps(3)*(1.E-3)*d
    obs(5) = seps(4)*(1.E-3)*d
    obs(6) = pa(4)-pa(3)

    !Error
    er = 1.

    !Errors
    err(1) = obs(1) * SQRT((seps_err(1)/seps(1))**2. + (d_err/d)**2.)*er
    err(2) = obs(2) * SQRT((seps_err(2)/seps(2))**2. + (d_err/d)**2.)*er
    err(3) = (pa_err(1)+pa_err(2))*er
    err(4) = obs(4) * SQRT((seps_err(3)/seps(3))**2. + (d_err/d)**2.)*er
    err(5) = obs(5) * SQRT((seps_err(4)/seps(4))**2. + (d_err/d)**2.)*er
    err(6) = (pa_err(3)+pa_err(4))*er

    !Calculating the relative positions of the two companions
    CALL cosine(obs(1),obs(3),ABS(pa(1)-pa(3)),obs(7))
    err(7) = (pa(1)+pa(3))*er

    !Calculating the distance moved by each star using the cosine rule
    !Function cosine2 calculates the distances with errors
    CALL cosine2(obs(1),0.,obs(2),0.,pa(2),pa_err(2),pa(1),pa_err(1),dmoved(1))
    CALL cosine2(obs(4),0.,obs(5),0.,pa(4),pa_err(4),pa(3),pa_err(3),dmoved(2))

    !Time between observations
    dtobs=obsdata(1,2)-obsdata(1,1)

    !Separation minimum value
    amin(1) = MIN(obs(1),obs(2))*0.5
    amin(2) = MIN(obs(4),obs(5))*0.5

    !Upper limits for a1 and a2
    vobs(1) = (dmoved(1)*au)/(dtobs*year)
    vobs(2) = (dmoved(2)*au)/(dtobs*year)
    amax(1) = 200.*G*(M(1)+M(2))*msun/(vobs(1)*vobs(1)*au)
    amax(2) = 200.*G*(M(1)+M(2)+M(3))*msun/(vobs(2)*vobs(2)*au)

    !Creating output files to write to
    WRITE(fname1,"(A,A100)") "params-",name
    OPEN(unit=1,file=fpath//trim(fname1),status="new")

    !Metadata
    CALL DATE_AND_TIME(DATE=dat,TIME=tim,VALUES=dt)
    do i=1,3
        WRITE(i,'(A)') "This file contains the results of testing the FOBOS-triples code (Houghton & Goodwin, in prep"
        WRITE(i,'(A)') "on the binary/exoplanetary system: ",trim(name)
        WRITE(i,'(A,I2,A,I2,A,I4,A,I2,A,I2,A,I2)') "This simulation ran on  ",dt(3),"/",dt(2),"/",dt(1), & 
                "  at  ",dt(5),":",dt(6),":",dt(7)
        WRITE(i,'(A,F4.2,A,F4.2,A,F4.2)') "Masses: m1=",M(1),", m2=",M(2),", and m3=",M(3)
        WRITE(i,'(A,F4.2,A)') "Time between observations:",dtobs,"years"
        WRITE(i,'(A,E8.1,A,F4.1,A)') "Tested ", REAL(n),"fake systems per loop with a match precision of ",err*100.,"%"
        WRITE(i,'(A,I3,A,I2,A)') "The random seed value was ",kdum, " and it ran on ",nt," threads."
    end do

    !Different random number seeds for threads
    CALL omp_set_num_threads(nt)
    ALLOCATE(rinit(nt))
    do i=1,nt
        rinit(i)=int(ran3(kdum)*1000000)
    end do

    !MAIN DO LOOP
1   CONTINUE
    !$OMP PARALLEL PRIVATE(my_id,idum,params,match,sep,coords,pr)
    if (loop==0) then
        my_id = omp_get_thread_num()+1
        idum  = rinit(my_id)
    end if
    !$OMP DO
    do i=1,n
        CALL fakesystem(idum,amin,amax,M,dtobs,obs,vobs,err,params,match,pr)
        if (match==1) then
            !$OMP CRITICAL
            count=count+1
            WRITE(1,*) params,pr
            !$OMP END CRITICAL
        end if
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    loop = loop + 1
    loopten=loopten+1

    WRITE(6,*) "Loop:",loop,", Matches:",count

    !Exit the simulation with no matches if it takes too long
    CALL cpu_time(t_run)

    if (t_run<texit) then
        if (count<nmatch) goto 1
    end if

    !Admin
    DEALLOCATE(rinit)
    WRITE(1,*) "Number of loops:",loop, " and final number of matches:",count
    WRITE(1,"(A,F8.2,A)") "Run time =", t_run/nt, "s"
    
    CLOSE(1)

END PROGRAM triples


SUBROUTINE fakesystem(idum,amin,amax,M,dtobs,obs,vobs,err,params,match,pr)
    IMPLICIT NONE

    REAL, INTENT(in) :: amin(1:2),amax(1:2),M(1:3),dtobs,obs(1:7),vobs(1:2),err(1:7)
    REAL, INTENT(out):: params(1:10)
    INTEGER, INTENT(in) :: idum
    INTEGER, INTENT(out):: match,pr
    INTEGER:: i,j,n,match1,match2,match3,is(1:2)
    REAL:: a1,e1,inc1,phi1,mean1i,mean1f
    REAL:: a2,e2,inc2,phi2,mean2i,mean2f
    REAL:: pi,msun,au,G,T1,T2,dM1,dM2,rt,a2min,v
    REAL:: sky1i(1:2),sky1f(1:2),sky2i(1:2),sky2f(1:2)
    REAL:: stab,qout,sep(1:7)
    REAL:: ran3
    EXTERNAL:: ran3

    !Constants
    pi=4.*ATAN(1.)
    G=6.67E-11
    msun=1.989E30
    au=1.496E11

    n=1000
    match=0
    match1=0
    match2=0

    !Generate inner binary parameters
    a1=amin(1)+(ran3(idum)*amax(1))
    e1=ran3(idum)*0.99
    inc1=ASIN(2.*ran3(idum)-1.) !HERE 2.*ran3(idum)-1.
    phi1=ran3(idum)*2.*pi
    mean1i=ran3(idum)*2.*pi

    !Calculate period and movement
    T1=SQRT((a1*a1*a1)/(M(1)+M(2)))
    dM1=(dtobs/T1)*2.*pi

    !Project the binary on the sky
    CALL skyproject(a1,e1,inc1,phi1,mean1i,rt,sky1i)
    !True distance can't be smaller than separation
    if (rt<obs(1)) goto 10 
    !Velocity check
    v=G*(M(1)+M(2))*msun*((2./(rt*au))-(1./(a1*au)))
    if (v<vobs(1)) goto 10
    !Get first secondary separation
    sep(1)=SQRT(sky1i(1)*sky1i(1)+sky1i(2)*sky1i(2))
    if (ABS(obs(1)-sep(1))>err(1)) goto 10

    !Move the binary in both directions
    match1=0
    do i=1,2
        mean1f=mean1i+(i*2.-3.)*dM1
        CALL skyproject(a1,e1,inc1,phi1,mean1f,rt,sky1f)
        if (rt<obs(2)) cycle
        !Velocity check
        v=G*(M(1)+M(2))*msun*((2./(rt*au))-(1./(a1*au)))
        if (v<vobs(1)) goto 10
        !Check second secondary separation
        sep(2)=SQRT(sky1f(1)*sky1f(1)+sky1f(2)*sky1f(2))
        if (ABS(obs(2)-sep(2))>err(2)) cycle
        !Calculate angle moved
        sep(3)=ATAN2(sky1f(2),sky1f(1))-ATAN2(sky1i(2),sky1i(1))
        if (ABS(obs(3)-sep(3))>err(3)) cycle
        match1=1
        is(1)=i*2-3
        exit
    end do

    if (match1==0) goto 10
    

    !--------------------------------------------------------------
    !Generate tertiary parameters

    !If taken out of j loop, switch cycles back to gotos
    do j=1,n

        !Mardling-Aarseth stability 
        e2=ran3(idum)*0.99
        inc2=ASIN(2.*ran3(idum)-1.) !HERE

        !Checking if the masses are too small for the stability condition
        if (M(1)>0.08.AND.M(2)>0.08.AND.M(3)>0.08) then
            qout=M(3)/(M(1)+M(2))
            stab=a1*(2.8/(1.-e2)) * (1-(0.3*ABS(inc2-inc1)/pi)) * ((1.+qout)*(1+e2)/SQRT(1-e2))**0.4
            if (amin(2)<stab) then
                a2min=stab
            else
                a2min=amin(2)
            end if
            if (a2min>=amax(2)) cycle!goto 10

        else
            a2min=amin(2)
        end if

        a2=a2min+(ran3(idum)*(amax(2)-a2min))
        phi2=ran3(idum)*2.*pi
        mean2i=ran3(idum)*2.*pi

        !Calculate period and movement
        !Assumption that tertiary is far enough away that M=M(1:3)
        T2=SQRT((a2*a2*a2)/(M(1)+M(2)+M(3)))
        dM2=(dtobs/T2)*2.*pi
        
        !Project the binary on the sky
        CALL skyproject(a2,e2,inc2,phi2,mean2i,rt,sky2i)
        if (rt<obs(4)) cycle
        ! Velocity check
        v=G*(M(1)+M(2)+M(3))*msun*((2./(rt*au))-(1./(a2*au)))
        if (v<vobs(2)) cycle
        ! Get first tertiary separation
        sep(4)=SQRT(sky2i(1)*sky2i(1)+sky2i(2)*sky2i(2))
        if (ABS(obs(4)-sep(4))>err(4)) cycle !Bottleneck

        !Move the binary in both directions and see if it matches
        match2=0
        do i=1,2
            mean2f=mean2i+(i*2.-3.)*dM2
            CALL skyproject(a2,e2,inc2,phi2,mean2f,rt,sky2f)
            if (rt<obs(5)) cycle
            !Velocity check
            v=G*(M(1)+M(2)+M(3))*msun*((2./(rt*au))-(1./(a2*au)))
            if (v<vobs(2)) cycle
            !Check second tertiary separation
            sep(5)=SQRT(sky2f(1)*sky2f(1)+sky2f(2)*sky2f(2))
            if (ABS(obs(5)-sep(5))>err(5)) cycle
            !Calculate angle moved
            sep(6)=ATAN2(sky2f(2),sky2f(1))-ATAN2(sky2i(2),sky2i(1))
            if (ABS(obs(6)-sep(6))>err(6)) cycle
            match2=1
            is(2)=i*2-3
            exit
        end do

        if (match2==1) exit
    end do

    !Additional angle relating secondary and tertiary positions
    !Doesn't matter which way this is defined as long as it's consistent with observations
    ! match3=0
    ! sep(7)=ATAN2(sky2i(2),sky2i(1))-ATAN2(sky1i(2),sky1i(1))
    ! if (ABS(obs(7)-sep(7))>err(7)) goto 10
    match3=1

    !Calculating separations between secondary and tertiary
    ! match3=0
    ! sep(7)=SQRT((sky1i(1)-sky2i(1))**2+(sky1i(2)-sky2i(2))**2)
    ! sep(8)=SQRT((sky1f(1)-sky2f(1))**2+(sky1f(2)-sky2f(2))**2)
    ! if (ABS(obs(7)-sep(7))>err(7)) goto 10
    ! if (ABS(obs(8)-sep(8))>err(8)) goto 10
    ! match3=1

    if (match1==1.AND.match2==1.AND.match3==1) then
        match=1
        if (is(1)+is(2) == 0) then
            pr=1
        else
            pr=0
        end if
        params=(/a1,a2,e1,e2,inc1,inc2,phi1,phi2,mean1i,mean2i/)
    end if

10 CONTINUE 
CONTAINS
SUBROUTINE skyproject(a,ecc,inc,phi,mean,rt,sky)
    IMPLICIT NONE

    REAL, INTENT(in) :: a,ecc,inc,phi,mean
    REAL, INTENT(out):: rt,sky(1:2)
    REAL:: Ee,old,theta,r(1:3),x,y

    !Numerically solve for eccentric anomaly
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

!----------------------------------------------------------------
SUBROUTINE cosine2(a,a_err,b,b_err,t1,t1_err,t2,t2_err,c)
    IMPLICIT NONE

    REAL, INTENT(in) :: a,a_err,b,b_err,t1,t2,t1_err,t2_err
    REAL, INTENT(out):: c
    REAL:: theta

    theta=(t2-t2_err)-(t1+t1_err)
    c=SQRT((a-a_err)**2 + (b-b_err)**2 - 2.*(a-a_err)*(b-b_err)*COS(theta))

END SUBROUTINE cosine2

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