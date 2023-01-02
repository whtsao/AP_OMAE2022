!===============================================================================
!	SOLVE 2D WATER WAVE PROBLEM WITH ONE SIDE WAVEMAKER AND ONESIDE FREE BY BEM
!	WEN-HAUI TSAO 2021 LSU
!===============================================================================
      PROGRAM PMTLD_IN_WAVE
      IMPLICIT NONE
	  INTEGER I,J,IR,IL,NT,ITR,JTR,NGA,NTIM,NFIELD,NPL,NFL,NBT,WTYP,WMKTYP,NNODE,NELM,ICON,NITER,NWG,NDOF,E1LOC,CTYPE
	  INTEGER,ALLOCATABLE::NELEM(:),ME(:),NS(:),BUNTYP(:),LN(:,:)
	  REAL*8 PI,TIME,DELTTIME,P_ATM,WC,E1,E2,ETOL,EA(3),EATOL,D_OUT,GRAV
      REAL*8 T_START,T_INITIAL,T_STEP,T_END
	  REAL*8,ALLOCATABLE::WT(:),RT(:),SHA1(:),SHA2(:),SH(:,:)
      
!---CLAIM FOR WAVE
      REAL*8 DEP,MU,WIDTH,THO,AMP,OMEGA,PSI,DIS,VEL,ACC,WGX(10),WGY(10)
	  REAL*8,ALLOCATABLE::COOR(:,:),SIDE_L(:)
	  REAL*8,ALLOCATABLE::NODE(:,:),NORM(:,:),JCB(:),LENG(:),PHI(:),PPHI(:),PHIT(:),PPHIT(:)
	  REAL*8,ALLOCATABLE::KER1(:,:),KER2(:,:),DP(:,:),DPDS(:),D2PDT(:),D2P(:,:)
      REAL*8,ALLOCATABLE::PR(:),DPDT(:),ACCMO(:),PHIT_TEMP(:)
      REAL*8,ALLOCATABLE::G1(:,:),H1(:,:),EYE(:)
      REAL*8,ALLOCATABLE::NODE_PC(:,:),PHI_PC(:)
	  REAL*8,ALLOCATABLE::DI(:),VE(:),AC(:)
!---CLAIM FOR PMTLD
      INTEGER TLD_YN,NNODE2,NELM2,NELEM2(4),ME2(4),NS2(4),BUNTYP2(4),ICON2
      INTEGER,ALLOCATABLE::LN2(:,:)
      REAL*8 KE,PE,TE,WORK,DTDT,DVDT,DEDT,DWDT
      REAL*8 COOR2(4,2),FOR2(3),MU1,MU2,POR,TL,TH,TWN,TLDLOC(2),WG1_TLD,WG2_TLD
      REAL*8,ALLOCATABLE::NODE2(:,:),PHI2(:)
!---CLAIM FOR BODY
      INTEGER IDOF(3)
      REAL*8 GAMMA,BETA,PC,BLX,BLY,THOB,MB,IB
	  REAL*8 RBC(5,2),COG(2),G_COOR(2),CLV(2),CRV(2),Z0_IC(3)
      REAL*8 FOR(3),FOR_PC(3),FOR_TEMP(3)
	  REAL*8,ALLOCATABLE::Z0(:),Z1(:),Z2(:),Z0_PC(:),Z1_PC(:),Z2_PC(:),Z2_TEMP(:)
!      REAL*8,ALLOCATABLE::AD(:,:),ED(:),BD(:)

      PI=DACOS(-1.D0)
      
!---SET INITIAL WALL TIME
!      CALL get_walltime(T_START)
      
      OPEN(UNIT=1,FILE='1.ipt',STATUS='OLD')
      OPEN(UNIT=2,FILE='2.ipt',STATUS='OLD')
      OPEN(UNIT=3,FILE='3.ipt',STATUS='OLD')
      OPEN(UNIT=4,FILE='4.ipt',STATUS='OLD')
      OPEN(UNIT=5,FILE='IO.DAT')
      OPEN(UNIT=61,FILE='S_wave.DAT')
      OPEN(UNIT=62,FILE='S_tld.DAT')
      OPEN(UNIT=71,FILE='WG_wave.DAT')
      OPEN(UNIT=72,FILE='WG_tld.DAT')
!      OPEN(UNIT=8,FILE='P.DAT')
      OPEN(UNIT=9,FILE='F.DAT')
      OPEN(UNIT=10,FILE='Z.DAT')
      OPEN(UNIT=11,FILE='RB.DAT')
!      OPEN(UNIT=12,FILE='DOMAIN.DAT')
      OPEN(UNIT=13,FILE='E.DAT')
	  OPEN(UNIT=20,FILE='code.log')      
	  OPEN(UNIT=21,FILE='ITER.DAT')
!	  OPEN(UNIT=23,FILE='CFL.DAT')
!	  OPEN(UNIT=24,FILE='PATM.DAT')
	  OPEN(UNIT=31,FILE='IC_config.DAT')
      OPEN(UNIT=32,FILE='IC_bound.DAT')
      OPEN(UNIT=99,FILE='TEST.TXT')

!---TOPOGRAGHY AND WAVE TYPE AND BODY CONFIGURATION
	CALL INPUT_2(NPL,NFL,NBT,WTYP,WMKTYP,NWG,WGX,NDOF,BLX,BLY,THOB,MB,IB,Z0_IC,IDOF)
	ALLOCATE(NELEM(NPL),ME(NPL),NS(NPL),BUNTYP(NPL),COOR(NPL,2),SIDE_L(NPL))
    ALLOCATE(Z0(3*NDOF),Z1(3*NDOF),Z2(3*NDOF),Z0_PC(3*NDOF),Z1_PC(3*NDOF),Z2_PC(3*NDOF),Z2_TEMP(3*NDOF)) !,AD(6*NDOF,6*NDOF),ED(6*NDOF),BD(6*NDOF)
	NELEM=0
	NS=0
	BUNTYP=0
	COOR=0.D0
	SIDE_L=0.D0
    Z0=0.D0
    Z1=0.D0
    Z2=0.D0
    Z0_PC=0.D0
    Z1_PC=0.D0
    Z2_PC=0.D0
    Z2_TEMP=0.D0
!    AD=0.D0
!    ED=0.D0
!    BD=0.D0
!    READ(4,*)  ((AD(I,J),J=1,6*NDOF),I=1,6*NDOF)
!    READ(4,*)  ((ED(I),I=1,6*NDOF)
    
!---INPUT PARAMETERS OF THE WAVE CHANNEL
	CALL INPUT_1(NPL,NFL,COOR,NFIELD,NNODE,NELM,NELEM,ME,NS,BUNTYP,NGA,&
               &GRAV,MU,WIDTH,THO,AMP,OMEGA,PSI,NTIM,DELTTIME,GAMMA,BETA,PC,ICON,&
               &NITER,ETOL,EATOL,CTYPE,DEP,BLX,BLY,G_COOR,CLV,CRV)
    ALLOCATE(LN(NELM,2),NODE(NNODE,2),NODE_PC(NNODE,2),NORM(NELM,2),JCB(NELM),LENG(NELM))
    ALLOCATE(PHI(NNODE),PHI_PC(NNODE),PPHI(NNODE),PHIT(NNODE),PPHIT(NNODE))
	ALLOCATE(KER1(NNODE,NNODE),KER2(NNODE,NNODE))
    ALLOCATE(DP(NNODE,2),D2PDT(NNODE),D2P(NNODE,2))
    ALLOCATE(G1(NNODE,NNODE),H1(NNODE,NNODE),EYE(NNODE))
	ALLOCATE(DPDS(NNODE),DPDT(NNODE),PR(NNODE),ACCMO(NNODE),PHIT_TEMP(NNODE))
	ALLOCATE(WT(NGA),RT(NGA),SHA1(NGA),SHA2(NGA),SH(2,NGA))
	ALLOCATE(DI(NTIM),VE(NTIM),AC(NTIM))
	LN=0
	NODE=0.D0
    NODE_PC=0.D0
	NORM=0.D0
	JCB=0.D0
	LENG=0.D0
	PHI=0.D0
    PHI_PC=0.D0
	PPHI=0.D0
	PHIT=0.D0
	PPHIT=0.D0
	KER1=0.D0
	KER2=0.D0
	DP=0.D0
	DPDS=0.D0
	DPDT=0.D0
    D2PDT=0.D0
    D2P=0.D0
	PR=0.D0
	ACCMO=0.D0
	PHIT_TEMP=0.D0
	WT=0.D0
	RT=0.D0
	SHA1=0.D0
	SHA2=0.D0
	SH=0.D0
	DI=0.D0
	VE=0.D0
	AC=0.D0
    G1=0.D0
    H1=0.D0
    EYE=1.D0

!---INPUT PARAMETERS OF THE SLOSHING TANK
    CALL INPUT_3(TLD_YN,COOR2,NNODE2,NELM2,NELEM2,ME2,NS2,BUNTYP2,MU1,MU2,POR,ICON2,TL,TH,TWN,TLDLOC,GRAV)
    ALLOCATE(LN2(NELM2,2),NODE2(NNODE2,2),PHI2(NNODE2))
    LN2=0
    NODE2=0.D0
    PHI2=0.D0
    WORK=0.D0

!---GAUSDSIN QUADRATURES AND SHAPE FUNCTION
    CALL GAUSS(WT,RT,NGA)
    CALL SHAP(SHA1,SHA2,SH,NGA,RT)
    
!---GIVE BODY'S INITIAL DISPLACEMENT AND VELOCITY
    Z0(1:2)=G_COOR+Z0_IC(1:2)
    Z0(3)=Z0_IC(3)
    
!---GENERATE WAVES: 1 = PERIODIC WAVE, 2 = SOLITARY WAVE (HALF-SINE IN A SHORT TIME)
	SELECT CASE (WTYP)
	CASE(1)
		DO I=1,NTIM
		 TIME=(I-1)*DELTTIME
		 DI(I)=AMP*DSIN(OMEGA*TIME+PSI)
		 VE(I)=AMP*OMEGA*DCOS(OMEGA*TIME+PSI)
		 AC(I)=-AMP*OMEGA**2*DSIN(OMEGA*TIME+PSI)
		 END DO
	CASE(2)
		NT=OMEGA/DELTTIME
		DO I=1,NT+1
		 TIME=(I-1)*DELTTIME
		 DI(I)=0.5D0*AMP-0.5D0*AMP*DCOS(PI/OMEGA*TIME)
		 VE(I)=PI/OMEGA*0.5*AMP*DSIN(PI/OMEGA*TIME)
		 AC(I)=(PI/OMEGA)**2*0.5D0*AMP*DCOS(PI/OMEGA*TIME)
		END DO
		 DI(NT+2:NTIM)=DI(NT+1)
    END SELECT
    
!---INITIALIZE WAVE FLOW PATTERN
	CALL LENGTH(NPL,COOR,SIDE_L)
	CALL MESH(NPL,NNODE,NELM,NELEM,LN,COOR,SIDE_L,NODE,Z0,CLV,CRV,BLX,BLY)
    CALL PRESSURE(ICON,THO,GRAV,DEP,NPL,NFL,NNODE,NS,NODE,PHIT,DP,PR,P_ATM)
    CALL KERNEL(KER1,KER2,NODE,NORM,JCB,LENG,LN,NNODE,NELM,NGA,SHA1,SHA2,SH,WT,EYE) ! CALL BCZ NEED JCB
    CALL FORCE_WAVE(TIME,GRAV,MB,NPL,NFL,NNODE,NELM,ME,LN,NODE,NORM,JCB,PR,Z0,FOR)
    FOR_TEMP=FOR    
!    DO I=1,2
!        Z2(I)=FOR(I)/MB
!    END DO
!    Z2(3)=FOR(3)/IB
    WRITE(20,*) 'PASS WAVE INITIALIZATION'
    
!---INITIALIZE PMTLD FLOW PATTERN
    IF (TLD_YN==1)THEN
	    CALL MESH_TLD(NNODE2,NELM2,NELEM2,ME2,NS2,LN2,TL,TH,TLDLOC,BLX,BLY,COOR2,NODE2,Z0,DEP)
        CALL PMTLD(NNODE2,NELM2,NGA,ICON2,NELEM2,ME2,NS2,LN2,BUNTYP2,MU1,MU2,POR,&
            &GRAV,MU,THO,DELTTIME,TL,TH,DEP,TLDLOC,COOR2,Z0,Z1,Z2,WT,SHA1,SHA2,SH,&
            &NODE2,PHI2,FOR2,KE,PE,TE,WORK,DTDT,DVDT,DEDT,DWDT)
        WRITE(20,*) 'PASS PMTLD INITIALIZATION'
    ELSE
        FOR2=0.D0
        WRITE(20,*) 'NO PMTLD INVOLVED'
    END IF
    
!---OUTPUT INITIAL STATE STRUCTURE
    WRITE(31,*) 'body dimensions = ',BLX,'x',BLY,'m'
    WRITE(31,*) 'body density = ',THOB,'kg/m^3'
    WRITE(31,*) 'body mass = ',MB,'kg'
    WRITE(31,*) 'body moment of inertia = ',IB,'kg-m^2'
    WRITE(31,*) 'PMTLD dimensions = ',TL,'x',TH,'m'
    WRITE(31,*) 'PMTLD mass = ',POR*THO*TL*TH,'kg'
    WRITE(31,*) 'PMTLD 1st mass = ',8.D0*POR*THO*DTANH(PI*TH/TL)/(PI/TL)**3/TL,'kg'
    WRITE(31,*) 'PMTLD 1st frequency = ',TWN,'rad/s'
    WRITE(31,*) 'PMTLD 1st damping ratio = ',MU1/(2.D0*TWN)*100.D0,'%'
    
    
!---OUTPUT INITIAL STATE OF WATER
    DO I=1,NNODE
    WRITE(32,"(I15,4(F15.8,2X))") I,NODE(I,:) !,PR(I)
    END DO
    DO I=1,NNODE2
    WRITE(32,"(I15,4(F15.8,2X))") I,NODE2(I,:) !,PR(I)
    END DO

!---WALL TIME FOR INITIALIZATION
!      CALL get_walltime(T_INITIAL)
!      WRITE(20,"(A1,F5.1,A16)") '[',T_INITIAL-T_START,'] INITIALIZATION'
    
!******************************** SIMULATION STARTS ********************************
DO NT=1,NTIM
    WRITE(*,*) NT,'TH'
    WRITE(20,*) NT,'TH'
    TIME=(NT-1)*DELTTIME    
!---GET STEP CONSUMED WALL TIME
!    CALL get_walltime(T_STEP)
!    WRITE(20,"(A1,F5.1,A1,I6,A2)") '[',T_STEP-T_START,']',NT,'TH'

!    WRITE(*,*) '***********'
!    WRITE(*,*) z0
!    WRITE(*,*) z1
!    WRITE(*,*) z2
!    WRITE(*,*) for
!    WRITE(*,*) for2
!    WRITE(*,*) '***********'
    
!---WAVEMAKER DIS, VEL, ACC, PHASE LAG (IN RADIUS)
    DIS=DI(NT)
    VEL=VE(NT)
    ACC=AC(NT)

!---CALCULATE WAVE SPEED OF THE OUTLET FOR RADIATION CONDITION
    D_OUT=NODE(NS(NFL+2),2)-NODE(NS(NPL-NBT-1),2)    
    CALL WAVE_SPD(GRAV,OMEGA,D_OUT,WC)
        
!************** WAVE AND BODY PROBLEM ******************        
!---GUESS INITIAL ACC/FOR AND NODE AND PHI
    FOR_PC=FOR
    Z0_PC=Z0
    Z1_PC=Z1
    Z2_PC=Z2
    NODE_PC=NODE
    PHI_PC=PHI
        
DO JTR=1,NITER

!---STORE THE ACC AND FOR OF THE KTH ITERATION TO COMPARE WITH (K+1)TH SOLUTION
    Z2_TEMP=Z2_PC
    FOR_TEMP=FOR_PC
      
!---BUILD KERNEL
    CALL KERNEL(KER1,KER2,NODE_PC,NORM,JCB,LENG,LN,NNODE,NELM,NGA,SHA1,SHA2,SH,WT,EYE)

!---IMPLICIT SOLVER FOR CONVERGENCE OF PHIT ON RADIATION OUTLET
!DO ITR=1,NITER
    PHIT_TEMP=PHIT

!---APPLY BC, SOLVING FREE-SURFACE PPHI, GET FIRST-ORDER TAYLOR SERIES EXPANSION
    CALL BOUND(WMKTYP,NDOF,NPL,NFL,NNODE,NELM,NELEM,ME,NS,LN,BUNTYP,NODE_PC,NORM,&
              &PHI_PC,PPHI,PHIT,VEL,Z0_PC,Z1_PC,WC)
!   CALL SOLVE_BACK(NPL,PHI_PC,PPHI,KER1,KER2,NNODE,NELEM,BUNTYP)        
    CALL SOLVE_LAPACK(NPL,PHI_PC,PPHI,KER1,KER2,NNODE,NELEM,BUNTYP)
!   CALL SOLVE_LAPACK2_1(NPL,PHI_PC,PPHI,KER1,KER2,G1,H1,NNODE,NELEM,BUNTYP)
    CALL TAYES1(NPL,NFL,NBT,NNODE,NELM,ME,NS,LN,NODE_PC,NORM,JCB,PHI_PC,PPHI,&
               &DPDS,DP,PHIT,DPDT,DEP,GRAV,MU,VEL,WC)

!---APPLY BC, SOLVING FREE-SURFACE PPHIT
    CALL ACCBC_KO(NPL,NNODE,NELM,NELEM,ME,LN,NORM,JCB,PHI_PC,PPHI,DP,ACCMO)    
!   CALL ACCBC_TSAO(NPL,NNODE,NELM,NELEM,ME,LN,NORM,JCB,PHI,PPHI,DP,ACCMO)
!   CALL ACCBC_GRILLI(NPL,NFL,NNODE,NELM,ME,LN,NODE_PC,NORM,JCB,PHI_PC,PPHI,Z0_PC,Z1_PC,ACCMO)        
!   CALL ACCBC_GUERBER(NPL,NFL,NNODE,NELM,ME,LN,NODE_PC,NORM,JCB,PHI_PC,PPHI,Z0_PC,Z1_PC,ACCMO)
!   CALL ACCBC_CAO(NPL,NFL,NNODE,NELM,ME,LN,NODE_PC,NORM,JCB,PHI_PC,PPHI,Z0_PC,Z1_PC,ACCMO)        
    CALL BOUNDT(WMKTYP,NDOF,NPL,NFL,NNODE,NELM,NELEM,ME,NS,LN,BUNTYP,NODE_PC,NORM,&
               &PHIT,PPHIT,DPDS,JCB,ACC,ACCMO,Z0_PC,Z1_PC,Z2_PC,WC)
!   CALL SOLVE_BACK(NPL,PHIT,PPHIT,KER1,KER2,NNODE,NELEM,BUNTYP)        
    CALL SOLVE_LAPACK(NPL,PHIT,PPHIT,KER1,KER2,NNODE,NELEM,BUNTYP)
!   CALL SOLVE_LAPACK2_2(NPL,PHIT,PPHIT,G1,H1,NNODE,NELEM,BUNTYP)

!    PHIT=PHIT+PC*(PHIT-PHIT_TEMP)
!   E1=MAXVAL(DABS(PHIT_TEMP(NS(NFL+2)+1:NS(NFL+3))-PHIT(NS(NFL+2)+1:NS(NFL+3))))
!   WRITE(*,*) JTR,ITR,E1
!END DO
!**************************************************************************************

!205 CONTINUE        

!---CALCULATE PRESSURE AND FORCE/MOMENTON ON ALL  BOUNDARIES (SOLID + FREE SURFACE)
    CALL PRESSURE(ICON,THO,GRAV,DEP,NPL,NFL,NNODE,NS,NODE_PC,PHIT,DP,PR,P_ATM)
    CALL FORCE_WAVE(TIME,GRAV,MB,NPL,NFL,NNODE,NELM,ME,LN,NODE_PC,NORM,JCB,PR,Z0_PC,FOR_PC)

!************** CORRECT FOR_PC ******************    
    FOR_PC=FOR_PC+PC*(FOR_TEMP-FOR_PC) ! v1 converged version looks normal
!    FOR_PC=FOR_PC+PC*(FOR_PC-FOR_TEMP) ! v2 peter seems converge faster
!    FOR_PC=FOR_TEMP+PC*(FOR-FOR_TEMP)  ! NG, CHECK Grilli
        
!************** FLOATING BODY DYNAMICS ******************
    CALL STR(IDOF,GAMMA,BETA,NDOF,TIME,DELTTIME,MB,IB,G_COOR,FOR2,FOR_PC,Z0,Z1,Z2,Z0_PC,Z1_PC,Z2_PC)

!---GET SECOND-ORDER TAYLOR SERIES EXPANSION
    CALL TAYES2(NPL,NFL,NNODE,NELM,NS,ME,LN,NODE_PC,NORM,JCB,PHI_PC,PPHI,PHIT,PPHIT,&
                &DP,DPDS,DPDT,D2PDT,D2P,DELTTIME,GRAV,ACC,Z2_PC)
        
!---REMESH LOCATION AND POTENTIAL OF FREE-SURFACE NODES
    CALL REMESH(NDOF,NPL,NFL,NBT,NELM,NNODE,NELEM,ME,NS,NODE,NODE_PC,DP,D2P,PHI,PHI_PC,DPDT,D2PDT,NORM,&
                &TIME,DELTTIME,AMP,NWG,WGX,CLV,CRV,Z0,RBC,COG)
        
!---CONVERGE OR CORRECT
!   E2=MAXVAL(DABS(PHIT_TEMP(NS(NPL-1)+1:NNODE)-PHIT(NS(NPL-1)+1:NNODE)))
!   IF(E1<=ETOL)THEN !.AND.E2<ETOL)THEN
!       WRITE(*,*) 'PASS RADIATION',ITR,'ITR'
!       WRITE(21,"(F9.5,I5,E15.8,E15.8)") TIME,ITR,E1,E2
!       GOTO 205
!   ELSE IF(ITR>=NITER)THEN
!       WRITE(20,"(A18,F9.5,I5,E15.8,E15.8)") 'RADIATION DIVERGE',TIME,E1LOC,E1,E2
!       WRITE(*,*) 'RADIATION DIVERGE'
!       STOP
!   END IF
!END DO

    SELECT CASE(CTYPE)
        CASE (1)
        EA=DABS(Z2_PC-Z2_TEMP)
        CASE (2)
        EA=DABS(FOR_PC-FOR_TEMP)
    END SELECT
    
    IF(MAXVAL(EA)<=EATOL)THEN
        WRITE(*,*) 'PASS FSI',JTR,'ITR'
        WRITE(20,*) 'PASS FSI',JTR,'ITR'
        WRITE(21,"(F9.5,I5,E15.8,E15.8)") TIME,ITR,E1,E2
        FOR=FOR_PC
        Z0=Z0_PC
        Z1=Z1_PC
        Z2=Z2_PC
        NODE=NODE_PC
        PHI=PHI_PC
        GOTO 215
    ELSE IF(JTR>=NITER)THEN
        WRITE(*,*) 'FSI ITERATION FAIL'
        WRITE(20,*) 'FSI ITERATION FAIL'
        WRITE(21,"(F9.5,I5,E15.8,E15.8)") TIME,ITR,E1,E2
        STOP
    END IF    
        
END DO
215 CONTINUE
    
!---CALCULATE FORCE AND MOMENT FROM PMTLD IF PMTLD IS INVOLVED
IF (TLD_YN==1)THEN
!***USED ONLY FOR TEST TLD***
!        z0(1)=5.D0+AMP*DSIN(OMEGA*TIME+PSI)
!        Z0(2)=2.D0
!        Z0(3)=0.D0
!        Z1=0.D0
!        Z2=0.D0
!        z1(1)=AMP*OMEGA*DCOS(OMEGA*TIME+PSI)
!        z2(1)=-AMP*OMEGA**2*DSIN(OMEGA*TIME+PSI)
!***USED ONLY FOR TEST TLD***
    CALL PMTLD(NNODE2,NELM2,NGA,ICON2,NELEM2,ME2,NS2,LN2,BUNTYP2,MU1,MU2,POR,&
              &GRAV,MU,THO,DELTTIME,TL,TH,DEP,TLDLOC,COOR2,Z0,Z1,Z2,WT,SHA1,SHA2,SH,&
              &NODE2,PHI2,FOR2,KE,PE,TE,WORK,DTDT,DVDT,DEDT,DWDT)
    WRITE(*,*) 'PASS PMTLD'
    WRITE(20,*) 'PASS PMTLD'
END IF
    
!---CHECK IF UNREASONABLY LARGE SOLUTION OCCUR (STOP IF YES)
    CALL STABLE(TIME,NNODE,NODE,NNODE2,NODE2)
    
!****************************** OUTPUT ******************************
WRITE(5,"(4(E15.8,1X))") TIME,DIS,VEL,ACC
WRITE(61,"(5000(E15.8,1X))") NODE(:,1)
WRITE(61,"(5000(E15.8,1X))") NODE(:,2)
WRITE(62,"(1000(E15.8,1X))") NODE2(:,1) 
WRITE(62,"(1000(E15.8,1X))") NODE2(:,2)
WRITE(9,"(10(E15.8,1X))") TIME,FOR(1),FOR(2)+MB*GRAV,FOR(3),FOR2,FOR+FOR2
WRITE(10,"(10(E15.8,1X))") TIME,Z0(1)-G_COOR(1),Z0(2)-G_COOR(2),Z0(3),Z1,Z2
WRITE(11,"(6(E15.8,1X))") RBC(:,1),COG(1)
WRITE(11,"(6(E15.8,1X))") RBC(:,2),COG(2)
WRITE(13,"(9(E15.8,1X))") TIME,KE,PE,TE,-WORK,DTDT,DVDT,DEDT,-DWDT

DO I=1,NWG
    CALL BWLOC(-WGX(I),NS(NFL+2),-NODE(1:NS(NFL+2),1),0,IR,IL)
    WGY(I)=NODE(IL,2)+(WGX(I)-NODE(IL,1))/(NODE(IR,1)-NODE(IL,1))*(NODE(IR,2)-NODE(IL,2))-DEP
END DO
WRITE(71,"(11(E15.8,1X))") TIME,WGY(1:NWG)*100.D0

WG1_TLD=DSQRT((NODE2(1,1)-NODE2(NS2(3),1))**2+(NODE2(1,2)-NODE2(NS2(3),2))**2)-TH
WG2_TLD=DSQRT((NODE2(NS2(1),1)-NODE2(NS2(2),1))**2+(NODE2(NS2(1),2)-NODE2(NS2(2),2))**2)-TH
WRITE(72,"(3(E15.8,1X))") TIME,WG1_TLD*100.D0,WG2_TLD*100.D0
!****************************** OUTPUT ******************************

END DO

!---GET PROGRAM CONSUMED WALL TIME
!    CALL get_walltime(T_END)
!    WRITE(20,"(A1,F5.1,A13)") '[',T_END-T_START,'] PROGRAM END'

STOP
END
!******************** END OF THE PROGRAM *************** 
!---CHECK THE CFL NUMBER	
!	  CALL COURANT(TIME,DELTTIME,NNODE,NELM,LN,NODE,DP,JCB)!++++

!---CALCULATE PRESSURE AND VELOCITY IN THE DOMAIN
!	  CALL DOMAIN(NPL,NGA,NFIELD,NNODE,NELM,NELEM,NS,LN,NODE,NORM,JCB,PHI,PPHI,PHIT,PPHIT,SHA1,SHA2,SH,WT,THO,GRAV,DEP,P_ATM,DP,PR) !++++++++++++++++++++++++++++++++++++++++++++++++

!**********************************************************************
      SUBROUTINE HEADLINE(ID,IREAD)
!**********************************************************************
      CHARACTER*2 ID
      ID =  '*'
      DO WHILE (ID .EQ. '*')
      READ(IREAD,'(A1)') ID
      END DO
      RETURN
      END
!**********************************************************************
      SUBROUTINE INPUT_1(NPL,NFL,COOR,NFIELD,NNODE,NELM,NELEM,ME,NS,BUNTYP,NGA,&
                      &GRAV,MU,WIDTH,THO,AMP,OMEGA,PSI,NTIM,DELTTIME,GAMMA,BETA,PC,ICON,&
                      &NITER,ETOL,EATOL,CTYPE,DEP,BLX,BLY,G_COOR,CLV,CRV)
!**********************************************************************
      IMPLICIT NONE
	  INTEGER I,J,NPL,NFL,NFIELD,NNODE,NELM,NGA,NTIM,ICON,NITER,CTYPE
	  INTEGER NELEM(NPL),ME(NPL),NS(NPL),BUNTYP(NPL)
	  REAL*8 ENDTIME,DELTTIME,GRAV,MU,WIDTH,THO,AMP,OMEGA,PSI,ETOL,EATOL
      REAL*8 DEP,BLX,BLY,GAMMA,BETA,PC,ET
      REAL*8 COOR(NPL,2),G_COOR(2),CLV(2),CRV(2)
      CHARACTER*2 ID
      ET=1.E-12
         ID = '*'
!---READ VERTEX NODE COORDINATES
         CALL HEADLINE(ID,1)
         READ(1,*)  ((COOR(I,J),J=1,2),I=1,NPL)
!---GIVE WATER DEPTH
         DEP=COOR(NPL,2)
!---READ ELEMENT MESH NUMBER
        CALL HEADLINE(ID,1)
        READ(1,*) (NELEM(I),I=1,NPL)
!---CALCULATE NUMBERS OF TOTAL NODE, TOTAL ELEMENT, ACCUMULATE NODE ON PLANE, ACCUMULATE ELEMENT ON PLANE
        NNODE=0
        NELM=0
        DO I=1,NPL
          NELM = NELM+NELEM(I)
          NNODE = NNODE+(NELEM(I)+1)
        END DO
		ME(1)=NELEM(1)
        NS(1)=NELEM(1)+1        
		DO I=2,NPL
		  ME(I)=ME(I-1)+NELEM(I)
          NS(I)=NS(I-1)+NELEM(I)+1          
		END DO
!---CALCULATE NUMBER OF NODE IN DOMAIN
	   NFIELD=(NELEM(1)-1)*(NELEM(2)-1)
!---DEFINE NODE INDEX FOR BERNOULLI CONSTANT
        IF(MOD(NS(1),2)==1)THEN
          ICON=(NS(1)+1)/2
        ELSE
          ICON=NS(1)/2
        END IF
!---BOUNDARY TYPE
         CALL HEADLINE(ID,1)
         READ(1,*) (BUNTYP(I),I=1,NPL)
!---READ THE GAUSSIAN INTEGRATION POINT NO.
         CALL HEADLINE(ID,1)
         READ(1,*)  NGA
!---READ gravity, mu, channel width, tho
         CALL HEADLINE(ID,1)
         READ(1,*)  GRAV,MU,WIDTH,THO
!---READ wavemaker parameters
         CALL HEADLINE(ID,1)
         READ(1,*)  AMP,OMEGA,PSI
!---READ endtime, delttime, niter, etol, etol, converge type
         CALL HEADLINE(ID,1)
         READ(1,*)  ENDTIME,DELTTIME,NITER,ETOL,EATOL,CTYPE
!---GIVE NUMBER OF TIME STEP
		 NTIM=ENDTIME/DELTTIME+1
!---READ gamma beta of Newmark method, coefficient of predictor
         CALL HEADLINE(ID,1)
         READ(1,*) GAMMA,BETA,PC
!---CHECK GEOMETRY CONFLICT        
     IF (DABS(BLX-DABS(COOR(NPL-1,1)-COOR(NPL-NFL-1,1)))>=ET)THEN
		WRITE(*,*) "FLOATING BODY LENGTH ERR"
        WRITE(20,*) "FLOATING BODY LENGTH ERR"
        STOP
     ELSE IF (BLY<DABS(COOR(NPL-1,2)-COOR(NPL-2,2)))THEN
        WRITE(*,*) "FLOATING BODY TOO SHORT"
        WRITE(20,*) "FLOATING BODY TOO SHORT"
        STOP
     END IF
!---CALCULATE COG COORDINATE
    G_COOR(1)=COOR(NPL-2,1)+0.5D0*BLX
    G_COOR(2)=COOR(NPL-2,2)+0.5D0*BLY
!---CALCULATE VECTOR OF CORNER-COG 
    CRV=COOR(NPL-3,:)-G_COOR
    CLV=COOR(NPL-2,:)-G_COOR

    RETURN
    END
!**********************************************************************
      SUBROUTINE INPUT_3(TLD_YN,COOR,NNODE,NELM,NELEM,ME,NS,BUNTYP,MU1,MU2,POR,ICON,TL,TH,WN,TLDLOC,GRAV)
!**********************************************************************
      IMPLICIT NONE
      INTEGER I,J,IREAD,NNODE,NELM,ICON,TLD_YN
      INTEGER NELEM(4),BUNTYP(4),ME(4),NS(4)
      REAL*8 COOR(4,2),MU1,MU2,POR,WN,GRAV,TL,TH,TLDLOC(2)
      CHARACTER*2 ID
      ID = '*'
      IREAD = 3
!---do you have PMTLD?
      CALL HEADLINE(ID,IREAD)
      READ(IREAD,*)  TLD_YN
!---input the coordinate of vertex of PMTLD domain
      CALL HEADLINE(ID,IREAD)
      READ(IREAD,*)  ((COOR(I,J),J=1,2),I=1,4)
      TL=COOR(2,1)-COOR(1,1)
      TH=COOR(3,2)-COOR(2,2)
!--- READ ELEMENT MESH NUMBER
      CALL HEADLINE(ID,IREAD)
      READ(IREAD,*) (NELEM(I),I=1,4)
!---CALCULATE NUMBERS OF TOTAL NODE, TOTAL ELEMENT, ACCUMULATE NODE ON PLANE, ACCUMULATE ELEMENT ON PLANE
      NNODE = 0
      NELM  = 0
      DO I=1,4
         NELM = NELM+NELEM(I)
         NNODE = NNODE+(NELEM(I)+1)
      END DO
      ME(1) = NELEM(1)
      NS(1) = NELEM(1)+1
      DO I=2,4
        ME(I)=ME(I-1)+NELEM(I)
        NS(I)=NS(I-1)+NELEM(I)+1          
      END DO
!---DEFINE NODE INDEX FOR BERNOULLI CONSTANT
        IF(MOD(NS(1),2)==1)THEN
          ICON=(NS(1)+1)/2
        ELSE
          ICON=NS(1)/2
        END IF      
!---BOUNDARY TYPE
      CALL HEADLINE(ID,IREAD)
      READ(IREAD,*) (BUNTYP(I),I=1,4)
!---READ linear damping ratio, mu2, porosity
      CALL HEADLINE(ID,IREAD)
      READ(IREAD,*)  MU1,MU2,POR
!---CALCULATE NATURAL FREQUENCY AND TRANSFER DAMPING EFFECT COEFFICIENT
      WN=SQRT(DACOS(-1.D0)*GRAV*TANH(DACOS(-1.D0)*COOR(3,2)/COOR(3,1))/COOR(3,1))
      MU1=MU1*2.D0*WN
!---*INPUT THE POSITION OF THE TLD RELATIVE TO COG OF THE BODY
      CALL HEADLINE(ID,IREAD)
      READ(IREAD,*) TLDLOC
      
      RETURN
    END
!**********************************************************************
      SUBROUTINE INPUT_4(NDOF,AD,ED,BD)
!**********************************************************************
      IMPLICIT NONE
      INTEGER I,J,IREAD,NDOF
      REAL*8 AD(6*NDOF,6*NDOF),ED(6*NDOF),BD(6*NDOF)
      CHARACTER*2 ID
      ID = '*'
      IREAD = 4
!---input state space matrices for equivalent mechanical PMTLD system
      CALL HEADLINE(ID,IREAD)
      READ(IREAD,*)  ((AD(I,J),J=1,6*NDOF),I=1,6*NDOF)
      READ(IREAD,*)  (ED(I),I=1,6*NDOF)
      READ(IREAD,*)  (BD(I),I=1,6*NDOF)

      RETURN
      END
!**********************************************************************
      SUBROUTINE INPUT_2(NPL,NFL,NBT,WTYP,WMKTYP,NWG,WGX,NDOF,BLX,BLY,THOB,MB,IB,Z0_IC,IDOF)
!**********************************************************************
      IMPLICIT NONE
	  INTEGER NPL,NFL,NBT,WTYP,WMKTYP,NWG,NDOF,IDOF(3)
	  REAL*8 BLX,BLY,THOB,MB,IB,WGX(10),Z0_IC(3)
      CHARACTER*2 ID
      ID = '*'
!---NUMBER OF PLANES
         CALL HEADLINE(ID,2)
         READ(2,*) NPL,NFL,NBT
!---WAVE GENERATION: 0=NO WAVE; 1=PERIODIC; 2=SOLITARY
         CALL HEADLINE(ID,2)
         READ(2,*) WTYP
!---WAVEMAKER TYPE: 1=PISTON; 2=FLAP
         CALL HEADLINE(ID,2)
         READ(2,*) WMKTYP
!---WAVE GENERATION: 1=PERIODIC; 2=SOLITARY
         CALL HEADLINE(ID,2)
         READ(2,*) NWG
!---CHECK WAVE GAUGE ALLOCATION
		 IF (NWG>10)THEN
		 WRITE(*,*) "NEED MORE ALLOCATION FOR WAVE GAUGE"
		 WRITE(20,*) "NEED MORE ALLOCATION FOR WAVE GAUGE"
		 STOP
         END IF
!---READ WAVE GAUGE X-LOCATION
         CALL HEADLINE(ID,2)
         READ(2,*) WGX(1:NWG)
!---NUMBER OF MASS OF THE FLOATING BODY
         CALL HEADLINE(ID,2)
         READ(2,*) NDOF
!---READ LENGTH, HEIGHT, AND DENSITY OF THE FLOATING BODY
         CALL HEADLINE(ID,2)
         READ(2,*) BLX,BLY,THOB
!---CALCULATE THE MASS AND MOMENT OF INERTIA OF THE RECTANGULAR FLOATING BODY
         MB=BLX*BLY*THOB
         IB=MB*(BLX**2+BLY**2)/12.D0
!---INPUT INITIAL DISPLACEMENT (NOTE THE ORIGIN LOACTES ON FREE SURFACE)
         CALL HEADLINE(ID,2)
         READ(2,*) Z0_IC(1:3)
         Z0_IC(3)=Z0_IC(3)*DACOS(-1.D0)/180.D0
!---INPUT DEGREE OF FREEDOM OF THE BODY ALONG 3 DIRECTIONS (1 = MOVABLE ; 0 = FIXED)
         CALL HEADLINE(ID,2)
         READ(2,*) IDOF
         
      RETURN
      END
!**********************************************************************
      SUBROUTINE LENGTH(NPL,COOR1,SIDE_L1)
!**********************************************************************
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
	  INTEGER NPL
      REAL*8  COOR1(NPL,2),SIDE_L1(NPL)
      DO I=1,NPL-1
        SIDE_L1(I)=DSQRT((COOR1(I+1,1)-COOR1(I,1))**2+(COOR1(I+1,2)-COOR1(I,2))**2)
      END DO
        SIDE_L1(NPL)=DSQRT((COOR1(NPL,1)-COOR1(1,1))**2+(COOR1(NPL,2)-COOR1(1,2))**2)

      RETURN
      END
!**********************************************************************
      SUBROUTINE MESH(NPL,NNODE,NELM,NELEM,LN,COOR,LENG,NODE,Z0,CLV,CRV,BLX,BLY)
!********************************************************************
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NPL,NEL(NPL),NNODE,NELM,NELEM(NPL),LN(NELM,2)
      REAL*8 SX,SY,NORM,DELT,BLX,BLY,LENG(NPL),COOR(NPL,2),NODE(NNODE,2)
      REAL*8 Z0(3),Q(2,2),CLV(2),CRV(2),RL(2),RR(2)

!------SET ROTATIONAL MATRIX
    Q(1,1)=DCOS(Z0(3))
    Q(1,2)=-DSIN(Z0(3))
    Q(2,1)=-Q(1,2)
    Q(2,2)=Q(1,1)
    
!------SET THE CORNER POINTS IN INERTIA FRAME (THIS BODY HAS 2 CORNERS)
    RL=MATMUL(Q,CLV)
    RR=MATMUL(Q,CRV)

    DO I=1,2
    COOR(NPL-2,I)=RL(I)+Z0(I)
    COOR(NPL-3,I)=RR(I)+Z0(I)
    END DO
    
    COOR(NPL-1,1)=Z0(1)-0.5*BLX/DCOS(Z0(3))
    COOR(NPL-4,1)=Z0(1)+0.5*BLX/DCOS(Z0(3))
    
!------MESH
    K=0
    DO I=1,NPL-1
	J=NPL-I
    NEL(I) = NELEM(I)+1
!    DELT=LENG(J)/NELEM(I)   
	SX=COOR(J,1)-COOR(J+1,1)
	SY=COOR(J,2)-COOR(J+1,2)
	NORM=DSQRT(SX**2+SY**2)
    DELT=NORM/NELEM(I)    
	SX=SX/NORM
	SY=SY/NORM
    DO L=1,NELEM(I)+1
       NODE(L+K,1)=COOR(J+1,1)+(L-1)*DELT*SX
       NODE(L+K,2)=COOR(J+1,2)+(L-1)*DELT*SY
    END DO
    K=K+NEL(I)
    END DO

    NEL(NPL) = NELEM(NPL)+1
    DELT=LENG(NPL)/NELEM(NPL)
	SX=COOR(NPL,1)-COOR(1,1)
	SY=COOR(NPL,2)-COOR(1,2)
	NORM=DSQRT(SX**2+SY**2)
	SX=SX/NORM
	SY=SY/NORM
    DO I=1,NELEM(NPL)+1
      NODE(I+K,1)=COOR(1,1)+(I-1)*DELT*SX
      NODE(I+K,2)=COOR(1,2)+(I-1)*DELT*SY
    END DO

!----TO GIVE THE LOCAL ELEMENT NODE NUMBER
      L=1
	  N=1
      DO I=1,NPL
       DO J=1,NELEM(I)
        LN(N,1)=L
        LN(N,2)=L+1
        L=L+1
		N=N+1
       END DO
       L=L+1
      END DO

      RETURN
    END
!**********************************************************************
      SUBROUTINE MESH_TLD(NNODE,NELM,NELEM,ME,NS,LN,TL,TH,TLDLOC,BLX,BLY,COOR,NODE,Z0,DEP)
!********************************************************************
      IMPLICIT NONE
      INTEGER I,J,L,N,NNODE,NELM,NELEM(4),ME(4),NS(4),LN(NELM,2)
      REAL*8 TL,TH,BLX,BLY,SX,SY,RL,LL,DELT,DEP,COOR(4,2),NODE(NNODE,2),Z0(3),TLDLOC(2),TEMP(2)

!---CHECK WAVE ELEVATION ON SIDE WALLS > 0      
      RL=TH-0.5D0*TL*DTAN(Z0(3))
      LL=TH+0.5D0*TL*DTAN(Z0(3))
      IF (RL<=0.D0.OR.LL<=0.D0)THEN
        WRITE(20,*) 'PMTLD SIDE DEPTH < 0'
        WRITE(*,*) 'PMTLD SIDE DEPTH < 0'
        STOP
      END IF

!------MESH FS
    DELT=TL/DCOS(Z0(3))/NELEM(1) ! FS REMAIN FLAT
    DO I=1,NS(1)
       NODE(I,1)=COOR(4,1)+(I-1)*DELT
       NODE(I,2)=COOR(4,2)
    END DO

!------MESH RIGHT WALL
    DELT=RL/NELEM(2)
    SX=DSIN(Z0(3))
    SY=-DCOS(Z0(3))
    J=1
    DO I=NS(1)+1,NS(2)
       NODE(I,1)=NODE(NS(1),1)+(J-1)*DELT*SX
       NODE(I,2)=NODE(NS(1),2)+(J-1)*DELT*SY
       J=J+1
    END DO

!------MESH BOTTOM
    DELT=TL/NELEM(3)
    SX=-DCOS(Z0(3))
    SY=-DSIN(Z0(3))
    J=1
    DO I=NS(2)+1,NS(3)
       NODE(I,1)=NODE(NS(2),1)+(J-1)*DELT*SX
       NODE(I,2)=NODE(NS(2),2)+(J-1)*DELT*SY
       J=J+1
    END DO   
    
!------MESH LEFT WALL
    DELT=LL/NELEM(4)
    SX=-DSIN(Z0(3))
    SY=DCOS(Z0(3))
    J=1
    DO I=NS(3)+1,NS(4)
       NODE(I,1)=NODE(NS(3),1)+(J-1)*DELT*SX
       NODE(I,2)=NODE(NS(3),2)+(J-1)*DELT*SY
       J=J+1
    END DO
    
!---SET THE TLD ON THE BOX BY LOCATING LEFT-BOTTOM CORNER
    TEMP=NODE(NS(3),:)
    DO I=1,NNODE
        NODE(I,1)=Z0(1)+(TLDLOC(1)-0.5D0*TL)*DCOS(Z0(3))-DSIN(Z0(3))*TLDLOC(2)+NODE(I,1)-TEMP(1)
        NODE(I,2)=Z0(2)+(TLDLOC(1)-0.5D0*TL)*DSIN(Z0(3))+DCOS(Z0(3))*TLDLOC(2)+NODE(I,2)-TEMP(2)
    END DO
    
!---SET Y0 ELEVATION FOR THE TLD
    DEP=NODE(1,2)
    
!----TO GIVE THE LOCAL ELEMENT NODE NUMBER
      L=1
	  N=1
      DO I=1,4
       DO J=1,NELEM(I)
        LN(N,1)=L
        LN(N,2)=L+1
        L=L+1
		N=N+1
       END DO
       L=L+1
      END DO

      RETURN
      END
!********************************************************************
      SUBROUTINE SHAP(SHA1,SHA2,SH,NGA,RT)
!********************************************************************
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NGA
      REAL*8 RT(NGA),SHA1(NGA),SHA2(NGA),SH(2,NGA)
      DO M=1,NGA
        SHA1(M)=0.5D0*(1-RT(M))
        SHA2(M)=0.5D0*(1+RT(M))

        SH(1,M)=SHA1(M)
        SH(2,M)=SHA2(M)
      END DO
	RETURN
	END 
!**********************************************************************
      SUBROUTINE REMESH(NDOF,NPL,NFL,NBT,NELM,NNODE,NELEM,ME,NS,NODE,NODE_PC,DP,D2P,PHI,PHI_PC,DPDT,D2PDT,NORM,&
                      &TIME,DELTTIME,AMP,NWG,WGX,CLV,CRV,Z0,RBC,COG)
!**********************************************************************
      IMPLICIT INTEGER(I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NDOF,NPL,NFL,NBT,NELM,NNODE,ME(NPL),NELEM(NPL),NS(NPL),NWG,IL,IR
      REAL*8 LENG,SX,SY,SNORM,HEL,HER,HE,DELTTIME,TEMP
      REAL*8 Q(2,2),CLV(2),CRV(2),CL(2),CR(2),RL(2),RR(2),RBC(5,2),COG(2)
      REAL*8 TIME,AMP,WGX(10),WGY(10),Z0(3*NDOF)
      REAL*8 NODE(NNODE,2),NODE_PC(NNODE,2),NORM(NELM,2),DP(NNODE,2),D2P(NNODE,2)
      REAL*8 PHI(NNODE),PHI_PC(NNODE),DPDT(NNODE),D2PDT(NNODE)

!-----REMESH FREE SURFACE NODE-----
    DO I=1,NS(1)
	NODE_PC(I,1)=NODE(I,1)+DP(I,1)*DELTTIME+0.5D0*D2P(I,1)*DELTTIME**2
    NODE_PC(I,2)=NODE(I,2)+DP(I,2)*DELTTIME+0.5D0*D2P(I,2)*DELTTIME**2
	PHI_PC(I)=PHI(I)+DELTTIME*DPDT(I)+0.5D0*D2PDT(I)*DELTTIME**2
    END DO

	DO I=NS(NFL+1)+1,NS(NFL+2)
	NODE_PC(I,1)=NODE(I,1)+DP(I,1)*DELTTIME+0.5D0*D2P(I,1)*DELTTIME**2
    NODE_PC(I,2)=NODE(I,2)+DP(I,2)*DELTTIME+0.5D0*D2P(I,2)*DELTTIME**2
	PHI_PC(I)=PHI(I)+DELTTIME*DPDT(I)+0.5D0*D2PDT(I)*DELTTIME**2
	END DO
      
!------ENSURE DUPLICATE POINT ON SIDE END NODE OF FREE SURFACE-----
    NODE_PC(NNODE,:)=NODE_PC(1,:)
	NODE_PC(NS(NFL+2)+1,:)=NODE_PC(NS(NFL+2),:)

!------BOTTOM END NODE GOES WITH FREE SURFACE -----
	NODE_PC(NS(NPL-1),1)=NODE_PC(1,1)
	NODE_PC(NS(NFL+3),1)=NODE_PC(NS(NFL+2),1)

!------SET ROTATIONAL MATRIX
    Q(1,1)=DCOS(Z0(3))
    Q(1,2)=-DSIN(Z0(3))
    Q(2,1)=-Q(1,2)
    Q(2,2)=Q(1,1)

!------SET THE CORNER POINTS IN INERTIA FRAME (THIS BODY HAS 2 CORNERS)
    RL=MATMUL(Q,CLV)
    RR=MATMUL(Q,CRV)
    DO I=1,2
    CL(I)=RL(I)+Z0(I)
    CR(I)=RR(I)+Z0(I)
    
    RBC(1,I)=CL(I)
    RBC(2,I)=CR(I)
    RBC(3,I)=-RL(I)+Z0(I)
    RBC(4,I)=-RR(I)+Z0(I)
    RBC(5,I)=RBC(1,I)    
    END DO

    COG=0.25D0*(RBC(1,:)+RBC(2,:)+RBC(3,:)+RBC(4,:))

!------REMESH FOR FLOATING BODY LEFT SIDE
    I=2 !THE Ith PLANE
	SNORM=DSQRT((CR(1)-CL(1))**2+(CR(2)-CL(2))**2)
	SX=-(CR(2)-CL(2))/SNORM
	SY=(CR(1)-CL(1))/SNORM
    HL=(NODE_PC(NS(I-1),1)-CL(1))*SX+(NODE_PC(NS(I-1),2)-CL(2))*SY
	DELT=HL/NELEM(I)
    K=1
    DO L=NS(I-1)+1,NS(I)
        NODE_PC(L,1)=CL(1)+HL*SX-DELT*(K-1)*SX
        NODE_PC(L,2)=CL(2)+HL*SY-DELT*(K-1)*SY
        K=K+1
    END DO

!---MAKE THE FS NODE STILL ON THE RIGID BODY PLANE
    NODE_PC(NS(I-1),:)=NODE_PC(NS(I-1)+1,:) ! KEEP PHI THE SAME FOR THE NEW LOCATION
    
!---CHECK THE FS NODE BESIDE THE INTERSECTION WENT INTO THE RIGID BODY OR NOT
    TEMP=(NODE_PC(NS(1)-1,1)-NODE_PC(NS(1),1))*NORM(ME(1)+1,1)+&
        &(NODE_PC(NS(1)-1,2)-NODE_PC(NS(1),2))*NORM(ME(1)+1,2)
    IF (TEMP>=0.D0)THEN
        WRITE(20,*) 'WAVE NODE PENETRATE - LEFT'
        WRITE(*,*) 'WAVE NODE PENETRATE - LEFT'
        STOP
    END IF
 
!------REMESH FOR FLOATING BODY BOTTOM
    I=3 !THE Ith PLANE
    SX=CR(1)-CL(1)
	SY=CR(2)-CL(2)
	SNORM=DSQRT(SX**2+SY**2)
	SX=SX/SNORM
	SY=SY/SNORM
	DELT=SNORM/NELEM(I)
    K=1
    DO L=NS(I-1)+1,NS(I)
        NODE_PC(L,1)=CL(1)+DELT*(K-1)*SX
        NODE_PC(L,2)=CL(2)+DELT*(K-1)*SY
        K=K+1
    END DO

!------REMESH FOR FLOATING BODY RIGHT SIDE
    I=4 !THE Ith PLANE
	SNORM=DSQRT((CR(1)-CL(1))**2+(CR(2)-CL(2))**2)
	SX=-(CR(2)-CL(2))/SNORM
	SY=(CR(1)-CL(1))/SNORM
    HR=(NODE_PC(NS(I)+1,1)-CR(1))*SX+(NODE_PC(NS(I)+1,2)-CR(2))*SY
	DELT=HR/NELEM(I)
    K=1
    DO L=NS(I-1)+1,NS(I)
        NODE_PC(L,1)=CR(1)+DELT*(K-1)*SX
        NODE_PC(L,2)=CR(2)+DELT*(K-1)*SY
        K=K+1
    END DO

!---MAKE THE FS NODE STILL ON THE RIGID BODY PLANE
    NODE_PC(NS(I)+1,:)=NODE_PC(NS(I),:) ! KEEP PHI THE SAME FOR THE NEW LOCATION
    
!---CHECK THE FS NODE BESIDE THE INTERSECTION WENT INTO THE RIGID BODY OR NOT
    TEMP=(NODE_PC(NS(1+NFL)+2,1)-NODE_PC(NS(1+NFL)+1,1))*NORM(ME(1+NFL),1)+&
        &(NODE_PC(NS(1+NFL)+2,2)-NODE_PC(NS(1+NFL)+1,2))*NORM(ME(1+NFL),2)
    IF (TEMP>=0.D0)THEN
        WRITE(20,*) 'WAVE NODE PENETRATE - RIGHT'
        WRITE(*,*) 'WAVE NODE PENETRATE - RIGHT'
        STOP
    END IF
    
!------REMESH FOR ALL PLANES EXCEPT FLOATING BODY AND FREE SURFACE
DO I=NFL+3,NPL
	SX=NODE_PC(NS(I),1)-NODE_PC(NS(I-1),1)
	SY=NODE_PC(NS(I),2)-NODE_PC(NS(I-1),2)
	SNORM=DSQRT(SX**2+SY**2)
	SX=SX/SNORM
	SY=SY/SNORM
	DELT=SNORM/NELEM(I)
	  K=1
      DO L=NS(I-1)+1,NS(I)
        NODE_PC(L,1)=NODE_PC(NS(I-1),1)+DELT*(K-1)*SX
        NODE_PC(L,2)=NODE_PC(NS(I-1),2)+DELT*(K-1)*SY
		K=K+1
      END DO
END DO

    RETURN
    END
!********************************************************************
SUBROUTINE BOUND(WMKTYP,NDOF,NPL,NFL,NNODE,NELM,NELEM,ME,NS,LN,BUNTYP,NODE,NORM,PHI,PPHI,PHIT,VEL,Z0,Z1,WC)
!********************************************************************
IMPLICIT NONE
INTEGER I,J,K,N,WMKTYP,NDOF,NPL,NFL,NELM,NNODE,NELEM(NPL),ME(NPL),NS(NPL),LN(NELM,2),BUNTYP(NPL)
REAL*8 VX,VY,DX,DY,R,VEL,WC
REAL*8 NODE(NNODE,2),NORM(NELM,2),PHI(NNODE),PPHI(NNODE),PHIT(NNODE)
REAL*8 Z0(3*NDOF),Z1(3*NDOF)

N=0
DO I=1,NPL ! PLANE NUMBER
    DO J=N+1,ME(I) ! ELEMENT NUMBER
        DO K=1,2   ! LN(J,K) = NODE NUMBER
            IF (BUNTYP(I) .EQ. 1) THEN          ! BC FOR FREE SURFACE
                PHI(LN(J,K))=PHI(LN(J,K))
                PPHI(LN(J,K))=0.D0
            ELSE
                IF (I==(NFL+3))THEN               ! BC FOR RADIATION SIDE
                    PPHI(LN(J,K))=0.D0 !-PHIT(LN(J,K))/WC !
                    PHI(LN(J,K))=0.D0
                ELSE IF (I==NPL)THEN              ! BC FOR WAVE MAKER
                    IF (WMKTYP==1)THEN
                        PPHI(LN(J,K))=-VEL !0.D0 !-PHIT(LN(J,K))/WC !
                        PHI(LN(J,K))=0.D0
                    ELSE IF (WMKTYP==2)THEN
                        PPHI(LN(J,K))=-VEL*DSQRT((NODE(LN(J,K),1)-NODE(NS(NPL-1),1))**2+&
                        &(NODE(LN(J,K),2)-NODE(NS(NPL-1),2))**2) !0.D0 !-PHIT(LN(J,K))/WC !
                        PHI(LN(J,K))=0.D0
                    END IF
                ELSE IF (I>1.AND.I<=(NFL+1))THEN  ! BC FOR RIGID BODY
                    DX=NODE(LN(J,K),1)-Z0(1)
                    DY=NODE(LN(J,K),2)-Z0(2)
                    VX=Z1(1)-Z1(3)*DY
                    VY=Z1(2)+Z1(3)*DX
                    PPHI(LN(J,K))=VX*NORM(J,1)+VY*NORM(J,2)
                    PHI(LN(J,K))=0.D0
                ELSE                              ! BC FOR BOTTOM
                    PPHI(LN(J,K))=0.D0
                    PHI(LN(J,K))=0.D0
                END IF
            END IF
        END DO
    END DO
    N=ME(I)
END DO

RETURN
END
!********************************************************************
      SUBROUTINE KERNEL(KER1,KER2,NODE,NORM,JCB,LENG,LN,NNODE,NELM,NGA,SHA1,SHA2,SH,WT,EYE)
!********************************************************************
      IMPLICIT NONE
      INTEGER I,J,K,L,M,N,NNODE,NELM,NGA,LN(NELM,2)
      REAL*8 RD,SIGMA(NNODE),EYE(NNODE)
      REAL*8  KER1(NNODE,NNODE),KER2(NNODE,NNODE)
      REAL*8  NX,NY,H(2),G(2),XFUNC(10),YFUNC(10),PXI1(2)
      REAL*8  LENG(NELM),NORM(NELM,2),JCB(NELM),NODE(NNODE,2)
      REAL*8  WT(NGA),SHA1(NGA),SHA2(NGA),SH(2,NGA)

        KER1=0.D0
        KER2=0.D0
!**** CALCULATE THE JACOBIAN
      DO J=1,NELM
      LENG(J)=DSQRT((NODE(LN(J,1),1)-NODE(LN(J,2),1))**2+(NODE(LN(J,1),2)-NODE(LN(J,2),2))**2)
		DO L=1,2
          PXI1(L)=(-0.5D0)*NODE(LN(J,1),L)+ 0.5D0*NODE(LN(J,2),L)
        END DO
        NX=PXI1(2)
        NY=PXI1(1)
        JCB(J)=DSQRT(NX**2+NY**2)
        NORM(J,1)=-NX/JCB(J)
        NORM(J,2)=NY/JCB(J)
      END DO
      
!$omp parallel do private(I,J,K,M,XFUNC,YFUNC,RD,G,H)
!***THE SURFACE KERNELS
      DO I = 1,NNODE
       DO J=1,NELM
       DO M=1,NGA
          XFUNC(M)=SHA1(M)*NODE(LN(J,1),1)+SHA2(M)*NODE(LN(J,2),1)
          YFUNC(M)=SHA1(M)*NODE(LN(J,1),2)+SHA2(M)*NODE(LN(J,2),2)
       END DO
        DO K=1,2
         G(K)=0.D0
         H(K)=0.D0
         RD=DSQRT((NODE(I,1)-NODE(LN(J,K),1))**2+(NODE(I,2)-NODE(LN(J,K),2))**2)
        IF (RD .LE. 0.000001D0) THEN
        H(K)=0.D0
	    G(K)=LENG(J)/2*(1.5D0-DLOG(LENG(J)))
		ELSE !---NON DSINGULER TERM
         G(K)=0.D0
         DO M=1,NGA
            H(K)=H(K)+(-1.D0)/((XFUNC(M)-NODE(I,1))**2+(YFUNC(M)-NODE(I,2))**2)*&
                    &((XFUNC(M)-NODE(I,1))*NORM(J,1)+(YFUNC(M)-NODE(I,2))*NORM(J,2))&
                    &*JCB(J)*SH(K,M)*WT(M)
            G(K)=G(K)+DLOG(1.D0/((XFUNC(M)-NODE(I,1))**2+(YFUNC(M)-NODE(I,2))**2)**0.5D0)*JCB(J)*SH(K,M)*WT(M)
         END DO
      END IF
         KER1(I,LN(J,K))=KER1(I,LN(J,K))+H(K)
         KER2(I,LN(J,K))=KER2(I,LN(J,K))+G(K)
      END DO
      END DO
      END DO
!$omp end parallel do
      
!***DSINGULAR OF KER1
    DO I=1,NNODE
    KER1(I,I)=0.D0
    END DO      
      CALL DGEMM('N','N',NNODE,1,NNODE,1.D0,KER1,NNODE,EYE,NNODE,0.D0,SIGMA,NNODE)
    DO I=1,NNODE
    KER1(I,I)=-SIGMA(I)
    END DO  
!       DO I=1,NNODE
!          SIGMA=0.D0
!          DO J=1,NNODE
!             SIGMA=KER1(I,J)+SIGMA
!          END DO
!          KER1(I,I)=-SIGMA
!       END DO

      RETURN
      END
!**********************************************************************
       SUBROUTINE SOLVE_LAPACK(NPL,PHI,PPHI,KER1,KER2,NNODE,NELEM,BUNTYP)
!**********************************************************************
       IMPLICIT NONE
       INTEGER I,J,K,L,N,NPL,NNODE,NELEM(NPL),BUNTYP(NPL)
	   INTEGER INFO,IPIV(NNODE)
       REAL*8 PHI(NNODE),PPHI(NNODE),KER1(NNODE,NNODE),KER2(NNODE,NNODE)
       REAL*8 H1(NNODE,NNODE),Q1(NNODE),G1(NNODE,NNODE),P1(NNODE)

!-----MOVE PHI AND PPHI----         
       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,(NELEM(I)+1)+N
            IF (BUNTYP(I) .EQ. 1) THEN
               Q1(J)=PHI(J)
            ELSE
               Q1(J)=PPHI(J)
            END IF
          END DO
          N=N+(NELEM(I)+1)
       END DO

!-----MOVE KER1 AND KER2---- 
         DO I=1,NNODE
           N=0
           DO L=1,NPL
             DO J=K+N,(NELEM(L)+1)+N
             IF (BUNTYP(L) .EQ. 1) THEN
               H1(I,J)=KER1(I,J)  !G*UNKNOWN=H*KNOWN
               G1(I,J)=KER2(I,J)
              ELSE
               H1(I,J)=-KER2(I,J)
               G1(I,J)=-KER1(I,J)
             END IF
             END DO
            N=N+(NELEM(L)+1)
          END DO
        END DO

       P1=MATMUL(H1,Q1)

!*************SOLVE BY CALLING LAPACK*********
CALL DGESV(NNODE,1,G1,NNODE,IPIV,P1,NNODE,INFO)

       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,(NELEM(I)+1)+N
            IF (BUNTYP(I) .EQ. 1) THEN
               PPHI(J)=P1(J)
            ELSE
               PHI(J)=P1(J)
            END IF
          END DO
          N=N+NELEM(I)+1
       END DO

      RETURN
    END
!**********************************************************************
       SUBROUTINE SOLVE_LAPACK2_1(NPL,PHI,PPHI,KER1,KER2,G1,H1,NNODE,NELEM,BUNTYP)
!**********************************************************************
       IMPLICIT NONE
       INTEGER I,J,K,L,N,NPL,NNODE,NELEM(NPL),BUNTYP(NPL)
	   INTEGER INFO,IPIV(NNODE)
       REAL*8 PHI(NNODE),PPHI(NNODE),KER1(NNODE,NNODE),KER2(NNODE,NNODE)
       REAL*8 H1(NNODE,NNODE),Q1(NNODE),G1(NNODE,NNODE),P1(NNODE)
	   CHARACTER*1 TRANS
       TRANS = 'N'
       
!-----MOVE PHI AND PPHI----         
       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,(NELEM(I)+1)+N
            IF (BUNTYP(I) .EQ. 1) THEN
               Q1(J)=PHI(J)
            ELSE
               Q1(J)=PPHI(J)
            END IF
          END DO
          N=N+(NELEM(I)+1)
       END DO
       
!-----MOVE KER1 AND KER2---- 
         DO I=1,NNODE
           N=0
           DO L=1,NPL
             DO J=K+N,(NELEM(L)+1)+N
             IF (BUNTYP(L) .EQ. 1) THEN
               H1(I,J)=KER1(I,J)  !G*UNKNOWN=H*KNOWN
               G1(I,J)=KER2(I,J)
              ELSE
               H1(I,J)=-KER2(I,J)
               G1(I,J)=-KER1(I,J)
             END IF
             END DO
            N=N+(NELEM(L)+1)
          END DO
        END DO

       P1=MATMUL(H1,Q1)
       
!*************SOLVE BY CALLING LAPACK*********
CALL DGETRF(NNODE,NNODE,G1,NNODE,IPIV,INFO)
CALL DGETRS(TRANS,NNODE,1,G1,NNODE,IPIV,P1,NNODE,INFO)
       
       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,(NELEM(I)+1)+N
            IF (BUNTYP(I) .EQ. 1) THEN
               PPHI(J)=P1(J)
            ELSE
               PHI(J)=P1(J)
            END IF
          END DO
          N=N+NELEM(I)+1
       END DO
                     
      RETURN
    END
!**********************************************************************
       SUBROUTINE SOLVE_LAPACK2_2(NPL,PHI,PPHI,G1,H1,NNODE,NELEM,BUNTYP)
!**********************************************************************
       IMPLICIT NONE
       INTEGER I,J,K,L,N,NPL,NNODE,NELEM(NPL),BUNTYP(NPL)
	   INTEGER INFO,IPIV(NNODE)
       REAL*8 PHI(NNODE),PPHI(NNODE)
       REAL*8 H1(NNODE,NNODE),Q1(NNODE),G1(NNODE,NNODE),P1(NNODE)
	   CHARACTER*1 TRANS
       TRANS = 'N'
       
!-----MOVE PHI AND PPHI----         
       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,(NELEM(I)+1)+N
            IF (BUNTYP(I) .EQ. 1) THEN
               Q1(J)=PHI(J)
            ELSE
               Q1(J)=PPHI(J)
            END IF
          END DO
          N=N+(NELEM(I)+1)
       END DO

       P1=MATMUL(H1,Q1)

!*************SOLVE BY CALLING LAPACK*********
CALL DGETRS(TRANS,NNODE,1,G1,NNODE,IPIV,P1,NNODE,INFO)

       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,(NELEM(I)+1)+N
            IF (BUNTYP(I) .EQ. 1) THEN
               PPHI(J)=P1(J)
            ELSE
               PHI(J)=P1(J)
            END IF
          END DO
          N=N+NELEM(I)+1
       END DO

      RETURN
      END
!********************************************************************
      SUBROUTINE TAYES1(NPL,NFL,NBT,NNODE,NELM,ME,NS,LN,NODE,NORM,JCB,&
                       &PHI,PPHI,DPDS,DP,PHIT,DPDT,DEP,GRAV,MU,VEL,WC)
!********************************************************************
    IMPLICIT NONE
    INTEGER I,J,K,NPL,NFL,NBT,NNODE,NELM,ME(NPL),NS(NPL),LN(NELM,2)
    REAL*8 WC,DEP,GRAV,MU,VEL,TEMP1,TEMP2
	REAL*8 NODE(NNODE,2),NORM(NELM,2),JCB(NELM),PHI(NNODE),PPHI(NNODE)
    REAL*8 DP(NNODE,2),DPDS(NNODE),PHIT(NNODE),DPDT(NNODE)

	DPDS=0.D0
!*********************ON RIGID BODY*********************
    DO I=ME(1)+1,ME(NFL+1)
    DO J=1,2
        DPDS(LN(I,J))=(-0.5D0*PHI(LN(I,1))+0.5D0*PHI(LN(I,2)))/JCB(I)
        DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
	    DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
    END DO
    END DO

!********************* ON FREE SURFACE 1 *********************
    DO I=1,ME(1)
    DO J=1,2
		IF(LN(I,J).EQ.1) THEN
            DP(LN(I,J),1)=-PPHI(NNODE)
            DPDS(LN(I,J))=(DP(LN(I,J),1)-PPHI(LN(I,J))*NORM(I,1))/NORM(I,2)
            DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
        ELSE IF(LN(I,J).EQ.NS(1)) THEN
            TEMP1=PPHI(NS(1)+1)-PPHI(NS(1))*(NORM(I,1)*NORM(I+1,1)+NORM(I,2)*NORM(I+1,2))
            TEMP2=NORM(I,2)*NORM(I+1,1)-NORM(I,1)*NORM(I+1,2)
            DPDS(LN(I,J))=TEMP1/TEMP2
			DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE
			DPDS(LN(I,J))=(-0.5D0*PHI(LN(I,1))+0.5D0*PHI(LN(I,2)))/JCB(I)
			DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		END IF
    END DO
    END DO

    DO I=1,NS(1)
	  DPDT(I)=0.5D0*(DP(I,1)**2+DP(I,2)**2)-GRAV*(NODE(I,2)-DEP)-MU*PHI(I) !-GRAV*NODE(I,2)-MU*PHI(I) !
      PHIT(I)=DPDT(I)-(DP(I,1)**2+DP(I,2)**2)
    END DO

!********************* ON FREE SURFACE 2 *********************
    DO I=ME(NFL+1)+1,ME(NFL+2)
    DO J=1,2
		IF(LN(I,J).EQ.NS(NFL+1)+1) THEN
            TEMP1=PPHI(NS(NFL+1))-PPHI(NS(NFL+1)+1)*(NORM(I,1)*NORM(I-1,1)+NORM(I,2)*NORM(I-1,2))
            TEMP2=NORM(I,2)*NORM(I-1,1)-NORM(I,1)*NORM(I-1,2)
            DPDS(LN(I,J))=TEMP1/TEMP2
			DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE IF(LN(I,J).EQ.NS(NFL+2)) THEN
		    DP(LN(I,J),1)=PPHI(NS(NFL+2)+1)
            DPDS(LN(I,J))=(DP(LN(I,J),1)-PPHI(LN(I,J))*NORM(I,1))/NORM(I,2)
            DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE
			DPDS(LN(I,J))=(-0.5D0*PHI(LN(I,1))+0.5D0*PHI(LN(I,2)))/JCB(I)
			DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		END IF
    END DO
    END DO

    DO I=NS(NFL+1)+1,NS(NFL+2)
	  DPDT(I)=0.5D0*(DP(I,1)**2+DP(I,2)**2)-GRAV*(NODE(I,2)-DEP)-MU*PHI(I) !-GRAV*NODE(I,2)-MU*PHI(I) !
      PHIT(I)=DPDT(I)-(DP(I,1)**2+DP(I,2)**2)
    END DO

!*********************ON RIGHT WALL*********************
    DO I=ME(NFL+2)+1,ME(NFL+3)
    DO J=1,2
		IF(LN(I,J).EQ.NS(NFL+2)+1) THEN
			DPDS(LN(I,J))=-DP(NS(NFL+2),2)
            DP(LN(I,J),1)=PPHI(LN(I,J))
            DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE IF(LN(I,J).EQ.NS(NFL+3)) THEN
			DPDS(LN(I,J))=0.D0
            DP(LN(I,J),1)=PPHI(LN(I,J))
            DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE
			DPDS(LN(I,J))=(-0.5D0*PHI(LN(I,1))+0.5D0*PHI(LN(I,2)))/JCB(I)
			DP(LN(I,J),1)=PPHI(LN(I,J))
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		END IF
    END DO
    END DO

!*********************LEFT WALL*********************
    DO I=ME(NPL-1)+1,ME(NPL)
    DO J=1,2
		IF(LN(I,J).EQ.NS(NPL-1)+1) THEN
			DPDS(LN(I,J))=0.D0
			DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE IF(LN(I,J).EQ.NNODE) THEN
			DPDS(LN(I,J))=DP(1,2)
			DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE
			DPDS(LN(I,J))=(-0.5D0*PHI(LN(I,1))+0.5D0*PHI(LN(I,2)))/JCB(I)
			DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		END IF
    END DO
    END DO

!*********************BOTTOM*********************
    DO I=ME(NPL-NBT-1)+1,ME(NPL-1)
    DO J=1,2
		IF(LN(I,J).EQ.NS(NPL-NBT-1)+1) THEN
          DPDS(LN(I,J))=-PPHI(NS(NPL-NBT-1))
		  DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
          DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE IF(LN(I,J).EQ.NS(NPL-1)) THEN
          DPDS(LN(I,J))=PPHI(NS(NPL-1)+1)
		  DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
          DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE
		  DPDS(LN(I,J))=(-0.5D0*PHI(LN(I,1))+0.5D0*PHI(LN(I,2)))/JCB(I)
		  DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
		  DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		END IF
    END DO
    END DO

      RETURN
    END
!********************************************************************
      SUBROUTINE ACCBC_TSAO(NPL,NNODE,NELM,NELEM,ME,LN,NORM,JCB,PHI,PPHI,DP,ACCMO)
!********************************************************************
    IMPLICIT NONE
    INTEGER I,J,K,N,NPL,NNODE,NELM,NELEM(NPL),ME(NPL),LN(NELM,2)
    REAL*8 DPNDX,DPNDY,DPNDS,DPDS
	REAL*8 NORM(NELM,2),JCB(NELM),PHI(NNODE),PPHI(NNODE),DP(NNODE,2),ACCMO(NNODE)
N=0
DO K=1,NPL
  DO I=N+1,ME(K)
    DO J=1,2
		DPNDS=(-0.5D0*PPHI(LN(I,1))+0.5D0*PPHI(LN(I,2)))/JCB(I)
		DPDS=(-0.5D0*PHI(LN(I,1))+0.5D0*PHI(LN(I,2)))/JCB(I)
        ACCMO(LN(I,J))=DPNDS*DPDS/JCB(I)-DPNDS*DPDS ! PHINN=PHISS=0 FOR LINEAR ELEMENT
!        ACCMO(LN(I,J))=DPNDS*DPDS/JCB(I)-DPNDS*(DP(LN(I,J),1)*NORM(I,2)-DP(LN(I,J),2)*NORM(I,1))
    END DO
  END DO
  N=ME(K)
END DO

      RETURN
    END
!********************************************************************
      SUBROUTINE ACCBC_KO(NPL,NNODE,NELM,NELEM,ME,LN,NORM,JCB,PHI,PPHI,DP,ACCMO)
!********************************************************************
    IMPLICIT NONE
    INTEGER I,J,K,N,NPL,NNODE,NELM,NELEM(NPL),ME(NPL),LN(NELM,2)
    REAL*8 DPNDX,DPNDY,DPNDS,DPDS
	REAL*8 NORM(NELM,2),JCB(NELM),PHI(NNODE),PPHI(NNODE),DP(NNODE,2),ACCMO(NNODE)
N=0
DO K=1,NPL
  DO I=N+1,ME(K)
    DO J=1,2
		DPNDS=(-0.5D0*PPHI(LN(I,1))+0.5D0*PPHI(LN(I,2)))/JCB(I)
		DPDS=(-0.5D0*PHI(LN(I,1))+0.5D0*PHI(LN(I,2)))/JCB(I)
        ACCMO(LN(I,J))=-DPNDS*DPDS ! PHINN=PHISS=0 FOR LINEAR ELEMENT
!		DPNDX=DPNDS*NORM(I,2)	
!		DPNDY=-DPNDS*NORM(I,1)
!		ACCMO(LN(I,J))=DP(LN(I,J),1)*DPNDX+DP(LN(I,J),2)*DPNDY
    END DO
  END DO
  N=ME(K)
END DO

      RETURN
    END
!********************************************************************
      SUBROUTINE ACCBC_GUERBER(NPL,NFL,NNODE,NELM,ME,LN,NODE,NORM,JCB,PHI,PPHI,Z0,Z1,ACCMO)
!********************************************************************
    IMPLICIT NONE
    INTEGER I,J,NPL,NFL,NNODE,NELM,ME(NPL),LN(NELM,2)
    REAL*8 DX,DY,VX,VY,DPNDS,DPDS,Z0(3),Z1(3)
	REAL*8 NODE(NNODE,2),NORM(NELM,2),JCB(NELM),PHI(NNODE),PPHI(NNODE),ACCMO(NNODE)

  DO I=ME(2)+1,ME(NFL+1)
    DO J=1,2
		DPNDS=(-0.5D0*PPHI(LN(I,1))+0.5D0*PPHI(LN(I,2)))/JCB(I)
		DPDS=(-0.5D0*PHI(LN(I,1))+0.5D0*PHI(LN(I,2)))/JCB(I)
        DX=NODE(LN(I,J),1)-Z0(1)
        DY=NODE(LN(I,J),2)-Z0(2)
        VX=Z1(1)-Z1(3)*DY
        VY=Z1(2)+Z1(3)*DX
        ACCMO(LN(I,J))=Z1(3)*(VX*NORM(I,2)-VY*NORM(I,1)-DPDS)-DPNDS*(VX*NORM(I,2)-VY*NORM(I,1))
    END DO
  END DO

      RETURN
    END
!********************************************************************
      SUBROUTINE ACCBC_GRILLI(NPL,NFL,NNODE,NELM,ME,LN,NODE,NORM,JCB,PHI,PPHI,Z0,Z1,ACCMO)
!********************************************************************
    IMPLICIT NONE
    INTEGER I,J,NPL,NFL,NNODE,NELM,ME(NPL),LN(NELM,2)
    REAL*8 DX,DY,VX,VY,DPNDS,DPDS,Z0(3),Z1(3)
	REAL*8 NODE(NNODE,2),NORM(NELM,2),JCB(NELM),PHI(NNODE),PPHI(NNODE),ACCMO(NNODE)

  DO I=ME(2)+1,ME(NFL+1)
    DO J=1,2
		DPNDS=(-0.5D0*PPHI(LN(I,1))+0.5D0*PPHI(LN(I,2)))/JCB(I)
		DPDS=(-0.5D0*PHI(LN(I,1))+0.5D0*PHI(LN(I,2)))/JCB(I)
        
        DX=NODE(LN(I,J),1)-Z0(1)
        DY=NODE(LN(I,J),2)-Z0(2)
        VX=Z1(1)-Z1(3)*DY
        VY=Z1(2)+Z1(3)*DX
        ACCMO(LN(I,J))=-DPNDS*(VX*NORM(I,2)-VY*NORM(I,1))*(JCB(I)+1.D0)/JCB(I)+DPNDS*DPDS/JCB(I)
    END DO
  END DO

      RETURN
    END 
!********************************************************************
      SUBROUTINE ACCBC_CAO(NPL,NFL,NNODE,NELM,ME,LN,NODE,NORM,JCB,PHI,PPHI,Z0,Z1,ACCMO)
!********************************************************************
    IMPLICIT NONE
    INTEGER I,J,NPL,NFL,NNODE,NELM,ME(NPL),LN(NELM,2)
    REAL*8 DX,DY,VX,VY,DPNDS,DPDS,Z0(3),Z1(3)
	REAL*8 NODE(NNODE,2),NORM(NELM,2),JCB(NELM),PHI(NNODE),PPHI(NNODE),ACCMO(NNODE)

  DO I=ME(2)+1,ME(NFL+1)
    DO J=1,2
		DPNDS=(-0.5D0*PPHI(LN(I,1))+0.5D0*PPHI(LN(I,2)))/JCB(I)
        
		DPDS=(-0.5D0*PHI(LN(I,1))+0.5D0*PHI(LN(I,2)))/JCB(I)
        DX=NODE(LN(I,J),1)-Z0(1)
        DY=NODE(LN(I,J),2)-Z0(2)
        VX=Z1(1)-Z1(3)*DY
        VY=Z1(2)+Z1(3)*DX
        ACCMO(LN(I,J))=0.D0
        
    END DO
  END DO

      RETURN
    END    
!********************************************************************
SUBROUTINE BOUNDT(WMKTYP,NDOF,NPL,NFL,NNODE,NELM,NELEM,ME,NS,LN,BUNTYP,&
                &NODE,NORM,PHIT,PPHIT,DPDS,JCB,ACC,ACCMO,Z0,Z1,Z2,WC)
!********************************************************************
IMPLICIT NONE
INTEGER I,J,K,N,WMKTYP,NDOF,NPL,NFL,NNODE,NELM,NS(NPL),ME(NPL),NELEM(NPL),BUNTYP(NPL),LN(NELM,2)
REAL*8 AX,AY,DX,DY,R,ACC,WC
REAL*8 NODE(NNODE,2),NORM(NELM,2),JCB(NELM),PHIT(NNODE),PPHIT(NNODE)
REAL*8 DPDS(NNODE),DPDSS(NNODE),ACCMO(NNODE)
REAL*8 Z0(3*NDOF),Z1(3*NDOF),Z2(3*NDOF)

!---PREPARE RADIATION CONDITION
    DO I=ME(NFL+2)+1,ME(NFL+3)
    DO J=1,2
		IF(LN(I,J).EQ.NS(NFL+3)) THEN
		  DPDSS(LN(I,J))=0.D0 ! ASSUME A FLAT GROUND SO VN=0
		ELSE
		  DPDSS(LN(I,J))=0.D0 ! DPDSS=0 FOR 2D LINEAR BEM !0.5D0*(-DPDS(LN(I,1))+DPDS(LN(I,2)))/JCB(I)
		END IF
    END DO
    END DO
    
N=0
DO I=1,NPL ! PLANE NUMBER
    DO J=N+1,ME(I) ! ELEMENT NUMBER
        DO K=1,2   ! LN(J,K) = NODE NUMBER
            IF (BUNTYP(I) .EQ. 1)THEN           ! BC FOR FREE SURFACE
                PHIT(LN(J,K))=PHIT(LN(J,K))
                PPHIT(LN(J,K))=0.D0
            ELSE
                IF (I==(NFL+3))THEN               ! BC FOR RADIATION SIDE
                    PPHIT(LN(J,K))=0.D0 !WC*DPDSS(LN(J,K)) !
                    PHIT(LN(J,K))=0.D0
                ELSE IF (I==NPL)THEN              ! BC FOR WAVE MAKER
                    IF(WMKTYP==1)THEN
                        PPHIT(LN(J,K))=-ACC+ACCMO(LN(J,K)) !0.D0 !WC*DPDSS(LN(J,K)) !
                        PHIT(LN(J,K))=0.D0
                    ELSE IF(WMKTYP==2)THEN
                        PPHIT(LN(J,K))=-ACC**DSQRT((NODE(LN(J,K),1)-NODE(NS(NPL-1),1))**2+&
                                &(NODE(LN(J,K),2)-NODE(NS(NPL-1),2))**2)+ACCMO(LN(J,K)) !0.D0 !WC*DPDSS(LN(J,K)) !
                        PHIT(LN(J,K))=0.D0
                    END IF
                ELSE IF (I>1.AND.I<=(NFL+1))THEN  ! BC FOR RIGID BODY
                    DX=NODE(LN(J,K),1)-Z0(1)
                    DY=NODE(LN(J,K),2)-Z0(2)
                    AX=Z2(1)-Z2(3)*DY-Z1(3)*Z1(3)*DX
                    AY=Z2(2)+Z2(3)*DX-Z1(3)*Z1(3)*DY
                    PPHIT(LN(J,K))=AX*NORM(J,1)+AY*NORM(J,2)+ACCMO(LN(J,K))                 
                    PHIT(LN(J,K))=0.D0
                ELSE                              ! BC FOR BOTTOM
                    PPHIT(LN(J,K))=0.D0
                    PHIT(LN(J,K))=0.D0
                END IF
            END IF
        END DO
    END DO
    N=ME(I)
END DO

RETURN
END
!**********************************************************************
      SUBROUTINE TAYES2(NPL,NFL,NNODE,NELM,NS,ME,LN,NODE,NORM,JCB,PHI,PPHI,PHIT,PPHIT,&
                       &DP,DPDS,DPDT,D2PDT,D2P,DELTTIME,GRAV,ACC,Z2)
!**********************************************************************
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER  NPL,NFL,NNODE,NELM,NS(NPL),ME(NPL),LN(NELM,2)
      REAL*8 DELTTIME,GRAV,DPTDS,DPDNDS,Z2(3)
      REAL*8 NODE(NNODE,2),NORM(NELM,2),JCB(NELM)
      REAL*8 PHI(NNODE),PPHI(NNODE),PHIT(NNODE),PPHIT(NNODE)
      REAL*8 DP(NNODE,2),DPDS(NNODE),DPDT(NNODE),D2PDT(NNODE),D2P(NNODE,2)

!********************* ON FREE SURFACE 1 *********************
	DO I=1,ME(1)
	  DO J=1,2
		IF(LN(I,J).EQ.1)THEN
            D2P(LN(I,J),1)=ACC
            DPDNDS=(-0.5*PPHI(LN(I,1))+0.5*PPHI(LN(I,2)))/JCB(I)
            DPTDS=(D2P(LN(I,J),1)-PPHI(1)*DPDNDS*NORM(I,2)-(DPDS(1)*DPDNDS+PPHIT(1))*NORM(1,1))/NORM(1,2)
            D2P(LN(I,J),2)=(PPHIT(LN(I,J))+DPDS(LN(I,J))*DPDNDS)*NORM(I,2)-(DPTDS+PPHI(LN(I,J))*DPDNDS)*NORM(I,1)
        ELSE IF (LN(I,J).EQ.NS(1))THEN
!+++++++++++++++++++++++may be revised for continu
            DPTDS=(-0.5*PHIT(LN(I,1))+0.5*PHIT(LN(I,2)))/JCB(I)
			DPDNDS=(-0.5*PPHI(LN(I,1))+0.5*PPHI(LN(I,2)))/JCB(I)
			D2P(LN(I,J),1)=(DPTDS+PPHI(LN(I,J))*DPDNDS)*NORM(I,2)-(-DPDS(LN(I,J))*DPDNDS-PPHIT(LN(I,J)))*NORM(I,1)
			D2P(LN(I,J),2)=(PPHIT(LN(I,J))+DPDS(LN(I,J))*DPDNDS)*NORM(I,2)-(DPTDS+PPHI(LN(I,J))*DPDNDS)*NORM(I,1)
		ELSE
			DPTDS=(-0.5*PHIT(LN(I,1))+0.5*PHIT(LN(I,2)))/JCB(I)
			DPDNDS=(-0.5*PPHI(LN(I,1))+0.5*PPHI(LN(I,2)))/JCB(I)
			D2P(LN(I,J),1)=(DPTDS+PPHI(LN(I,J))*DPDNDS)*NORM(I,2)-(-DPDS(LN(I,J))*DPDNDS-PPHIT(LN(I,J)))*NORM(I,1)
			D2P(LN(I,J),2)=(PPHIT(LN(I,J))+DPDS(LN(I,J))*DPDNDS)*NORM(I,2)-(DPTDS+PPHI(LN(I,J))*DPDNDS)*NORM(I,1)
		END IF
		D2PDT(LN(I,J))=DP(LN(I,J),1)*D2P(LN(I,J),1)+DP(LN(I,J),2)*D2P(LN(I,J),2)-GRAV*DP(LN(I,J),2)
	  END DO
    END DO

!********************* ON FREE SURFACE 2 *********************
	DO I=ME(NFL+1)+1,ME(NFL+2)
	  DO J=1,2
		IF(LN(I,J).EQ.NS(NFL+1)+1)THEN
!+++++++++++++++++++++++may be revised for continu
			DPTDS=(-0.5*PHIT(LN(I,1))+0.5*PHIT(LN(I,2)))/JCB(I)
			DPDNDS=(-0.5*PPHI(LN(I,1))+0.5*PPHI(LN(I,2)))/JCB(I)
			D2P(LN(I,J),1)=(DPTDS+PPHI(LN(I,J))*DPDNDS)*NORM(I,2)-(-DPDS(LN(I,J))*DPDNDS-PPHIT(LN(I,J)))*NORM(I,1)
			D2P(LN(I,J),2)=(PPHIT(LN(I,J))+DPDS(LN(I,J))*DPDNDS)*NORM(I,2)-(DPTDS+PPHI(LN(I,J))*DPDNDS)*NORM(I,1)
        ELSE IF(LN(I,J).EQ.NS(NFL+2))THEN
            D2P(LN(I,J),1)=0.D0 !PPHIT(NS(NFL+2)+1)+ACCMO(NS(NFL+2)+1)
            DPDNDS=(-0.5*PPHI(LN(I,1))+0.5*PPHI(LN(I,2)))/JCB(I)
            DPTDS=(D2P(LN(I,J),1)-PPHI(1)*DPDNDS*NORM(I,2)-(DPDS(1)*DPDNDS+PPHIT(1))*NORM(1,1))/NORM(1,2)
            D2P(LN(I,J),2)=(PPHIT(LN(I,J))+DPDS(LN(I,J))*DPDNDS)*NORM(I,2)-(DPTDS+PPHI(LN(I,J))*DPDNDS)*NORM(I,1)    
        ELSE
			DPTDS=(-0.5*PHIT(LN(I,1))+0.5*PHIT(LN(I,2)))/JCB(I)
			DPDNDS=(-0.5*PPHI(LN(I,1))+0.5*PPHI(LN(I,2)))/JCB(I)
			D2P(LN(I,J),1)=(DPTDS+PPHI(LN(I,J))*DPDNDS)*NORM(I,2)-(-DPDS(LN(I,J))*DPDNDS-PPHIT(LN(I,J)))*NORM(I,1)
			D2P(LN(I,J),2)=(PPHIT(LN(I,J))+DPDS(LN(I,J))*DPDNDS)*NORM(I,2)-(DPTDS+PPHI(LN(I,J))*DPDNDS)*NORM(I,1)
		END IF
		D2PDT(LN(I,J))=DP(LN(I,J),1)*D2P(LN(I,J),1)+DP(LN(I,J),2)*D2P(LN(I,J),2)-GRAV*DP(LN(I,J),2)
	  END DO
    END DO

      RETURN
      END
!********************************************************************
SUBROUTINE PRESSURE(ICON,THO,GRAV,DEP,NPL,NFL,NNODE,NS,NODE,PHIT,DP,PR,P_ATM)
!********************************************************************
      IMPLICIT NONE
      INTEGER  I,ICON,NPL,NFL,NNODE,NS(NPL)
      REAL*8 DEP,THO,GRAV,P1,P2,P3,P_ATM
      REAL*8 NODE(NNODE,2),PHIT(NNODE),DP(NNODE,2),PR(NNODE)
	  REAL*8 CP1(NS(NFL+2)),CP2(NS(NFL+2)),CP3(NS(NFL+2)),CP(NS(NFL+2))

!----ATOM PRESSURE (BERNOULLI CONSTANT) ON THE FREE SURFACE
	DO I=ICON,ICON !1,NS(1) !
	CP1(I)=THO*PHIT(I)
	CP2(I)=THO*0.5D0*(DP(I,1)**2+DP(I,2)**2)
	CP3(I)=THO*GRAV*(NODE(I,2)-DEP) !NODE(I,2) !
	CP(I)=CP1(I)+CP2(I)+CP3(I)
    ENDDO
!	DO I=NS(1+NFL)+1,NS(2+NFL) !ICON,ICON !
!	CP1(I)=THO*PHIT(I)
!	CP2(I)=THO*0.5D0*(DP(I,1)**2+DP(I,2)**2)
!	CP3(I)=THO*GRAV*(NODE(I,2)-DEP)
!	CP(I)=CP1(I)+CP2(I)+CP3(I)
!    ENDDO   
!    WRITE(24,"(500(E15.8,1X))") TIME,CP(1:NS(1)),CP(NS(1+NFL)+1:NS(2+NFL))
    
	P_ATM=CP(ICON)

!----PRESSURE ON FLOATING BODY
	DO I=NS(1)+1,NS(NFL+1)
	P1=THO*PHIT(I)
	P2=THO*0.5D0*(DP(I,1)**2+DP(I,2)**2)
	P3=THO*GRAV*(NODE(I,2)-DEP) !NODE(I,2) !
	PR(I)=P_ATM-(P1+P2+P3)
    END DO
    
!----PRESSURE ON WALL AND BOTTOM
	DO I=NS(NFL+2)+1,NNODE
	P1=THO*PHIT(I)
	P2=THO*0.5D0*(DP(I,1)**2+DP(I,2)**2)
	P3=THO*GRAV*(NODE(I,2)-DEP) !NODE(I,2) !
	PR(I)=P_ATM-(P1+P2+P3)
	END DO

    RETURN
    END
!********************************************************************
	SUBROUTINE FORCE_WAVE(TIME,GRAV,MB,NPL,NFL,NNODE,NELM,ME,LN,NODE,NORM,JCB,PR,Z0,FOR)
!********************************************************************
	IMPLICIT NONE
	INTEGER I,NPL,NFL,NNODE,NELM,ME(NPL),LN(NELM,2)
	REAL*8 P1(2),P2(2),R1(2),R2(2),M1,M2
	REAL*8 TIME,GRAV,MB,Z0(3),FOR(3)
	REAL*8 NODE(NNODE,2),NORM(NELM,2),JCB(NELM),PR(NNODE)

	FOR=0.D0
	DO I=ME(1)+1,ME(1+NFL)
!--- FORCE ALONG X- AND Y-DIRECTION
        P1(1)=PR(LN(I,1))*NORM(I,1)
        P1(2)=PR(LN(I,1))*NORM(I,2)
        P2(1)=PR(LN(I,2))*NORM(I,1)
        P2(2)=PR(LN(I,2))*NORM(I,2)     
		FOR(1)=FOR(1)+JCB(I)*(P1(1)+P2(1)) ! NOTE JCB=0.5*LENG
		FOR(2)=FOR(2)+JCB(I)*(P1(2)+P2(2))
!--- MOMENT ALONG Z-DIRECTION
		R1(1)=NODE(LN(I,1),1)-Z0(1)
		R1(2)=NODE(LN(I,1),2)-Z0(2)
		R2(1)=NODE(LN(I,2),1)-Z0(1)
		R2(2)=NODE(LN(I,2),2)-Z0(2)
        M1=R1(1)*P1(2)-R1(2)*P1(1)      ! moment = p(r x n)
        M2=R2(1)*P2(2)-R2(2)*P2(1)
		FOR(3)=FOR(3)+JCB(I)*(M1+M2)
    END DO
    
	FOR(2)=FOR(2)-MB*GRAV
    
	RETURN
    END
!********************************************************************
	SUBROUTINE STR(IDOF,GAMMA,BETA,NDOF,TIME,DELTTIME,MB,IB,G_COOR,FOR2,FOR,Z0,Z1,Z2,Z0_PC,Z1_PC,Z2_PC)
!********************************************************************
! THIS SIMPLE CALCULATION IS ONLY AVAILABLE WHEN NDOF = 1     
	IMPLICIT NONE
	INTEGER NDOF,NNODE,NMTYP,IDOF(3)
    REAL*8 TIME,GAMMA,BETA,DELTTIME,MB,IB,FOR(3),FOR2(3),G_COOR(2)
	REAL*8 Z0(3*NDOF),Z1(3*NDOF),Z2(3*NDOF)
    REAL*8 Z0_PC(3*NDOF),Z1_PC(3*NDOF),Z2_PC(3*NDOF)
        
!--- SOLVE TRANSLATIONAL MOTION    
    IF (IDOF(1)==1)THEN
        Z2_PC(1)=(FOR(1)+FOR2(1))/MB
        Z1_PC(1)=Z1(1)+DELTTIME*((1-GAMMA)*Z2(1)+GAMMA*Z2_PC(1))
        Z0_PC(1)=Z0(1)+DELTTIME*Z1(1)+DELTTIME**2*((0.5D0-BETA)*Z2(1)+BETA*Z2_PC(1))
    ELSE
        Z2_PC(1)=Z2(1)
        Z1_PC(1)=Z1(1)
        Z0_PC(1)=Z0(1)
    END IF
    
    IF (IDOF(2)==1)THEN
        Z2_PC(2)=(FOR(2)+FOR2(2))/MB
        Z1_PC(2)=Z1(2)+DELTTIME*((1-GAMMA)*Z2(2)+GAMMA*Z2_PC(2))
        Z0_PC(2)=Z0(2)+DELTTIME*Z1(2)+DELTTIME**2*((0.5D0-BETA)*Z2(2)+BETA*Z2_PC(2))
    ELSE
        Z2_PC(2)=Z2(2)
        Z1_PC(2)=Z1(2)
        Z0_PC(2)=Z0(2)
    END IF
    
!--- SOLVE ROTATIONAL MOTION
    IF (IDOF(3)==1)THEN
        Z2_PC(3)=(FOR(3)+FOR2(3))/IB
        Z1_PC(3)=Z1(3)+DELTTIME*((1-GAMMA)*Z2(3)+GAMMA*Z2_PC(3))
        Z0_PC(3)=Z0(3)+DELTTIME*Z1(3)+DELTTIME**2*((0.5D0-BETA)*Z2(3)+BETA*Z2_PC(3))
    ELSE
        Z2_PC(3)=Z2(3)
        Z1_PC(3)=Z1(3)
        Z0_PC(3)=Z0(3)
    END IF

	RETURN
	END
!******************************************************
      SUBROUTINE DOMAIN(NPL,NGA,NFIELD,NNODE,NELM,NELEM,NS,LN,NODE,NORM,JCB,PHI,PPHI,PHIT,PPHIT,SHA1,SHA2,SH,WT,THO,GRAV,DEP,P_ATM,DP,PR)
!********************************************************************
      IMPLICIT NONE
      INTEGER  I,J,K,L,M,IL,IR,NPL,NGA,NFIELD,NNODE,NELM,NELEM(NPL),NS(NPL),LN(NELM,2)
	  REAL*8 THO,GRAV,DEP,P_ATM
	  REAL*8 HB,DX,DY,PI2,TEMP,P1,P2,P3
	  REAL*8 NODE(NNODE,2),NORM(NELM,2),JCB(NELM)
	  REAL*8 PHI(NNODE),PPHI(NNODE),PHIT(NNODE),PPHIT(NNODE),DP(NNODE,2),PR(NNODE)
	  REAL*8 DNODE(NFIELD,2),DVX(NFIELD),DVY(NFIELD),DPHIT(NFIELD),DPR(NFIELD)
	  REAL*8 KER1(NFIELD,NNODE),KER2(NFIELD,NNODE)
      REAL*8 H(2),G(2),XFUNC(10),YFUNC(10),PXI1(2)
	  REAL*8 WT(NGA),SHA1(NGA),SHA2(NGA),SH(2,NGA)
	
	PI2=2.D0*DACOS(-1.D0)

!----CREATE DOMAIN POINT
	L=1
	DO I=2,NS(1)-1
	  CALL BWLOC(NODE(I,1),NS(NPL-1)-NS(2),NODE(NS(2)+1:NS(NPL-1),1),NS(2),IL,IR)
	  HB=0.5D0*(NODE(IL,2)+NODE(IR,2))+0.015D0 ! keep it a little far away from the boundary
	  DY=-NODE(I,2)/NELEM(2)
		DO J=2,NELEM(2)
		  TEMP=NODE(I,2)+DY*(J-1)
			IF(TEMP>HB)THEN
			DNODE(L,1)=NODE(I,1)
			DNODE(L,2)=NODE(I,2)+DY*(J-1)
			L=L+1
			END IF
		END DO
	END DO

!--- SET A DUMMY NODE TO USELESS DNODE
	DO I=L,NFIELD
	DNODE(I,:)=DNODE(1,:)
	END DO

	KER1=0.D0
	KER2=0.D0
	DVX=0.D0
!----CALCULATE X VELOCITY BY BIE
	DO I = 1,NFIELD
       DO J=1,NELM
		 DO M=1,NGA
			XFUNC(M)=SHA1(M)*NODE(LN(J,1),1)+SHA2(M)*NODE(LN(J,2),1)
			YFUNC(M)=SHA1(M)*NODE(LN(J,1),2)+SHA2(M)*NODE(LN(J,2),2)
		 END DO
		DO K=1,2
         G(K)=0.D0
         H(K)=0.D0
          DO M=1,NGA
			TEMP=JCB(J)*SH(K,M)*WT(M)
            H(K)=H(K)+(((YFUNC(M)-DNODE(I,2))**2-(XFUNC(M)-DNODE(I,1))**2)*NORM(J,1)-2.D0*(YFUNC(M)-DNODE(I,2))*(XFUNC(M)-DNODE(I,1))*NORM(J,2))/((XFUNC(M)-DNODE(I,1))**2+(YFUNC(M)-DNODE(I,2))**2)**2*TEMP
            G(K)=G(K)+(XFUNC(M)-DNODE(I,1))/((XFUNC(M)-DNODE(I,1))**2+(YFUNC(M)-DNODE(I,2))**2)*TEMP
         END DO
	     KER1(I,LN(J,K))=KER1(I,LN(J,K))+H(K)
         KER2(I,LN(J,K))=KER2(I,LN(J,K))+G(K)
		END DO
       END DO
      END DO
	DVX=(MATMUL(KER2,PPHI)-MATMUL(KER1,PHI))/PI2

	KER1=0.D0
	KER2=0.D0
	DVY=0.D0
!----CALCULATE Y VELOCITY BY BIE
	DO I = 1,NFIELD
       DO J=1,NELM
		 DO M=1,NGA
			XFUNC(M)=SHA1(M)*NODE(LN(J,1),1)+SHA2(M)*NODE(LN(J,2),1)
			YFUNC(M)=SHA1(M)*NODE(LN(J,1),2)+SHA2(M)*NODE(LN(J,2),2)
		 END DO
		DO K=1,2
         G(K)=0.D0
         H(K)=0.D0
          DO M=1,NGA
			TEMP=JCB(J)*SH(K,M)*WT(M)
            H(K)=H(K)+(((XFUNC(M)-DNODE(I,1))**2-(YFUNC(M)-DNODE(I,2))**2)*NORM(J,2)-2.D0*(YFUNC(M)-DNODE(I,2))*(XFUNC(M)-DNODE(I,1))*NORM(J,1))/((XFUNC(M)-DNODE(I,1))**2+(YFUNC(M)-DNODE(I,2))**2)**2*TEMP
            G(K)=G(K)+(YFUNC(M)-DNODE(I,2))/((XFUNC(M)-DNODE(I,1))**2+(YFUNC(M)-DNODE(I,2))**2)*TEMP
         END DO
	     KER1(I,LN(J,K))=KER1(I,LN(J,K))+H(K)
         KER2(I,LN(J,K))=KER2(I,LN(J,K))+G(K)
		END DO
       END DO
      END DO
	DVY=(MATMUL(KER2,PPHI)-MATMUL(KER1,PHI))/PI2
	
	KER1=0.D0
	KER2=0.D0
	DPHIT=0.D0
!----CALCULATE PARTIAL POTENTIAL OVER TIME BY BIE
      DO I = 1,NFIELD
       DO J=1,NELM
		 DO M=1,NGA
			XFUNC(M)=SHA1(M)*NODE(LN(J,1),1)+SHA2(M)*NODE(LN(J,2),1)
			YFUNC(M)=SHA1(M)*NODE(LN(J,1),2)+SHA2(M)*NODE(LN(J,2),2)
		 END DO
		DO K=1,2
         G(K)=0.D0
         H(K)=0.D0
          DO M=1,NGA
			TEMP=JCB(J)*SH(K,M)*WT(M)
            H(K)=H(K)+(-1)/((XFUNC(M)-DNODE(I,1))**2+(YFUNC(M)-DNODE(I,2))**2)*((XFUNC(M)-DNODE(I,1))*NORM(J,1)+(YFUNC(M)-DNODE(I,2))*NORM(J,2))*TEMP
            G(K)=G(K)+DLOG(1./((XFUNC(M)-DNODE(I,1))**2+(YFUNC(M)-DNODE(I,2))**2)**0.5)*TEMP	 
		 END DO
	     KER1(I,LN(J,K))=KER1(I,LN(J,K))+H(K)
         KER2(I,LN(J,K))=KER2(I,LN(J,K))+G(K)
		END DO
       END DO
      END DO
	DPHIT=(MATMUL(KER2,PPHIT)-MATMUL(KER1,PHIT))/PI2

!----CALCULATE PRESSURE DISTRIBUTION IN DOMAIN
	DO I=1,NFIELD
	P1=THO*DPHIT(I)
	P2=THO*0.5*(DVX(I)**2+DVY(I)**2)
	P3=THO*GRAV*(DNODE(I,2)-DEP)
	DPR(I)=P_ATM-(P1+P2+P3)
	END DO

	WRITE(12,'(5000(1X,F15.7))') NODE(:,1),DNODE(:,1)
	WRITE(12,'(5000(1X,F15.7))') NODE(:,2),DNODE(:,2)
	WRITE(12,'(5000(1X,F15.7))') DP(:,1),DVX
	WRITE(12,'(5000(1X,F15.7))') DP(:,2),DVY
	WRITE(12,'(5000(1X,F15.7))') PR,DPR

      RETURN
      END
!********************************************************************
      SUBROUTINE GAUSS(WT,RT,NGA)
!********************************************************************
      INTEGER NGA
      REAL*8 WT(NGA),RT(NGA)

      SELECT CASE(NGA)
       CASE(3)
        WT(1)=0.55555555
        WT(2)=0.88888889
        WT(3)=0.55555555
        RT(1)=0.77459667
        RT(2)=0.D0
        RT(3)=-0.77459667
       CASE(4)
        WT(1)=0.65214515
        WT(2)=0.34785484
        WT(3)=0.34785484
        WT(4)=0.65214515
        RT(1)=0.33998104
        RT(2)=0.86113631
        RT(3)=-0.86113631
        RT(4)=-0.33998104
       CASE(5)
        WT(1)=0.23692689
        WT(2)=0.47862867
        WT(3)=0.56888889
        WT(4)=0.47862867
        WT(5)=0.23692689
        RT(1)=0.90617985
        RT(2)=0.53846931
        RT(3)=0.D0
        RT(4)=-0.53846931
        RT(5)=-0.90617985
	 CASE(6)
	  WT(1)=0.17132449
	  WT(2)=0.36076157
	  WT(3)=0.46791393
	  WT(4)=0.46791393
	  WT(5)=0.36076157
	  WT(6)=0.17132449
	  RT(1)=0.93246951
	  RT(2)=0.66120938
	  RT(3)=0.23861918
	  RT(4)=-0.23861918
	  RT(5)=-0.66120938
	  RT(6)=-0.9346951
       CASE(8)
        WT(1)=0.1012285362903763D0
        WT(2)=0.2223810344533745D0	
        WT(3)=0.3137066458778873D0
        WT(4)=0.3626837833783620D0
        
        WT(8)=0.1012285362903763D0
        WT(7)=0.2223810344533745D0
        WT(6)=0.3137066458778873D0
        WT(5)=0.3626837833783620D0
        
        RT(1)=0.9602898564975363D0
        RT(2)=0.7966664774136267D0
        RT(3)=0.5255324099163290D0
        RT(4)=0.1834346424956498D0
        
        RT(8)=-0.9602898564975363D0
        RT(7)=-0.7966664774136267D0
        RT(6)=-0.5255324099163290D0
        RT(5)=-0.1834346424956498D0
       CASE(10)
        WT(1)=0.D06667134
        WT(2)=0.14945134
        WT(3)=0.21908636
        WT(4)=0.26926671
        WT(5)=0.29552422
        WT(10)=0.D06667134
        WT(9)=0.14945134
        WT(8)=0.21908636
        WT(7)=0.26926671
        WT(6)=0.29552422
        RT(1)=0.97390652
        RT(2)=0.86506336
        RT(3)=0.67940956
        RT(4)=0.43339539
        RT(5)=0.14887433
        RT(10)=-0.97390652
        RT(9)=-0.86506336
        RT(8)=-0.67940956
        RT(7)=-0.43339539
        RT(6)=-0.14887433
      END SELECT

      RETURN
      END
!********************************************************************
	SUBROUTINE CONVERGE(N,P1,P2,E1LOC,E1,E2)
!********************************************************************
	IMPLICIT NONE
	INTEGER I,N,E1LOC
	REAL*8 P1(N),P2(N),E1,E2

	E1=MAXVAL(DABS(P1-P2))
	E1LOC=MAXLOC(DABS(P1-P2),1)

	E2=0.D0
	DO I=1,N
	E2=E2+(P1(I)-P2(I))**2
	END DO
	E2=DSQRT(E2/N)

	RETURN
	END
!********************************************************************
	SUBROUTINE BWLOC(PX,N,X,IST,IL,IR)
!********************************************************************
	IMPLICIT NONE
	INTEGER I,N,IST,IL,IR
	REAL*8 PX,X(N)

	DO I=1,N-1
		IF(X(I)>PX.AND.X(I+1)<=PX)THEN
		IR=IST+I
		IL=IST+I+1
		GOTO 777
		END IF
	END DO
777 CONTINUE

	RETURN
	END
!********************************************************************
	SUBROUTINE WAVE_SPD(GRAV,OMEGA,D,C)
!********************************************************************
	IMPLICIT NONE
	INTEGER I
	REAL*8 T,K,K2,GRAV,OMEGA,D,C,PI,F0,F1
	PI=DACOS(-1.D0)

	K=1.D0
	DO I=1,100
	F0=K*DTANH(K*D)-OMEGA**2/GRAV
	F1=DTANH(K*D)+K-K*(DTANH(K*D)**2)
	K2=K-(F0)/(F1)
		IF((K2-K)/K<=0.000001D0) THEN
		GOTO 717
		END IF
	K=K2
	END DO
	717 CONTINUE

	C=DSQRT(GRAV*DTANH(K*D)/K)

	RETURN
	END
!********************************************************************
SUBROUTINE COURANT(TIME,DELTTIME,NNODE,NELM,LN,NODE,DP,JCB)
!********************************************************************
    IMPLICIT NONE
    INTEGER I,J,CFLOC,NNODE,NELM,LN(NELM,2)
    REAL*8 TIME,DELTTIME,U,V,VE
	REAL*8 NODE(NNODE,2),DP(NNODE,2),JCB(NELM),CN(NELM),CFL

      DO I=1,NELM
        U=DSQRT(DP(LN(I,1),1)**2+DP(LN(I,1),2)**2)
        V=DSQRT(DP(LN(I,2),1)**2+DP(LN(I,2),2)**2)
        VE=MAX(U,V)
        CN(I)=0.5D0*VE*DELTTIME/JCB(I)
      END DO
      CFL=MAXVAL(CN)
      CFLOC=MAXLOC(CN,1)
  
	!WRITE(23,*) TIME,CFL

    IF (CFL>=250.D0)THEN
    WRITE(20,*) "CFL TOO HIGH"
    WRITE(20,*) "TIME=",TIME
    WRITE(20,*) "CFL=",CFL
    WRITE(20,*) "ELEMENT #",CFLOC
    STOP
    END IF    

	RETURN
    END
!**********************************************************************
       SUBROUTINE SOLVE_BACK(NPL,PHI,PPHI,KER1,KER2,NNODE,NELEM,BUNTYP)
!**********************************************************************
!      TO SOLVE KER1*PHI=KER2*PPHI
!      PPHI=PARTIAL PHI OVER PARTIAL N
!      USING GAUSSIAN ELIMINATION WITH BACKSUBSTITUTION
!======================================================================
       IMPLICIT NONE
       INTEGER I,J,K,L,N,NPL,NNODE,NELEM(NPL),BUNTYP(NPL)
       REAL*8 PHI(NNODE),PPHI(NNODE)
       REAL*8 KER1(NNODE,NNODE),KER2(NNODE,NNODE)
       REAL*8 H1(NNODE,NNODE),Q1(NNODE),TEMP(NNODE)
       REAL*8 SUM,A,SIG,G1(NNODE,NNODE+1),P1(NNODE)

!-----PHI PPHI----         
       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,(NELEM(I)+1)+N
            IF (BUNTYP(I) .EQ. 1) THEN
               Q1(J)=PHI(J)
            ELSE
               Q1(J)=PPHI(J)
            END IF
          END DO
          N=N+(NELEM(I)+1)
       END DO
!-------------------
         DO I=1,NNODE
           N=0
           DO L=1,NPL
             DO J=K+N,(NELEM(L)+1)+N
             IF (BUNTYP(L) .EQ. 1) THEN
               H1(I,J)=KER1(I,J)  !G*UNKNOWN=H*KNOWN
               G1(I,J)=KER2(I,J)
              ELSE
               H1(I,J)=-KER2(I,J)
               G1(I,J)=-KER1(I,J)
             END IF
             END DO
            N=N+(NELEM(L)+1)
          END DO
        END DO

       DO I=1,NNODE
          TEMP(I)=0.0
          DO J=1,NNODE
          TEMP(I)=TEMP(I)+H1(I,J)*Q1(J)
          END DO
       END DO
!*************GAUSSIAN ELIMINATION WITH BACKSUBSTITUTION*********
      DO I=1,NNODE
        G1(I,NNODE+1)=TEMP(I)
      END DO
      DO I=1,NNODE
         SUM=0.0
         DO K=I,NNODE
            IF (G1(K,I) .NE. 0) THEN
               IF (K .NE. I) THEN
               IF (G1(I,I) .EQ. 0) THEN
                 WRITE(*,*) 'SOME OF THE DIAG-TERMS ARE ZERO'
                 WRITE(20,*) 'SOME OF THE DIAG-TERMS ARE ZERO'
                 STOP
               END IF
               A=G1(K,I)/G1(I,I)
               DO J=I,NNODE+1
                  G1(K,J)=G1(K,J)-A*G1(I,J)
               END DO
               END IF
            END IF
            SUM=SUM+G1(K,I)
          END DO
          IF (SUM .EQ. 0.) THEN
          WRITE(*,*) 'NO UNIQUE SOLUTION EXISTS STOP AT GAUSSELI'
          WRITE(20,*) 'NO UNIQUE SOLUTION EXISTS STOP AT GAUSSELI'
          STOP
          END IF
      END DO

      P1(NNODE)=G1(NNODE,NNODE+1)/G1(NNODE,NNODE)
      DO I=NNODE-1,1,-1
         SIG=0.0
         DO J=I+1,NNODE
            SIG=G1(I,J)*P1(J)+SIG
          END DO
         P1(I)=(G1(I,NNODE+1)-SIG)/G1(I,I)
      END DO

       K=1
       N=0
       DO I=1,NPL
          DO J=K+N,(NELEM(I)+1)+N
            IF (BUNTYP(I) .EQ. 1) THEN
               PPHI(J)=P1(J)
            ELSE
               PHI(J)=P1(J)
            END IF
          END DO
          N=N+NELEM(I)+1
       END DO

      RETURN
    END
!********************************************************************
	SUBROUTINE STABLE(TIME,NNODE,NODE,NNODE2,NODE2)
!********************************************************************
	IMPLICIT NONE
	INTEGER NNODE,NNODE2
	REAL*8 TIME,NODE(NNODE),NODE2(NNODE2)

    IF(MAXVAL(DABS(NODE))>=1000000.D0)THEN
        WRITE(20,*) 'WAVE TO HIGH'
        WRITE(*,*) 'WAVE TO HIGH'
        STOP
    ELSE IF(MAXVAL(DABS(NODE2))>=1000000.D0)THEN
        WRITE(20,*) 'PMTLD TO HIGH'
        WRITE(*,*) 'PMTLD TO HIGH'
        STOP
    END IF
    
	RETURN
    END
!**********************************************************************
!   2D PMTLD HAS 3 DOFs MOTIONON, PUT ON A FLOATING PLATFORM
!   BY WEN-HUAI TSAO AT LSU DEC. 2021
!**********************************************************************
SUBROUTINE PMTLD(NNODE,NELM,NGA,ICON,NELEM,ME,NS,LN,BUNTYP,MU1,MU2,POR,&
                &GRAV,MU,THO,DELTTIME,TL,TH,DEP,TLDLOC,COOR,Z0,Z1,Z2,WT,SHA1,SHA2,SH,&
                &NODE,PHI,FOR,KE,PE,TE,WORK,DTDT,DVDT,DEDT,DWDT)
IMPLICIT NONE
!-----PARAMETERS ONLY USED HERE
INTEGER I,J
REAL*8 P_ATM
REAL*8 KER1(NNODE,NNODE),KER2(NNODE,NNODE)
REAL*8 NORM(NELM,2),JCB(NELM),LENG(NELM)
REAL*8 PPHI(NNODE),PHIT(NNODE),PPHIT(NNODE)
REAL*8 DP(NNODE,2),D2PDT(NNODE),D2P(NNODE,2)
REAL*8 DPDS(NNODE),DPDT(NNODE),PR(NNODE),ACCMO(NNODE)
REAL*8 G1(NNODE,NNODE),H1(NNODE,NNODE),EYE(NNODE),NODE_TEMP(NNODE,2)

!-----INPUT PARAMETERS
INTEGER NNODE,NELM,NGA,ICON,NELEM(4),ME(4),NS(4),LN(NELM,2),BUNTYP(4)
REAL*8 MU1,MU2,POR,GRAV,MU,THO,DELTTIME,TL,TH,DEP,TLDLOC(2)
REAL*8 COOR(4,2),Z0(3),Z1(3),Z2(3)
REAL*8 WT(NGA),SHA1(NGA),SHA2(NGA),SH(2,NGA)

!-----OUTPUT PARAMETERS
REAL*8 KE,PE,TE,WORK,DTDT,DVDT,DEDT,DWDT
REAL*8 NODE(NNODE,2),PHI(NNODE),FOR(3)

EYE=1.D0
NODE_TEMP=NODE

!---SET PMTLD BIE
CALL KERNEL(KER1,KER2,NODE,NORM,JCB,LENG,LN,NNODE,NELM,NGA,SHA1,SHA2,SH,WT,EYE)

!---SOLVE 1ST-ORDER TERM
CALL BOUND_TLD(PHI,PPHI,NNODE,NELM,NS,ME,LN,NODE,NORM,Z0,Z1)
CALL SOLVE_LAPACK(4,PHI,PPHI,KER1,KER2,NNODE,NELEM,BUNTYP)
!CALL SOLVE_LAPACK2_1(4,PHI,PPHI,KER1,KER2,G1,H1,NNODE,NELEM,BUNTYP)
CALL TAYES1_TLD(NNODE,NELM,ME,NS,LN,NODE,NORM,JCB,PHI,PPHI,DPDS,DP,PHIT,DPDT,GRAV,MU1,MU2,POR,Z0,Z1,DEP)
 
!---SOLVE 2ND-ORDER TERM
CALL ACCBC_TLD(NNODE,NELM,NELEM,ME,LN,NORM,JCB,PHI,PPHI,DP,ACCMO)
CALL BOUNDT_TLD(NNODE,NELM,NS,ME,LN,NODE,NORM,PHIT,PPHIT,ACCMO,Z0,Z1,Z2)
CALL SOLVE_LAPACK(4,PHIT,PPHIT,KER1,KER2,NNODE,NELEM,BUNTYP)
!CALL SOLVE_LAPACK2_2(4,PHIT,PPHIT,G1,H1,NNODE,NELEM,BUNTYP)
CALL TAYES2_TLD(NNODE,NELM,NS,ME,LN,NODE,NORM,JCB,PHI,PPHI,PHIT,PPHIT,&
                &DP,DPDS,DPDT,D2PDT,D2P,DELTTIME,GRAV,Z2)

!---CALCULATE PRESSURE AND FORCE TO THE BODY
CALL PRESSURE_TLD(ICON,DEP,THO,GRAV,NNODE,NS,NODE,PHIT,DP,PR,P_ATM)
CALL FORCE_TLD(NNODE,NELM,ME,LN,NODE,NORM,JCB,PR,Z0,POR,FOR)

!---REMESH FREE SURFACE OF PMTLD
CALL REMESH_TLD(NNODE,NELM,NELEM,ME,NS,NODE,NORM,DP,D2P,PHI,DPDT,D2PDT,DELTTIME,Z0,TL,TLDLOC)

!---CALCULATE ENERGY OF PMTLD
CALL ENERGY_TLD(NNODE,NELEM,NELM,NS,LN,NGA,NODE,PHI,PPHI,PHIT,PPHIT,DP,PR,NODE_TEMP,&
                &JCB,WT,SH,SHA1,SHA2,DEP,GRAV,THO,Z0,KE,PE,TE,WORK,DTDT,DVDT,DEDT,DWDT)

RETURN
END
!**********************************************************************
SUBROUTINE ENERGY_TLD(NNODE,NELEM,NELM,NS,LN,NGA,NODE,PHI,PPHI,PHIT,PPHIT,DP,PR,NODE_TEMP,&
                    &JCB,WT,SH,SHA1,SHA2,DEP,GRAV,THO,Z0,KE,PE,TE,WORK,DTDT,DVDT,DEDT,DWDT)
!**********************************************************************
    IMPLICIT NONE
    INTEGER I,J,K,M,NNODE,NGA,NELM,NELEM(4),NS(4),LN(NELM,2)
    REAL*8 NODE(NNODE,2),PHI(NNODE),PPHI(NNODE),PHIT(NNODE),PPHIT(NNODE),DP(NNODE,2),PR(NNODE)
    REAL*8 JCB(NELM),WT(NGA),SHA1(NGA),SHA2(NGA),SH(2,NGA),Z0(3)
    REAL*8 XFUNC(NGA),YFUNC(NGA),PFUNC(NGA),PPFUNC(NGA)
    REAL*8 PTFUNC(NGA),PPTFUNC(NGA),DPXFUNC(NGA),DPYFUNC(NGA)
    REAL*8 TEMP,ITA,T1,T2,DN1,DN2,S1,S2,DDIS(NNODE,2),NODE_TEMP(NNODE,2)
    REAL*8 DEP,GRAV,THO,KE,PE,TE,WORK,DTDT,DVDT,DEDT,DWDT

!-----KE AND DT/DT
	KE=0.D0
	DTDT=0.D0
	DO J=1,NELM
		DO M=1,NGA
			PFUNC(M)=SHA1(M)*PHI(LN(J,1))+SHA2(M)*PHI(LN(J,2))
			PPFUNC(M)=SHA1(M)*PPHI(LN(J,1))+SHA2(M)*PPHI(LN(J,2))
			PTFUNC(M)=SHA1(M)*PHIT(LN(J,1))+SHA2(M)*PHIT(LN(J,2))
			PPTFUNC(M)=SHA1(M)*PPHIT(LN(J,1))+SHA2(M)*PPHIT(LN(J,2))
		END DO
		DO K=1,2
			DO M=1,NGA
			T1=PTFUNC(M)*PPFUNC(M)*JCB(J)*SH(K,M)*WT(M)
			T2=PFUNC(M)*PPTFUNC(M)*JCB(J)*SH(K,M)*WT(M)
			KE=KE+PFUNC(M)*PPFUNC(M)*JCB(J)*SH(K,M)*WT(M)
			DTDT=DTDT+T1+T2
			END DO
		END DO
	END DO
	KE=KE*0.5D0*THO
	DTDT=DTDT*0.5D0*THO

!-----PE AND DV/DT
	PE=0.D0
	DVDT=0.D0
	DO J=1,NELEM(1)
		DO M=1,NGA
			XFUNC(M)=SHA1(M)*NODE(LN(J,1),1)+SHA2(M)*NODE(LN(J,2),1)
			YFUNC(M)=SHA1(M)*NODE(LN(J,1),2)+SHA2(M)*NODE(LN(J,2),2)
			DPXFUNC(M)=SHA1(M)*DP(LN(J,1),1)+SHA2(M)*DP(LN(J,2),1)
			DPYFUNC(M)=SHA1(M)*DP(LN(J,1),2)+SHA2(M)*DP(LN(J,2),2)
		END DO
		DO K=1,2
			DO M=1,NGA
			ITA=YFUNC(M)-DEP
			PE=PE+ITA**2*JCB(J)*SH(K,M)*WT(M)
			DVDT=DVDT+2.D0*ITA*DPYFUNC(M)*JCB(J)*SH(K,M)*WT(M)
			END DO
		END DO
	END DO
	PE=PE*0.5D0*THO*GRAV
	DVDT=DVDT*0.5D0*THO*GRAV

!-----TOTAL ENERGY AND ITS TIME DERIVATIVE    
	TE=KE+PE
	DEDT=DTDT+DVDT
    
!-----INPUT WORK AND DWDT=SUM(P*VN)
	DO I=1,NNODE
	DDIS(I,:)=NODE(I,:)-NODE_TEMP(I,:)
    END DO

	S1=DCOS(Z0(3))
	S2=DSIN(Z0(3))
    
DWDT=0.D0
DO I=NS(1),NS(2)-1
	TEMP=SQRT((NODE(I,1)-NODE(I+1,1))**2+(NODE(I,2)-NODE(I+1,2))**2)
	DWDT=DWDT+(PR(I)*PPHI(I)+PR(I+1)*PPHI(I+1))*TEMP/2.D0
	DN1=DDIS(I,1)*S1+DDIS(I,2)*S2
	DN2=DDIS(I+1,1)*S1+DDIS(I+1,2)*S2
	WORK=WORK+(PR(I)*DN1+PR(I+1)*DN2)*TEMP/2.D0
END DO

DO I=NS(2),NS(3)-1
	TEMP=SQRT((NODE(I,1)-NODE(I+1,1))**2+(NODE(I,2)-NODE(I+1,2))**2)
	DWDT=DWDT+(PR(I)*PPHI(I)+PR(I+1)*PPHI(I+1))*TEMP/2.D0
	DN1=DDIS(I,1)*S2-DDIS(I,2)*S1
	DN2=DDIS(I+1,1)*S2-DDIS(I+1,2)*S1
	WORK=WORK+(PR(I)*DN1+PR(I+1)*DN2)*TEMP/2.D0
END DO

DO I=NS(3),NNODE-1
	TEMP=SQRT((NODE(I,1)-NODE(I+1,1))**2+(NODE(I,2)-NODE(I+1,2))**2)
	DWDT=DWDT+(PR(I)*PPHI(I)+PR(I+1)*PPHI(I+1))*TEMP/2.D0
	DN1=-DDIS(I,1)*S1-DDIS(I,2)*S2
	DN2=-DDIS(I+1,1)*S1-DDIS(I+1,2)*S2
	WORK=WORK+(PR(I)*DN1+PR(I+1)*DN2)*TEMP/2.D0
END DO
    


      RETURN
      END
!**********************************************************************
SUBROUTINE REMESH_TLD(NNODE,NELM,NELEM,ME,NS,NODE,NORM,DP,D2P,PHI,DPDT,D2PDT,DELTTIME,Z0,TL,TLDLOC)
!**********************************************************************
      IMPLICIT NONE
      INTEGER I,K,NELM,NNODE,NELEM(4),ME(4),NS(4)
      REAL*8 TL_TEMP,SX,SY,DELT,DELTTIME,TL,TLDLOC(2)
      REAL*8 Z0(3),Q(2,2),CR(2),CL(2),HL,HR,S1,S2
      REAL*8 NODE(NNODE,2),NORM(NELM,2),DP(NNODE,2),D2P(NNODE,2)
      REAL*8 PHI(NNODE),DPDT(NNODE),D2PDT(NNODE)

!-----REMESH FREE SURFACE NODE-----
    DO I=1,NS(1)
	NODE(I,1)=NODE(I,1)+DP(I,1)*DELTTIME+0.5D0*D2P(I,1)*DELTTIME**2
    NODE(I,2)=NODE(I,2)+DP(I,2)*DELTTIME+0.5D0*D2P(I,2)*DELTTIME**2
	PHI(I)=PHI(I)+DELTTIME*DPDT(I)+0.5D0*D2PDT(I)*DELTTIME**2
    END DO

!---CHECK IF THE UPDATED FREE SURFACE IS WIDER THAN TLD
    TL_TEMP=(NODE(NS(1),1)-NODE(1,1))*DCOS(Z0(3))+(NODE(NS(1),2)-NODE(1,2))*DSIN(Z0(3))
    IF(ABS((TL_TEMP-TL)/TL)>0.05D0)THEN
        WRITE(*,*) 'TLD GROW WIDE'
        WRITE(20,*) 'TLD GROW WIDE'
        STOP
    END IF
    
!------SET ROTATIONAL MATRIX
    S1=DCOS(Z0(3))
    S2=-DSIN(Z0(3))
    Q(1,1)=S1
    Q(1,2)=S2
    Q(2,1)=-Q(1,2)
    Q(2,2)=Q(1,1)
    
!------VECTOR OF TLD CORNER TO COG
    CR(1)=TLDLOC(1)+0.5D0*TL
    CR(2)=TLDLOC(2)
    CL(1)=TLDLOC(1)-0.5D0*TL
    CL(2)=TLDLOC(2)
    
    CR=MATMUL(Q,CR)+Z0(1:2)
    CL=MATMUL(Q,CL)+Z0(1:2)

!---CORRECT THE FREE-SURFACE END NODE REMAINS ON THE WALLS
    HL=(NODE(1,1)-CL(1))*S2+(NODE(1,2)-CL(2))*S1
    NODE(1,1)=CL(1)+HL*S2
    NODE(1,2)=CL(2)+HL*S1    
    HR=(NODE(NS(1),1)-CR(1))*S2+(NODE(NS(1),2)-CR(2))*S1
    NODE(NS(1),1)=CR(1)+HR*S2
    NODE(NS(1),2)=CR(2)+HR*S1  
    
!------REMESH FOR TLD RIGHT SIDE
    DELT=DSQRT((NODE(NS(1),1)-CR(1))**2+(NODE(NS(1),2)-CR(2))**2)/NELEM(2)
    SX=-S2
    SY=-S1
    K=1
    DO I=NS(1)+1,NS(2)
        NODE(I,1)=NODE(NS(1),1)+DELT*(K-1)*SX
        NODE(I,2)=NODE(NS(1),2)+DELT*(K-1)*SY
        K=K+1
    END DO

!------REMESH FOR TLD BOTTOM
    DELT=DSQRT((CL(1)-CR(1))**2+(CL(2)-CR(2))**2)/NELEM(3)
    SX=-S1
    SY=S2
    K=1
    DO I=NS(2)+1,NS(3)
        NODE(I,1)=NODE(NS(2),1)+DELT*(K-1)*SX
        NODE(I,2)=NODE(NS(2),2)+DELT*(K-1)*SY
        K=K+1
    END DO

!------REMESH FOR TLD LEFT SIDE
    DELT=DSQRT((NODE(1,1)-CL(1))**2+(NODE(1,2)-CL(2))**2)/NELEM(4)
    SX=S2
    SY=S1
    K=1
    DO I=NS(3)+1,NNODE
        NODE(I,1)=NODE(NS(3),1)+DELT*(K-1)*SX
        NODE(I,2)=NODE(NS(3),2)+DELT*(K-1)*SY
        K=K+1
    END DO
    
    RETURN
    END
!**********************************************************************
    SUBROUTINE BOUND_TLD(PHI,PPHI,NNODE,NELM,NS,ME,LN,NODE,NORM,Z0,Z1)
!**********************************************************************
    IMPLICIT NONE
    INTEGER I,J,K,NNODE,NELM,NS(4),ME(4),LN(NELM,2)
    REAL*8 DX,DY,VX,VY
    REAL*8 NODE(NNODE,2),NORM(NELM,2),PHI(NNODE),PPHI(NNODE),Z0(3),Z1(3)

      DO I=1,NS(1)
         PHI(I)=PHI(I)
         PPHI(I)=0.0
      END DO

      DO I=2,4
          DO J=ME(I-1)+1,ME(I)
              DO K=1,2
                  DX=NODE(LN(J,K),1)-Z0(1)
                  DY=NODE(LN(J,K),2)-Z0(2)
                  VX=Z1(1)-Z1(3)*DY
                  VY=Z1(2)+Z1(3)*DX
                  PPHI(LN(J,K))=VX*NORM(J,1)+VY*NORM(J,2)
                  PHI(LN(J,K))=0.D0     
              END DO
          END DO
      END DO

      RETURN
      END 
!********************************************************************
SUBROUTINE TAYES1_TLD(NNODE,NELM,ME,NS,LN,NODE,NORM,JCB,PHI,PPHI,DPDS,DP,PHIT,DPDT,GRAV,MU1,MU2,POR,Z0,Z1,DEP)
!********************************************************************
    IMPLICIT NONE
    INTEGER I,J,K,NNODE,NELM,ME(4),NS(4),LN(NELM,2)
    REAL*8 GRAV,MU1,MU2,POR,DEP,TEMP1,TEMP2,DX,DY,VX,VY
	REAL*8 NODE(NNODE,2),NORM(NELM,2),JCB(NELM),PHI(NNODE),PPHI(NNODE)
    REAL*8 DP(NNODE,2),DPDS(NNODE),PHIT(NNODE),DPDT(NNODE)
    REAL*8 Z0(3),Z1(3)

	DPDS=0.D0
!********************* ON FREE SURFACE *********************
    DO I=1,ME(1)
    DO J=1,2
		IF(LN(I,J).EQ.1) THEN
!            DP(LN(I,J),1)=-PPHI(NNODE)
!            DPDS(LN(I,J))=(DP(LN(I,J),1)-PPHI(LN(I,J))*NORM(I,1))/NORM(I,2)
!            DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
            TEMP1=PPHI(NNODE)-PPHI(1)*(NORM(I,1)*NORM(NELM,1)+NORM(I,2)*NORM(NELM,2))
            TEMP2=NORM(I,2)*NORM(NELM,1)-NORM(I,1)*NORM(NELM,2)
            DPDS(LN(I,J))=TEMP1/TEMP2
			DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)                       
        ELSE IF(LN(I,J).EQ.NS(1)) THEN
            TEMP1=PPHI(NS(1)+1)-PPHI(NS(1))*(NORM(I,1)*NORM(I+1,1)+NORM(I,2)*NORM(I+1,2))
            TEMP2=NORM(I,2)*NORM(I+1,1)-NORM(I,1)*NORM(I+1,2)
            DPDS(LN(I,J))=TEMP1/TEMP2
			DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		ELSE
			DPDS(LN(I,J))=(-0.5D0*PHI(LN(I,1))+0.5D0*PHI(LN(I,2)))/JCB(I)
			DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		END IF
    END DO
    END DO

DO I=1,NS(1)
!DPDT(I)=0.5D0*(DP(I,1)**2+DP(I,2)**2)-GRAV*NODE(I,2)-MU*PHI(I) !-GRAV*(NODE(I,2)-DEP)-MU*PHI(I) !
DPDT(I)=0.5*(DP(I,1)**2+DP(I,2)**2)-GRAV*(NODE(I,2)-DEP)& !-GRAV*NODE(I,2)&
&-MU1*PHI(I)& !-POR*MU1*(PHI(I)-VEL*NODE(I,1)) ! DARCY'S LAW
&-MU2*SQRT(DP(I,1)**2+DP(I,2)**2)*PHI(I)       ! FORCHHEMER LAW
! (PHI(I)-VEL*NODE(I,1))
PHIT(I)=DPDT(I)-(DP(I,1)**2+DP(I,2)**2)
END DO
    
!*********************ON RIGHT WALL*********************
    DO I=ME(1)+1,ME(2)
    DO J=1,2
		IF(LN(I,J).EQ.NS(1)+1) THEN
!			DPDS(LN(I,J))=-DP(NS(1),2)
!            DP(LN(I,J),1)=PPHI(LN(I,J))
!            DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
            DP(LN(I,J),:)=DP(NS(1),:)
        ELSE IF(LN(I,J).EQ.NS(2)) THEN
!			DPDS(LN(I,J))=0.D0
!            DP(LN(I,J),1)=PPHI(LN(I,J))
!            DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
            DX=NODE(LN(I,J),1)-Z0(1)
            DY=NODE(LN(I,J),2)-Z0(2)
            DP(LN(I,J),1)=Z1(1)-Z1(3)*DY
            DP(LN(I,J),2)=Z1(2)+Z1(3)*DX
		ELSE
			DPDS(LN(I,J))=(-0.5D0*PHI(LN(I,1))+0.5D0*PHI(LN(I,2)))/JCB(I)
			DP(LN(I,J),1)=PPHI(LN(I,J))
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		END IF
    END DO
    END DO

!*********************BOTTOM*********************
    DO I=ME(2)+1,ME(3)
    DO J=1,2
		  DPDS(LN(I,J))=(-0.5D0*PHI(LN(I,1))+0.5D0*PHI(LN(I,2)))/JCB(I)
!		  DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
!		  DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
            DX=NODE(LN(I,J),1)-Z0(1)
            DY=NODE(LN(I,J),2)-Z0(2)
            DP(LN(I,J),1)=Z1(1)-Z1(3)*DY
            DP(LN(I,J),2)=Z1(2)+Z1(3)*DX
    END DO
    END DO
    
!*********************LEFT WALL*********************
    DO I=ME(3)+1,ME(4)
    DO J=1,2
		IF(LN(I,J).EQ.NS(3)+1) THEN
!			DPDS(LN(I,J))=0.D0
!			DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
!			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
            DX=NODE(LN(I,J),1)-Z0(1)
            DY=NODE(LN(I,J),2)-Z0(2)
            DP(LN(I,J),1)=Z1(1)-Z1(3)*DY
            DP(LN(I,J),2)=Z1(2)+Z1(3)*DX            
        ELSE IF(LN(I,J).EQ.NNODE) THEN
!			DPDS(LN(I,J))=DP(1,2)
!			DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
!			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
            DP(LN(I,J),:)=DP(1,:)
		ELSE
			DPDS(LN(I,J))=(-0.5D0*PHI(LN(I,1))+0.5D0*PHI(LN(I,2)))/JCB(I)
			DP(LN(I,J),1)=DPDS(LN(I,J))*NORM(I,2)+PPHI(LN(I,J))*NORM(I,1)
			DP(LN(I,J),2)=-DPDS(LN(I,J))*NORM(I,1)+PPHI(LN(I,J))*NORM(I,2)
		END IF
    END DO
    END DO

    RETURN
    END
!**********************************************************************
SUBROUTINE BOUNDT_TLD(NNODE,NELM,NS,ME,LN,NODE,NORM,PHIT,PPHIT,ACCMO,Z0,Z1,Z2)
!**********************************************************************
       IMPLICIT NONE
       INTEGER I,J,K,NNODE,NELM,NS(4),ME(4),LN(NELM,2)
       REAL*8 DX,DY,AX,AY
       REAL*8 NODE(NNODE,2),NORM(NELM,2),PHIT(NNODE),PPHIT(NNODE),ACCMO(NNODE)
       REAL*8 Z0(3),Z1(3),Z2(3)

      DO I=1,NS(1)
       PHIT(I) =PHIT(I)
       PPHIT(I)=0.0
      END DO

      DO I=2,4
          DO J=ME(I-1)+1,ME(I)
              DO K=1,2
                  DX=NODE(LN(J,K),1)-Z0(1)
                  DY=NODE(LN(J,K),2)-Z0(2)
                  AX=Z2(1)-Z2(3)*DY-Z1(3)*Z1(3)*DX
                  AY=Z2(2)+Z2(3)*DX-Z1(3)*Z1(3)*DY
                  PPHIT(LN(J,K))=AX*NORM(J,1)+AY*NORM(J,2)+ACCMO(LN(J,K))                 
                  PHIT(LN(J,K))=0.D0
              END DO
          END DO
      END DO

      RETURN
      END
!**********************************************************************
      SUBROUTINE TAYES2_TLD(NNODE,NELM,NS,ME,LN,NODE,NORM,JCB,PHI,PPHI,PHIT,PPHIT,&
                       &DP,DPDS,DPDT,D2PDT,D2P,DELTTIME,GRAV,Z2)
!**********************************************************************
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER  NNODE,NELM,NS(4),ME(4),LN(NELM,2)
      REAL*8 DELTTIME,GRAV,DPTDS,DPDNDS,Z2(3)
      REAL*8 NODE(NNODE,2),NORM(NELM,2),JCB(NELM)
      REAL*8 PHI(NNODE),PPHI(NNODE),PHIT(NNODE),PPHIT(NNODE)
      REAL*8 DP(NNODE,2),DPDS(NNODE),DPDT(NNODE),D2PDT(NNODE),D2P(NNODE,2)

!********************* ON FREE SURFACE *********************
	DO I=1,ME(1)
	  DO J=1,2
!		IF(LN(I,J).EQ.1)THEN
!            D2P(LN(I,J),1)=0 !++++++++++++++++++++++++++++++
!            DPDNDS=(-0.5*PPHI(LN(I,1))+0.5*PPHI(LN(I,2)))/JCB(I)
!            DPTDS=(D2P(LN(I,J),1)-PPHI(1)*DPDNDS*NORM(I,2)-(DPDS(1)*DPDNDS+PPHIT(1))*NORM(1,1))/NORM(1,2)
!            D2P(LN(I,J),2)=(PPHIT(LN(I,J))+DPDS(LN(I,J))*DPDNDS)*NORM(I,2)-(DPTDS+PPHI(LN(I,J))*DPDNDS)*NORM(I,1)
            
!        ELSE IF (LN(I,J).EQ.NS(1))THEN !+++++++++++++++++++++++may be revised for continu
!            DPTDS=(-0.5*PHIT(LN(I,1))+0.5*PHIT(LN(I,2)))/JCB(I)
!			DPDNDS=(-0.5*PPHI(LN(I,1))+0.5*PPHI(LN(I,2)))/JCB(I)
!			D2P(LN(I,J),1)=(DPTDS+PPHI(LN(I,J))*DPDNDS)*NORM(I,2)-(-DPDS(LN(I,J))*DPDNDS-PPHIT(LN(I,J)))*NORM(I,1)
!			D2P(LN(I,J),2)=(PPHIT(LN(I,J))+DPDS(LN(I,J))*DPDNDS)*NORM(I,2)-(DPTDS+PPHI(LN(I,J))*DPDNDS)*NORM(I,1)
!		ELSE
			DPTDS=(-0.5*PHIT(LN(I,1))+0.5*PHIT(LN(I,2)))/JCB(I)
			DPDNDS=(-0.5*PPHI(LN(I,1))+0.5*PPHI(LN(I,2)))/JCB(I)
			D2P(LN(I,J),1)=(DPTDS+PPHI(LN(I,J))*DPDNDS)*NORM(I,2)-(-DPDS(LN(I,J))*DPDNDS-PPHIT(LN(I,J)))*NORM(I,1)
			D2P(LN(I,J),2)=(PPHIT(LN(I,J))+DPDS(LN(I,J))*DPDNDS)*NORM(I,2)-(DPTDS+PPHI(LN(I,J))*DPDNDS)*NORM(I,1)
!		END IF
		D2PDT(LN(I,J))=DP(LN(I,J),1)*D2P(LN(I,J),1)+DP(LN(I,J),2)*D2P(LN(I,J),2)-GRAV*DP(LN(I,J),2)
	  END DO
    END DO

      RETURN
    END
!********************************************************************
SUBROUTINE PRESSURE_TLD(ICON,DEP,THO,GRAV,NNODE,NS,NODE,PHIT,DP,PR,P_ATM)
!********************************************************************
      IMPLICIT NONE
      INTEGER  I,ICON,NNODE,NS(4)
      REAL*8 DEP,THO,GRAV,P1,P2,P3,P_ATM
      REAL*8 NODE(NNODE,2),PHIT(NNODE),DP(NNODE,2),PR(NNODE)
	  REAL*8 CP1(NS(1)),CP2(NS(1)),CP3(NS(1)),CP(NS(1))

!----ATOM PRESSURE (BERNOULLI CONSTANT) ON THE FREE SURFACE
	DO I=ICON,ICON !1,NS(1) !
	CP1(I)=THO*PHIT(I)
	CP2(I)=THO*0.5D0*(DP(I,1)**2+DP(I,2)**2)
	CP3(I)=THO*GRAV*(NODE(I,2)-DEP) !NODE(I,2) !
	CP(I)=CP1(I)+CP2(I)+CP3(I)
    END DO
!    WRITE(24,"(500(E15.8,1X))") TIME,CP(1:NS(1)),CP(NS(1+NFL)+1:NS(2+NFL))
    
	P_ATM=CP(ICON)

!----PRESSURE ON PMTLD
	DO I=NS(1)+1,NNODE
	P1=THO*PHIT(I)
	P2=THO*0.5D0*(DP(I,1)**2+DP(I,2)**2)
	P3=THO*GRAV*(NODE(I,2)-DEP) !NODE(I,2) !
	PR(I)=P_ATM-(P1+P2+P3)
    !write(*,*) P1,P2,P3
    END DO

    RETURN
    END
!********************************************************************
	SUBROUTINE FORCE_TLD(NNODE,NELM,ME,LN,NODE,NORM,JCB,PR,Z0,POR,FOR)
!********************************************************************
	IMPLICIT NONE
	INTEGER I,NPL,NFL,NNODE,NELM,ME(4),LN(NELM,2)
	REAL*8 POR,M1,M2,P1(2),P2(2),R1(2),R2(2)
	REAL*8 Z0(3),FOR(3)
	REAL*8 NODE(NNODE,2),NORM(NELM,2),JCB(NELM),PR(NNODE)

	FOR=0.D0
	DO I=ME(1)+1,ME(4)
!--- FORCE ALONG X- AND Y-DIRECTION
        P1(1)=PR(LN(I,1))*NORM(I,1)
        P1(2)=PR(LN(I,1))*NORM(I,2)
        P2(1)=PR(LN(I,2))*NORM(I,1)
        P2(2)=PR(LN(I,2))*NORM(I,2)     
		FOR(1)=FOR(1)+POR*JCB(I)*(P1(1)+P2(1)) ! NOTE JCB=0.5*LENG
		FOR(2)=FOR(2)+POR*JCB(I)*(P1(2)+P2(2))
!--- MOMENT ALONG Z-DIRECTION
		R1(1)=NODE(LN(I,1),1)-Z0(1)
		R1(2)=NODE(LN(I,1),2)-Z0(2)
		R2(1)=NODE(LN(I,2),1)-Z0(1)
		R2(2)=NODE(LN(I,2),2)-Z0(2)
        M1=R1(1)*P1(2)-R1(2)*P1(1)
        M2=R2(1)*P2(2)-R2(2)*P2(1)
		FOR(3)=FOR(3)+POR*JCB(I)*(M1+M2)
    END DO
    
	RETURN
    END
!********************************************************************
      SUBROUTINE ACCBC_TLD(NNODE,NELM,NELEM,ME,LN,NORM,JCB,PHI,PPHI,DP,ACCMO)
!********************************************************************
!-------------------------------------------
! THIS SUBROUTINE USE KO 2007 EQUATION
!-------------------------------------------
    IMPLICIT NONE
    INTEGER I,J,K,N,NPL,NNODE,NELM,NELEM(4),ME(4),LN(NELM,2)
    REAL*8 DPNDX,DPNDY,DPNDS,DPDS
	REAL*8 NORM(NELM,2),JCB(NELM),PHI(NNODE),PPHI(NNODE),DP(NNODE,2),ACCMO(NNODE)
N=0
DO K=1,4
  DO I=N+1,ME(K)
    DO J=1,2
		DPNDS=(-0.5D0*PPHI(LN(I,1))+0.5D0*PPHI(LN(I,2)))/JCB(I)
		DPDS=(-0.5D0*PHI(LN(I,1))+0.5D0*PHI(LN(I,2)))/JCB(I)
        ACCMO(LN(I,J))=-DPNDS*DPDS ! PHINN=PHISS=0 FOR LINEAR ELEMENT
!		DPNDX=DPNDS*NORM(I,2)	
!		DPNDY=-DPNDS*NORM(I,1)
!		ACCMO(LN(I,J))=DP(LN(I,J),1)*DPNDX+DP(LN(I,J),2)*DPNDY
    END DO
  END DO
  N=ME(K)
END DO

      RETURN
    END
subroutine get_walltime(wctime)
  use iso_fortran_env, only: int64
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp) :: wctime
  integer(int64) :: r, c
  call system_clock(c, r)
  wctime = real(c, dp) / r
end subroutine get_walltime