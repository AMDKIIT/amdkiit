!This code is taken from CPMD distributed under MIT License
!(https://github.com/CPMD-code)
!
!
!=== ==================================================================
SUBROUTINE SBT(F,G,L,RHOMIN,KAPMIN,H,N,SAVED,XA,TA,TRBE,WW,NDIM,DISC)
!     ==--------------------------------------------------------------==
!  THIS PROGRAM CALCULATES THE SPHERICAL BESSEL TRANSFORM OF ORDER      
!  L FOR A FUNCTION WHOSE VALUES ARE GIVEN IN THE ARRAY F FOR 2**N      
!  VALUES OF R.  THE VALUES OF R ARE GIVEN BY R(I)=EXP(RHOMIN+(I-1)*H)  
!  WHERE I=1,2,...,NT AND NT=2**N.  THE VALUES OF THE TRANSFORM         
!  ARE RETURNED IN THE ARRAY G FOR NT K VALUES.  THE VALUES OF K        
!  ARE EXP(KAPMIN+(I-1)*H), I=1,2,...,NT.                               
!                                                                       
!  TWO DIFFERENT SETS OF RESULTS ARE CALCULATED.  ONE SET IS ACCU-      
!  RATE FOR LARGE VALUES OF K AND THE OTHER IS ACCURATE FOR K CLOSE     
!  TO ZERO.  THE TWO SETS ARE JOINED AT THE K VALUE FOR WHICH THE       
!  DIFFERENCE BETWEEN THEM IS LEAST.  FOR L=0 OR L=1 THE SMALL K        
!  RESULTS ARE CALCULATED USING SIEGMANS METHOD.                        
!                                                                       
!  IF THE LOGICAL VARIABLE SAVED IS .FALSE. AUXILIARY DATA IN THE       
!  ARRAYS TA, TRBE, AND WW ARE CALCULATED ANEW.  THIS MUST BE DONE      
!  IF THE MESH PARAMETERS N, H, RHOMIN OR KAPMIN HAVE CHANGED FROM      
!  THE PREVIOUS CALL TO THE SUBROUTINE AND ON THE FIRST CALL TO         
!  THE SUBROUTINE.  OTHERWISE A SUBSTANTIAL IMPROVEMENT IN SPEED        
!  CAN BE OBTAINED BY TAKING SAVED=.TRUE..                              
!                                                                       
!  THE VARIABLE DISC IS THE DIFFERENCE BETWEEN THE SMALL K AND          
!  LARGE K VALUES FOR THE RESULT AT THE MATCHING POINT OF MINIMUM       
!  DISCREPANCY.  IT PROVIDES A ROUGH GUIDE TO THE ABSOLUTE ACC-         
!  URACY OF THE RESULTS AT THE MAXIMUM VALUE OF THE TRANSFORM.          
!  IT IS, HOWEVER, A GREAT OVERESTIMATE OF THE ABSOLUTE ERROR AT        
!  LARGE K VALUES, AND SMALL K VALUES IN THE CASE OF NON-ZERO L.        
!                                                                       
!  XA, TA, TRBE AND WW ARE COMPLEX WORKING AREA ARRAYS DIMENSIONED      
!  AS XA(NDIM), TA(NDIM), TRBE(NDIM,2) AND WW(NDIM) WHERE NDIM          
!  MUST BE AT LEAST AS LARGE AS 2**N.  AS NOTED, IT MAY BE USEFUL       
!  TO PRESERVE THE CONTENTS OF TA, TRBE AND WW BUT THE CONTENTS OF      
!  XA ARE NOT OF VALUE.                                                 
!                                                                       
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 XA(NDIM),TA(NDIM),TRBE(NDIM,2),WW(NDIM),Y1,Y2,W1,W2,ZZM
      DIMENSION F(NDIM),G(NDIM)
      REAL*8 KAPMIN
      LOGICAL SAVED
      NT=2**N
      IF (NT.GT.NDIM) GO TO 131
      NH=NT/2
      NHP=NH+1
      PI=4.D0*DATAN(1.D0)
      AN=NT
      DT=2.D0*PI/(H*AN)
      CM=DSQRT(2.D0*PI)/AN
      CL=0.25D0*H/AN
      IF (.NOT.SAVED) GO TO 100
!  CALCULATE RESULTS ACCURATE AT SMALL K VALUES                         
    1 IF (L.LT.2) GO TO 11
      AA=DEXP((L+1.5D0)*RHOMIN)
      BB=DEXP((L+1.5D0)*H)
      DO I=1,NT
        XA(I)=F(I)*AA
        AA=AA*BB
      ENDDO
      CALL NLOGN(N,XA,WW,NDIM)
      DO I=1,NH
        XA(I)=TA(I)*XA(I)
      ENDDO
      DO 4 JJ=1,L
       AA=2*JJ-0.5D0
       T=0.D0
       ZZM=1.d0/DCMPLX(AA,T)
       DO I=1,NH
         XA(I)=ZZM*XA(I)
         T=T+DT
       ENDDO
    4 CONTINUE
      DO I=NHP,NT
        XA(I)=0.D0
      ENDDO
      CALL NLOGN(N,XA,WW,NDIM)
      AA=DEXP((L-1.5D0)*KAPMIN)*CM
      BB=DEXP((L-1.5D0)*H)
      DO I=1,NT
        G(I)=DREAL(AA*XA(I))
        write(14,*)G(I)
        AA=AA*BB
      ENDDO
      GO TO 21                                                          
!  FOR L=0 OR 1 CALCULATE RESULTS ACCURATE AT SMALL K VALUES USING      
!  SIEGMANS METHOD                                                      
   11 LP=L+1
      DO I=NH,NT
        XA(I)=0.D0
      ENDDO
      AA=DEXP(3.D0*RHOMIN)
      BB=DEXP(3.D0*H)
      DO I=1,NH
        XX=AA*F(2*I-1)
        AA=AA*BB
        YY=AA*F(2*I)
        XA(I)=DCMPLX(XX,YY)
        AA=AA*BB
      ENDDO
      CALL NLOGN(N,XA,WW,NDIM)
      Y1=TRBE(1,LP)*XA(1)
      Y2=TRBE(1,LP)*DCONJG(XA(1))
      XA(1)=2.D0*(Y1+Y2+DCONJG(Y2-Y1))
      XA(NHP)=4.D0*DCONJG(XA(NHP)*TRBE(NHP,LP))
      DO I=2,NH
      IC=NT-I+2
        Y1=XA(I)
        Y2=DCONJG(XA(IC))
        W1=Y1+Y2
        W2=WW(I+NH)*(Y1-Y2)
        Y1=(W1-W2)*TRBE(I,LP)
        Y2=(W1+W2)*DCONJG(TRBE(IC,LP))
        W1=Y1+Y2
        W2=WW(I+NH)*(Y1-Y2)
        XA(I)=W1+W2
        XA(IC)=DCONJG(W1-W2)
      ENDDO
      CALL NLOGN(N,XA,WW,NDIM)
      DO I=1,NH
        G(2*I-1)=CL*DREAL(XA(I))
        G(2*I)=CL*DIMAG(XA(I))
      ENDDO
                                                                       
!  CALCULATE RESULTS ACCURATE AT LARGE K VALUES                         
                                                                       
   21 AA=DEXP(1.5D0*RHOMIN)
      BB=DEXP(1.5D0*H)
      DO I=1,NT
        XA(I)=AA*F(I)
        AA=AA*BB
      ENDDO
      CALL NLOGN(N,XA,WW,NDIM)
      IJ=MOD(L,2)
      IJK=MOD(L,4)
      DO I=1,NH
        Y1=XA(I)*TA(I+IJ*NH)
        IF (IJK.GT.1) Y1=-Y1
        XA(I)=Y1
      ENDDO
      IF (L.EQ.0) GO TO 24
      DO 25 JJ=1,L
        AA=2*JJ-L-0.5D0
        BB=JJ-0.5D0
        T=0.D0
        DO I=1,NH
          XA(I)=XA(I)*DCMPLX(BB,-T)/DCMPLX(AA,T)
          T=T+DT
        ENDDO
   25 CONTINUE
   24 DO I=NHP,NT
        XA(I)=0.D0
      ENDDO
      CALL NLOGN(N,XA,WW,NDIM)
      AA=DEXP(-1.5D0*KAPMIN)*CM
      BB=DEXP(-1.5D0*H)
      DO I=1,NT
        XA(I)=AA*XA(I)
        AA=AA*BB
      ENDDO
!                                                                       
!  FIND POINT OF MINIMUM DISCREPANCY BETWEEN SMALL K AND LARGE K VALUES 
!  OF TRANSFORM AND CONSTRUCT BEST APPROXIMATE RESULT.                  
!                                                                       
      D1=DABS(G(1)-DREAL(XA(1)))
      D2=DABS(G(2)-DREAL(XA(2)))
      TE=D1+D2
      II=2
      D1=D2
      DO 31 I=3,NT
        D2=DABS(G(I)-DREAL(XA(I)))
        AA=D1+D2
        IF (AA.GE.TE) GO TO 32
        II=I
        TE=AA
   32   D1=D2
   31 CONTINUE
      DO I=II,NT
        G(I)=DREAL(XA(I))
      ENDDO
      DISC=TE
      RETURN
    !                                                                       
!  INITIALIZE AUXILIARY DATA FOR TRANSFORM                              
!                                                                       
!  THE NEXT INSTRUCTIONS INITIALIZE THE ARRAY TA USED IN CALCULATING    
!  THE TRANSFORM AT LARGE K VALUES AND SMALL K VALUES FOR L GREATER     
!  THAN 1.                                                              
!                                                                       
  100 CONTINUE
      Y1=1.D0
      AA=(RHOMIN+KAPMIN)*DT
      Y2=DCMPLX(DCOS(AA),DSIN(AA))
      DO I=1,NH
        T=(I-1)*DT
        S=0.D0
        RR=DSQRT(110.25D0+T*T)
        PHI=DATAN(T/10.5D0)
        DO NA=1,10
          S=S+DATAN(T/(NA-0.5D0))
        ENDDO
        S=S-T*DLOG(RR)+T-10.D0*PHI+DSIN(PHI)/(12.D0*RR)
        S=S-DSIN(3.D0*PHI)/(360.D0*RR**3)&
     &     +DSIN(5.D0*PHI)/(1260.D0*RR**5)&
     &     -DSIN(7.D0*PHI)/(1680.D0*RR**7)
        XX=DEXP(PI*T)
        XX=DATAN((XX-1.D0)/(XX+1.D0))
        CC=S-XX
        TA(I)=Y1*DCMPLX(DCOS(CC),DSIN(CC))
        CC=S+XX
        TA(I+NH)=Y1*DCMPLX(DCOS(CC),DSIN(CC))
        Y1=Y1*Y2
      ENDDO
      TA(1)=TA(1)/2.D0
      TA(1+NH)=TA(1+NH)/2.D0
      TA(NH)=TA(NH)/2.D0
      TA(NT)=TA(NT)/2.D0
!                                                                       
!  THE NEXT INSTRUCTIONS INITIALIZE THE ARRAY WW USED BY                
!  THE NLOGN SUBROUTINE.  THE ELEMENTS IN THE SECOND HALF OF THE        
!  ARRAY WW ARE USED IN THE IMPLEMENTATION OF SIEGMANS METHOD.          
!                                                                       
      DO I=1,NH
        XX=(I-1)*PI/AN
        WW(I+NH)=DCMPLX(-DSIN(XX),DCOS(XX))
        XX=2.D0*XX
        WW(I)=DCMPLX(DCOS(XX),DSIN(XX))
      ENDDO
!                                                                       
!  THE NEXT INSTRUCTIONS INITIALIZE THE ARRAY TRBE USED IN THE          
!  IMPLEMENTATION OF SIEGMANS METHOD.                                   
!                                                                       
      DO I=1,NT
        XA(I)=0.D0
      ENDDO
      CF=DEXP(H)
      XX=DEXP(RHOMIN+KAPMIN)
      DO 113 I=1,NT
        AA=DSIN(XX)/XX
        XX=CF*XX
        BB=DSIN(XX)/XX
        XA(I)=DCMPLX(AA,BB)
        IF (XX.GT.1.D8) GO TO 114
        XX=XX*CF
  113 CONTINUE
  114 CALL NLOGN(N,XA,WW,NDIM)
      TRBE(1,1)=XA(1)
      TRBE(NHP,1)=DCONJG(XA(NHP))
      DO I=2,NH
        IC=NT-I+2
        Y1=XA(I)
        Y2=DCONJG(XA(IC))
        W1=Y1+Y2
        W2=WW(NH+I)*(Y1-Y2)
        Y1=W1-W2
        Y2=DCONJG(W1+W2)
        TRBE(I,1) =0.5D0*DCONJG(Y1)
        TRBE(IC,1)=0.5D0*DCONJG(Y2)
      ENDDO
      DO I=1,NT
        XA(I)=0.D0
      ENDDO
      XX=DEXP(RHOMIN+KAPMIN)
      DO 121 I=1,NT
       IF (XX.LT.0.1D0) GO TO 122
       AA=(DSIN(XX)/XX-DCOS(XX))/XX
       XX=CF*XX
       BB=(DSIN(XX)/XX-DCOS(XX))/XX
       XA(I)=DCMPLX(AA,BB)
       IF (XX.GT.1.D8) GO TO 123
       GO TO 121
  122  CC=XX*XX/2.D0
       CD=1.D0-CC/5.D0+CC*CC/70.D0-CC*CC*CC/1890.D0+CC**4/83160.D0
       AA=XX*CD/3.D0
       XX=CF*XX
       CC=XX*XX/2.D0
       CD=1.D0-CC/5.D0+CC*CC/70.D0-CC*CC*CC/1890.D0+CC**4/83160.D0
       BB=XX*CD/3.D0
       XA(I)=DCMPLX(AA,BB)
  121  XX=XX*CF
  123 CALL NLOGN(N,XA,WW,NDIM)
      TRBE(1,2)=XA(1)
      TRBE(NHP,2)=DCONJG(XA(NHP))
      DO I=2,NH
        IC=NT-I+2
        Y1=XA(I)
        Y2=DCONJG(XA(IC))
        W1=Y1+Y2
        W2=WW(NH+I)*(Y1-Y2)
        Y1=W1-W2
        Y2=DCONJG(W1+W2)
        TRBE(I,2) =0.5D0*DCONJG(Y1)
        TRBE(IC,2)=0.5D0*DCONJG(Y2)
      ENDDO
      GO TO 1
  131 WRITE (6,132) NDIM,NT
  132 FORMAT(1x,'SUBROUTINE LSFBTR NOT EXECUTED. DIMENSION',&
     &    ' NDIM TOO SMALL: NDIM =',i10,'  NT = 2**N =',i10)
      RETURN
      END
!     ==================================================================
      SUBROUTINE NLOGN(N,X,WW,NDIM)
!     ==--------------------------------------------------------------==
!   FAST FOURIER TRANSFORM ROUTINE                                      
!                                                                       
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(NDIM),MM(15)
      COMPLEX*16 WW(NDIM),X,WK,HOLD,Q
      DO I=1,N
        MM(I)=2**(N-I)
      ENDDO
      LX=2*MM(1)
      DO 4 L=1,N
      NBLOCK=2**(L-1)
      LBLOCK=LX/NBLOCK
      LBHALF=LBLOCK/2
      K=0
      DO 4 IBLOCK=1,NBLOCK
      WK=WW(K+1)
      ISTART=LBLOCK*(IBLOCK-1)
      DO 2 I=1,LBHALF
        J=ISTART+I
        JH=J+LBHALF
        Q=X(JH)*WK
        X(JH)=X(J)-Q
        X(J)=X(J)+Q
   2  CONTINUE
      DO 3 I=2,N
        II=I
        IF (K.LT.MM(I)) GO TO 4
        K=K-MM(I)
   3  CONTINUE   
   4  K=K+MM(II)
      K=0
      DO 7 J=1,LX
        IF (K.LT.J) GO TO 5
        HOLD= X(J)
        X(J)=X(K+1)
        X(K+1)= HOLD
   5    DO 6 I=1,N
          II=I
          IF (K.LT.MM(I)) GO TO 7
          K=K-MM(I)
   6    CONTINUE
   7  K=K+MM(II)
      RETURN
      END
