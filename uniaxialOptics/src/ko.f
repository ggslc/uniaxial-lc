c     ko.f
c     copyright (C) 2007  Steph Cornford
c     copyright (C) 1988  D Y K Ko, J R Sambles

c     This program is free software; you can redistribute it and/or
c     modify it under the terms of the GNU General Public License
c     as published by the Free Software Foundation; either version 2
c     of the License, or (at your option) any later version.

c     This program is distributed in the hope that it will be useful,
c     but WITHOUT ANY WARRANTY; without even the implied warranty of
c     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c     GNU General Public License for more details

c     You should have received a copy of the GNU General Public License
c     along with this program; if not, write to the Free Software
c     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


      SUBROUTINE AINTER(NLAY,NANG,ZE1,ZE3,TWIST,TILT
     .     ,DEPTH,LAMBDA,ICOEFF,ECOEFF,ANGLE,AZIMA)

C     R interface, angle and azimuth are pairs
C     doesn't really do much, just reduces the
C     overhead of repeated calls/copies 

      IMPLICIT NONE
c     scalar args
      REAL*8 LAMBDA
      INTEGER NLAY,NANG

c     vector args
      COMPLEX*16 ZE1(NLAY),ZE3(NLAY),ECOEFF(8,NANG)
      REAL*8    TWIST(NLAY),TILT(NLAY),DEPTH(NLAY),ANGLE(NANG),
     .     ICOEFF(8,NANG), AZIMA(NANG)

c     local scalars
      INTEGER I,J

c     local vectors
      REAL*8 IC(8)
      COMPLEX*16 EC(8)

      DO I = 1,NANG
         CALL CINTER(NLAY,1 ,ZE1,ZE3,TWIST,TILT
     .        ,DEPTH,LAMBDA,ICOEFF(1,I),ECOEFF(1,I),
     .        ANGLE(I),AZIMA(I))

c$$$         DO J = 1,8
c$$$            ICOEFF(J,I) = IC(J)
c$$$            ECOEFF(J,I) = EC(J)
c$$$         END DO

      END DO
      
      RETURN
      END


      SUBROUTINE CINTER(NLAY,NANG,ZE1,ZE3,TWIST,TILT
     .     ,DEPTH,LAMBDA,ICOEFF,ECOEFF,ANGLE,AZIMA)

C     R interface, run over all angles and one azima 

      IMPLICIT NONE
c     scalar args
      REAL*8 LAMBDA,AZIMA
      INTEGER NLAY,NANG

c     vector args
      COMPLEX*16 ZE1(NLAY),ZE3(NLAY),ECOEFF(8,NANG)
      REAL*8    TWIST(NLAY),TILT(NLAY),DEPTH(NLAY),ANGLE(NANG),
     .     ICOEFF(8,NANG)

c     local scalars
      REAL*8 RKVECT, ZF, PI, EOMEGA, C, MU0, E0, FACTOR
      INTEGER I, ERRFLAG

c     local vectors 
c      COMPLEX*16 ZE1(256),ZE3(256)


      IF (NLAY.GT.4096) THEN
         write(*,*) "rebuild fortran with more layers"
         RETURN
      END IF

c      DO I = 1,NLAY
c         ZE1(I) = DCMPLX(EPS1R(I),EPS1I(I))
c         ZE3(I) = DCMPLX(EPS3R(I),EPS3I(I))
c      END DO
       PI = 4.D0 * DATAN(1.0D0)
      C  = 2.99792458E+08
      MU0 = 4.0D-7*PI
      E0 = 1.0D0/(C*C*MU0)
      FACTOR = 1.0D0

     
      RKVECT  = 2.0*PI/LAMBDA
      ZF = RKVECT*C*MU0
      EOMEGA = E0*RKVECT*C

      ERRFLAG = 0

      CALL MODEL(NLAY,ZE1,ZE3,TWIST,TILT,DEPTH,1,NANG,AZIMA,FACTOR,
     .     RKVECT,ZF,EOMEGA,ICOEFF,ECOEFF,ANGLE,ERRFLAG)

      RETURN
      END



*     The model calculates optical intensities for system at input angles
      SUBROUTINE MODEL(NM,ZE1,ZE3,TWIN,TILT,D,NAZIM,NUMX,AZIMA,FACTOR,
     $           RKVECT,ZF,EOMEGA,ICOEFF,ECOEFF,ANGLE,ERRFLAG)

      IMPLICIT NONE

C     Define the variabls transferred in and out

      INTEGER     NM,NAZIM,NUMX,NMXB1,ERRFLAG
      REAL*8      TWIN,TILT,D,AZIMA,FACTOR,RKVECT,ZF,EOMEGA,ICOEFF
      REAL*8      ANGLE
C      REAL*8      NAN
      COMPLEX*16  ZE1,ZE3,ZZRR,ZZTT,ZQ,ZPHI,ECOEFF

      DIMENSION   TWIN(NM),TILT(NM),D(NM),ICOEFF(8,NUMX),ECOEFF(8,NUMX)
      DIMENSION   ZE1(NM),ZE3(NM),ANGLE(NUMX)
      DIMENSION   ZZRR(2,2),ZZTT(2,2),ZQ(4,4096)
      DIMENSION   ZPHI(4,4,4096)

C     Define the variables used internally

      INTEGER     AZ,I,J,N,M,A,B
      REAL*8      TWIST,TWM,TIM,STWI,CTWI,STIL,CTIL
      REAL*8      STWI2,CTWI2,STIL2,REF1,REFNM,THETA
      REAL*8      XCOMP,XCOMP2,RKX,FACTR,FACT1
      REAL*8      FACTIV,FACTT,TMP
      COMPLEX*16  ZEDIFF,ZDS1,ZDS2,ZQP1,ZQP2,ZAA,ZA,ZB,ZC
      COMPLEX*16  ZD,ZZE,ZG,ZQ3,ZQ4,ZSUMV,ZE33,INVZE33,ZMP

      DIMENSION   TWIST(4096),TWM(4096),TIM(4096),STWI(4096)
      DIMENSION   CTWI(4096),STIL(4096),CTIL(4096),STWI2(4096)
      DIMENSION   CTWI2(4096)
      DIMENSION   STIL2(4096),FACTR(2,2),FACTT(2,2),TMP(8),ZMP(8)
      DIMENSION   ZEDIFF(4096),ZE33(4096),INVZE33(4096)


      REAL* 8 TOLTIL,TOLTWI

C     Loops over the layers setting the correct twist values

          DO 162 M=1,NM

            TWIST(M)=TWIN(M)+AZIMA

C     ZE1 and ZE3 are the dielectric constants of the layers (ZE2 is
C       assumed to be equal to ZE1 in the calculations
C     TWIST and TILT are the Euler angles used to define the
C       coordinate system in which the dielectric tensor is diagonalised
C     D is the layer thickness

c    c        TOL
c$$$
c$$$	    IF (DABS(DABS(TWIST(M))-90.0).LT.1.0D-3) THEN
c$$$              IF (TWIST(M).GE.90.0) THEN
c$$$	        TWIST(M)=TWIST(M)+0.001
c$$$	      ELSE IF (TWIST(M).LT.0.0.AND.TWIST(M).GE.-90.0) THEN
c$$$	        TWIST(M)=TWIST(M)+0.001
c$$$	      ELSE
c$$$	        TWIST(M)=TWIST(M)-0.001
c$$$	      END IF
c$$$	    END IF
c$$$
c$$$	    IF (DABS(TWIST(M)).LT.1.0D-3) THEN
c$$$	      IF (TWIST(M).GE.0.00) THEN
c$$$	        TWIST(M)=TWIST(M)+0.001
c$$$	      ELSE
c$$$	        TWIST(M)=TWIST(M)-0.001
c$$$	      END IF
c$$$	    END IF
c$$$	
c$$$	    IF (DABS(DABS(TWIST(M))-180.0).LT.1.0D-3) THEN
c$$$	      IF (TWIST(M).GE.180.0) THEN
c$$$	        TWIST(M)=TWIST(M)+0.001
c$$$	      ELSE IF (TWIST(M).LT.0.0.AND.TWIST(M).GE.-180.0) THEN
c$$$	        TWIST(M)=TWIST(M)+0.001
c$$$	      ELSE
c$$$	        TWIST(M)=TWIST(M)-0.001
c$$$	      END IF
c$$$	    END IF
c$$$	
c$$$	    IF (DABS(DABS(TILT(M))-90.0).LT.1.0D-3) THEN
c$$$	      IF (TILT(M).GE.90.0) THEN
c$$$	        TILT(M)=TILT(M)+0.001
c$$$	      ELSE IF (TILT(M).LT.0.0.AND.TILT(M).GE.-90.0) THEN
c$$$	        TILT(M)=TILT(M)+0.001
c$$$	      ELSE
c$$$	        TILT(M)=TILT(M)-0.001
c$$$	      END IF
c$$$	    END IF
c$$$
c$$$	    IF (DABS(TILT(M)).LT.1.0D-3) THEN
c$$$	      IF (TILT(M).GE.0.00) THEN
c$$$	        TILT(M)=TILT(M)+0.001
c$$$	      ELSE
c$$$	        TILT(M)=TILT(M)-0.001
c$$$	      END IF
c$$$	    END IF
c$$$	
c$$$	    IF (DABS(DABS(TILT(M))-180.0).LT.1.0D-3) THEN
c$$$	      IF (TILT(M).GE.180.0) THEN
c$$$	        TILT(M)=TILT(M)+0.001
c$$$	      ELSE IF (TILT(M).LT.0.0.AND.TILT(M).GE.-180.0) THEN
c$$$	        TILT(M)=TILT(M)+0.001
c$$$	      ELSE
c$$$	        TILT(M)=TILT(M)-0.001
c$$$	      END IF
c$$$	    END IF



	
	    TWM(M)=TWIST(M)*FACTOR
	    TIM(M)=TILT(M)*FACTOR
	    STWI(M)=DSIN(TWM(M)) 
	    CTWI(M)=DCOS(TWM(M)) 
	    STIL(M)=DSIN(TIM(M)) 
	    CTIL(M)=DCOS(TIM(M)) 
 	    STWI2(M)=STWI(M)*STWI(M)
	    CTWI2(M)=CTWI(M)*CTWI(M)
	    STIL2(M)=STIL(M)*STIL(M)
	    ZEDIFF(M)=ZE3(M)-ZE1(M)
	    ZE33(M)=ZE3(M)-ZEDIFF(M)*STIL2(M)
            INVZE33(M)=1/ZE33(M)

 162      CONTINUE 

C
C     Refractive indices of both isotropic semi-infinite media
C
          REF1=DSQRT(DREAL(ZE1(1)))
          REFNM=DSQRT(DREAL(ZE1(NM)))

C     Angle loop begins, calculating only at the angle data points

          DO 163 N=1,NUMX

              THETA=ANGLE(N)*FACTOR
              

C     Calculate scattering matrix

	    XCOMP=(REF1*DSIN(THETA))
	    XCOMP2=XCOMP*XCOMP
	    RKX=RKVECT*XCOMP

	    DO 164 M=1,NM

	      ZQ(1,M)=RKVECT*CDSQRT(ZE1(M)-XCOMP2)
	      ZQ(2,M)=-ZQ(1,M)
	      ZDS1=ZEDIFF(M)*STIL(M)
	      ZDS2=ZDS1*STIL(M)
	      ZQP1=-XCOMP*ZDS1*CTIL(M)*STWI(M)
	      ZQP2=ZE33(M)
     $        -((1.0D0-(ZDS2/ZE3(M))*CTWI2(M))*XCOMP2)
	      ZQP2=CDSQRT(ZQP2*ZE1(M)*ZE3(M))
	      ZQ(3,M)=RKVECT*(ZQP1+ZQP2)*INVZE33(M)
	      ZQ(4,M)=RKVECT*(ZQP1-ZQP2)*INVZE33(M)

	      ZAA=-RKX*ZDS1*CTIL(M)*INVZE33(M)
	      ZA=ZAA*STWI(M)
	      ZB=ZF*(ZE33(M)-XCOMP2)*INVZE33(M)
	      ZC=-ZAA*CTWI(M)
	      ZD=ZE3(M)-ZDS2*CTWI2(M)
	      ZD=ZD*EOMEGA*ZE1(M)*INVZE33(M)
	      ZZE=ZE1(M)*ZDS2*CTWI(M)*STWI(M)
	      ZZE=-ZZE*EOMEGA*INVZE33(M)
	      ZG=ZE1(M)
     $        *(ZE3(M)-(ZDS2*STWI2(M)))-ZE33(M)*XCOMP2
              ZG=ZG*EOMEGA*INVZE33(M)

	      IF (CDABS(ZDS2*CTWI2(M)).LT.1.0D-16) THEN
                ZQ3=ZQ(3,M)
	        ZQ4=ZQ(4,M)
	        ZQ(3,M)=ZQ(1,M)
	        ZQ(4,M)=ZQ(2,M)
	        ZQ(1,M)=ZQ3
	        ZQ(2,M)=ZQ4
              END IF

	      IF (CDABS(ZDS2*CTWI2(M)).LT.1.0D-16.OR.CDABS(ZDS2*STWI2(M))
     $        .LT.1.0D-16.AND.CDABS(XCOMP*(ZE33(M)-ZE1(M)))
     $        .LE.1.0D-16) THEN

                ZQ3=ZQ(3,M)
	        ZQ4=ZQ(4,M)
	        ZQ(3,M)=ZQ(1,M)
	        ZQ(4,M)=ZQ(2,M)
	        ZQ(1,M)=ZQ3
	        ZQ(2,M)=ZQ4

                DO 167 I=1,2

	          ZPHI(1,I,M)=0.0D0
	          ZPHI(2,I,M)=0.0D0
	          ZPHI(3,I,M)=1.0D0
	          ZPHI(4,I,M)=ZQ(I,M)/ZF
            
	          ZSUMV=CDSQRT(1.0D0+ZPHI(4,I,M)*CONJG(ZPHI(4,I,M)))
 
	          DO 168 J=3,4
		    ZPHI(J,I,M)=ZPHI(J,I,M)/ZSUMV
 168              CONTINUE
 167            CONTINUE


	        DO 169 I=3,4

  	          ZPHI(1,I,M)=1.0D0
	          ZPHI(2,I,M)=(ZQ(I,M)-ZA)/ZB
	          ZPHI(3,I,M)=0.0D0
	          ZPHI(4,I,M)=0.0D0

	          ZSUMV=CDSQRT(1.0D0+ZPHI(2,I,M)*CONJG(ZPHI(2,I,M)))

	          DO 170 J=1,2
		    ZPHI(J,I,M)=ZPHI(J,I,M)/ZSUMV
 170              CONTINUE
 169            CONTINUE
              ELSE
                DO 172 I=1,4

	          ZPHI(1,I,M)=ZB*ZZE-(ZA-ZQ(I,M))*ZC
	          ZPHI(2,I,M)=(ZD*ZC-(ZA-ZQ(I,M))*ZZE)
	          ZPHI(3,I,M)=((ZA-ZQ(I,M))**2-ZB*ZD)
	          ZPHI(4,I,M)=ZPHI(3,I,M)*ZQ(I,M)/ZF

  	          ZSUMV=ZPHI(1,I,M)*CONJG(ZPHI(1,I,M))
     $            +ZPHI(2,I,M)*CONJG(ZPHI(2,I,M))
     $            +ZPHI(3,I,M)*CONJG(ZPHI(3,I,M))
     $            +ZPHI(4,I,M)*CONJG(ZPHI(4,I,M))
	          ZSUMV=CDSQRT(ZSUMV)

	          DO 174 J=1,4
		    ZPHI(J,I,M)=ZPHI(J,I,M)/ZSUMV
 174              CONTINUE
 172            CONTINUE

              END IF
         
 164        CONTINUE

          CALL SMAT(D,NM,ZZRR,ZZTT,ZQ,ZPHI,ERRFLAG)
          IF (ERRFLAG.NE.0) THEN
C            NAN=LOG(-1.0)
C            DO 500 A=1,8
C              DO 501 B=1,NUMX
C              OUT(A,B)=NAN  
C 501          CONTINUE
C 500        CONTINUE
            WRITE(*,*) 'Angle = ',ANGLE(N)
            RETURN
          END IF

C     COMPLEX REFLECTION COEFFICIENTS:
C       (p)in(p)out=ZZRR(1,1)
C       (p)in(s)out=ZZRR(2,1)
C       (s)in(p)out=ZZRR(1,2)
C       (s)in(s)out=ZZRR(2,2)
C     COMPLEX TRANSMISSION COEFFICIENTS:
C       (p)in(p)out=ZZTT(1,1)
C       (p)in(s)out=ZZTT(2,1)
C       (s)in(p)out=ZZTT(1,2)
C       (s)in(s)out=ZZTT(2,2)
C
C     Multiply all transmission coefficients by a factor
C       to compensate for change in beam area
C       between first and last media.
C     Another factor is also required to
C       compensate for the fact that ratios are
C       made to Ex which is always 1.
C     Hence this factor is required for Rsp, Rps, Tps, Tsp, Tss.

	  FACTR(1,1)=1.0D0
	  FACTR(2,2)=1.0D0
	  FACTR(1,2)=1.0D0/(DCOS(THETA))**2
	  FACTR(2,1)=1.0D0/FACTR(1,2)
 
	  FACT1=REF1/REFNM*DSIN(THETA)
	  FACT1=FACT1*FACT1
	  FACT1=1.0D0-FACT1
	  FACTIV=REFNM/REF1

	  IF (FACT1.LT.0.0D0) THEN
	    FACTT(2,2)=0.0D0
	    FACTT(1,1)=0.0D0
	    FACTT(1,2)=0.0D0
	    FACTT(2,1)=0.0D0
	  ELSE
	    FACT1=DSQRT(FACT1)
	    FACTT(2,2)=(FACTIV*FACT1)/DCOS(THETA)
	    FACTT(1,1)=FACTT(2,2)*(DCOS(THETA))**2/(FACT1**2)
	    FACTT(1,2)=FACTT(2,2)/(FACT1**2)
	    FACTT(2,1)=FACTT(2,2)*(DCOS(THETA)**2)
	  END IF

	  DO 176 I=1,2
	    DO 176 J=1,2
	      IF (DREAL(CONJG(ZZRR(I,J))*ZZRR(I,J)).LE.
     $        1.0E-16)ZZRR(I,J)=CMPLX(0.0,0.0)
	      IF (DREAL(CONJG(ZZTT(I,J))*ZZTT(I,J)).LE.
     $        1.0E-16)ZZTT(I,J)=CMPLX(0.0,0.0)
 176      CONTINUE
 
C     Ta-da! The answers! Elements one to eight are the conventional
C       Rpp, Rps,... Tss as always
       
          TMP(1)=DREAL(CONJG(ZZRR(1,1))*ZZRR(1,1))*FACTR(1,1) 
          TMP(2)=DREAL(CONJG(ZZRR(2,1))*ZZRR(2,1))*FACTR(2,1) 
          TMP(3)=DREAL(CONJG(ZZRR(1,2))*ZZRR(1,2))*FACTR(1,2) 
          TMP(4)=DREAL(CONJG(ZZRR(2,2))*ZZRR(2,2))*FACTR(2,2) 
          TMP(5)=DREAL(CONJG(ZZTT(1,1))*ZZTT(1,1))*FACTT(1,1) 
          TMP(6)=DREAL(CONJG(ZZTT(2,1))*ZZTT(2,1))*FACTT(2,1) 
          TMP(7)=DREAL(CONJG(ZZTT(1,2))*ZZTT(1,2))*FACTT(1,2) 
          TMP(8)=DREAL(CONJG(ZZTT(2,2))*ZZTT(2,2))*FACTT(2,2) 


          ZMP(1)=ZZRR(1,1)*DSQRT(FACTR(1,1)) 
          ZMP(2)=ZZRR(2,1)*DSQRT(FACTR(2,1)) 
          ZMP(3)=ZZRR(1,2)*DSQRT(FACTR(1,2)) 
          ZMP(4)=ZZRR(2,2)*DSQRT(FACTR(2,2)) 
          ZMP(5)=ZZTT(1,1)*DSQRT(FACTT(1,1)) 
          ZMP(6)=ZZTT(2,1)*DSQRT(FACTT(2,1)) 
          ZMP(7)=ZZTT(1,2)*DSQRT(FACTT(1,2)) 
          ZMP(8)=ZZTT(2,2)*DSQRT(FACTT(2,2)) 

C     Pop them into the relevant output array
       
          DO 177 I=1,8
            ICOEFF(I,N)=TMP(I)  
            ECOEFF(I,N)=ZMP(I)
 177      CONTINUE

C     End of the angle loop
      
 163    CONTINUE

      RETURN
      END

*     MODEL ends

*     ------------------

c$$$*     INVMAL
c$$$*     inverse of a 4*4 matrix using LAPACK
c$$$
c$$$      SUBROUTINE INVMAL(ZA,ZB,ZC,ERRFLAG)
c$$$      IMPLICIT NONE
c$$$
c$$$*     arguments
c$$$      COMPLEX*16 ZA(4,4),ZB(4,4),ZC(4,4)
c$$$      INTEGER ERRFLAG
c$$$
c$$$*     locals
c$$$      INTEGER IPIV(4,4), I, J
c$$$      COMPLEX*16 CWK(4096)
c$$$
c$$$*     compute the LU decomposition
c$$$      CALL ZGETRF( 4, 4, ZA, 4, IPIV, ERRFLAG )
c$$$*     compute the inverse
c$$$      CALL ZGETRI( 4, ZA, 4, IPIV, CWK, 4096, ERRFLAG )
c$$$
c$$$      DO I = 1,4
c$$$         DO J = 1,4
c$$$*            ZB(I,J) = 0.0D0
c$$$            ZC(I,J) = ZA(I,J)
c$$$         END DO
c$$$*         ZB(I,I) = 1.0D0
c$$$      END DO
c$$$
c$$$
c$$$      RETURN
c$$$      END


*     INVMAT

*     Calculates the inverse of a 4*4 matrix

      SUBROUTINE INVMAT(ZA,ZB,ZC,ERRFLAG)

      IMPLICIT NONE

      INTEGER    I,J,INDX,ERRFLAG
      COMPLEX*16 ZB,ZC,ZA

      DIMENSION  ZA(4,4),ZB(4,4),ZC(4,4),INDX(4)

      DO 211 I=1,4

	DO 212 J=1,4
	  ZB(I,J)=0.0D0
	  ZC(I,J)=0.0D0
 212     CONTINUE
	ZB(I,I)=1.0D0
	ZC(I,I)=1.0D0
 211   CONTINUE

      CALL LUDCMP(ZA,4,4,INDX,ERRFLAG)
      IF (ERRFLAG.NE.0) RETURN

      DO 213 J=1,4
	CALL LUBKSB(ZA,4,4,INDX,ZC(1,J))
 213  CONTINUE

      RETURN
      END

*     INVMAT ends

*     ------------------

*     LUBKSB

      SUBROUTINE LUBKSB(A,N,NP,INDX,B)

      IMPLICIT NONE

      INTEGER    I,N,II,LL,J,NP,INDX
      COMPLEX*16 A,B,SUM

      DIMENSION INDX(N),A(NP,NP),B(N)

      II=0
      DO 221 I=1,N
	LL=INDX(I)
	SUM=B(LL)
	B(LL)=B(I)
	IF (II.NE.0)THEN
	  DO 222 J=II,I-1
	    SUM=SUM-A(I,J)*B(J)
 222      CONTINUE
	ELSE IF (CDABS(SUM).GE.1.0D-12) THEN
	  II=I
	END IF
	B(I)=SUM
 221  CONTINUE
      DO 223 I=N,1,-1
	SUM=B(I)
	IF (I.LT.N)THEN
	  DO 224 J=I+1,N
	    SUM=SUM-A(I,J)*B(J)
 224      CONTINUE
	END IF
	B(I)=SUM/A(I,I)
 223  CONTINUE

      RETURN
      END

*     LUBKSB ends

*     ------------------

*     LUDCMP

*     Performs an L-U Decomposition

      SUBROUTINE LUDCMP(A,N,NP,INDX,ERRFLAG)

      IMPLICIT NONE

      INTEGER    I,N,J,K,IMAX,INDX,NP,ERRFLAG
      REAL*8     D,TINY
      COMPLEX*16 A,VV,AAMAX,SUM,DUM

      DIMENSION  INDX(N),A(NP,NP),VV(4096)

      PARAMETER  (TINY=1.0D-20)

      DO 231 I=1,N
	AAMAX=(0.0D00,0.0D00)
	DO 232 J=1,N
	  IF (CDABS(A(I,J)).GT.CDABS(AAMAX)) AAMAX=A(I,J)
 232    CONTINUE
	IF (CDABS(AAMAX).LT.TINY) THEN
          WRITE(*,*) 'Singular matrix.'
          ERRFLAG=1
          RETURN
        END IF
	VV(I)=1.0D0/AAMAX
 231  CONTINUE

      DO 233 J=1,N

	IF (J.GT.1) THEN
	  DO 234 I=1,J-1
	    SUM=A(I,J)
	    IF (I.GT.1)THEN
	      DO 235 K=1,I-1
		SUM=SUM-A(I,K)*A(K,J)
 235          CONTINUE
	      A(I,J)=SUM
	    END IF
 234      CONTINUE
	ENDIF

	AAMAX=(0.0D00,0.0D00)

	DO 236 I=J,N
	  SUM=A(I,J)
	  IF (J.GT.1)THEN
	    DO 237 K=1,J-1
	      SUM=SUM-A(I,K)*A(K,J)
 237        CONTINUE
	    A(I,J)=SUM
	  END IF
	  DUM=VV(I)*SUM
	  IF (CDABS(DUM).GE.CDABS(AAMAX)) THEN
	    IMAX=I
	    AAMAX=DUM
	  END IF
 236    CONTINUE

	IF (J.NE.IMAX)THEN
	  DO 238 K=1,N
	    DUM=A(IMAX,K)
	    A(IMAX,K)=A(J,K)
	    A(J,K)=DUM
 238      CONTINUE
	  VV(IMAX)=VV(J)
	END IF

	INDX(J)=IMAX
	IF (J.NE.N)THEN
	  IF (CDABS(A(J,J)).LT.TINY) A(J,J)=TINY
	  DUM=1.0D00/A(J,J)
	  DO 239 I=J+1,N
	    A(I,J)=A(I,J)*DUM
 239      CONTINUE
	END IF
 233  CONTINUE

      IF (CDABS(A(N,N)).LT.TINY) A(N,N)=TINY

      RETURN
      END

*     LUDCMP ends

*     ------------------

*     CMATXL
*     Multiplication of matrices using BLAS
      SUBROUTINE CMATXL(A,B,C,N)

*     returns C = A * B

      IMPLICIT NONE

*     arguments
      INTEGER N
      COMPLEX*16 A(N,N), B(N,N), C(N,N)

      CALL ZGEMM ( 'n', 'n', N, N, N, 
     .     DCMPLX(1.0D0,0.0D0), A, N, B, N,
     .     DCMPLX(0.0D0,0.0D0), C, N )



      RETURN
      END


*     CMATX1

*     Multiplication of matrices (ND1=4 only)

      SUBROUTINE CMATX1(ZS,ZQ,ZXY,ND1)

      IMPLICIT NONE

      INTEGER    I,J,K,ND1
      COMPLEX*16 ZS,ZQ,ZXY

      DIMENSION  ZS(ND1,ND1),ZQ(ND1,ND1),ZXY(ND1,ND1)

      DO 242 I=1,ND1
	  ZXY(I,1)=ZS(I,1)*ZQ(1,1)+ZS(I,2)*ZQ(2,1)+ZS(I,3)*ZQ(3,1)+
     $    ZS(I,4)*ZQ(4,1)
	  ZXY(I,2)=ZS(I,1)*ZQ(1,2)+ZS(I,2)*ZQ(2,2)+ZS(I,3)*ZQ(3,2)+
     $    ZS(I,4)*ZQ(4,2)
	  ZXY(I,3)=ZS(I,1)*ZQ(1,3)+ZS(I,2)*ZQ(2,3)+ZS(I,3)*ZQ(3,3)+
     $    ZS(I,4)*ZQ(4,3)
	  ZXY(I,4)=ZS(I,1)*ZQ(1,4)+ZS(I,2)*ZQ(2,4)+ZS(I,3)*ZQ(3,4)+
     $    ZS(I,4)*ZQ(4,4)
 242    CONTINUE

      DO 243 J=1,ND1
        DO 244 K=1,ND1

	  ZS(J,K)=ZXY(J,K)

 244    CONTINUE
 243  CONTINUE

      RETURN
      END

*     CMATX1 ends

*     ------------------

*     SORT

*     Heapsort routine.
*     Indexes an array A of length N i.e. outputs the array INDX such that
*       A(INDX(J)) is in descending order for J=1,2,3,...,N

      SUBROUTINE SORT(A,N,INDX)

      IMPLICIT NONE

      INTEGER   INDX,I,J,N,L,IR,INDXT
      REAL*8    A,Q

      DIMENSION INDX(N),A(N)

      DO 251 I=1,N
        INDX(I)=I
 251  CONTINUE
      L=N/2+1
      IR=N

 252    IF (L.GT.1) THEN
	L=L-1
	INDXT=INDX(L)
	Q=A(INDXT)
      ELSE
	INDXT=INDX(IR)
	Q=A(INDXT)
	INDX(IR)=INDX(1)
	IR=IR-1
	IF (IR.EQ.1) THEN
	  INDX(1)=INDXT
	  RETURN
        ENDIF
      ENDIF

      I=L
      J=L+L

 253  IF (J.LE.IR) THEN
	IF (J.LT.IR.AND.A(INDX(J)).GT.A(INDX(J+1))) J=J+1
	IF (Q.GT.A(INDX(J))) THEN
	  INDX(I)=INDX(J)
	  I=J
	  J=J+J
	ELSE
	  J=IR+1
	ENDIF
	GOTO 253
      ENDIF

      INDX(I)=INDXT
      GOTO 252

      END

*     SORT ends

*     ------------------

*     SORT2

*     Sorts out the solutions into the required order

      SUBROUTINE SORT2(A,B,N,INDX)

      IMPLICIT NONE

      INTEGER   INDX,N,I
      REAL*8    A,B,EPS

      DIMENSION A(N),B(N),INDX(N)
 
      DATA EPS/1.D-30/

      IF (DABS(B(INDX(1))).LT.DABS(B(INDX(2)))) THEN
	I=INDX(1)
	INDX(1)=INDX(2)
	INDX(2)=I
      ENDIF

      IF (DABS(B(INDX(3))).LT.DABS(B(INDX(4)))) THEN
	I=INDX(3)
	INDX(3)=INDX(4)
	INDX(4)=I
      ENDIF
 
      IF (DABS(B(INDX(2))).LT.EPS) THEN
	IF (DABS(B(INDX(1))).LT.EPS) THEN
	  CALL SORT(A,4,INDX)
	  I=INDX(3)
	  INDX(3)=INDX(4)
	  INDX(4)=I
	  RETURN
	ENDIF
      ELSE
	IF (A(INDX(2)).LT.A(INDX(4))) THEN
	  I=INDX(2)
	  INDX(2)=INDX(4)
	  INDX(4)=I
	ENDIF
      ENDIF

      RETURN
      END

*     SORT2 ends

*     ------------------

*     SMAT

*     Calculates scattering matrix

      SUBROUTINE SMAT(D,NM,RMAT,TMAT,ZQ,ZPHI,ERRFLAG)

      IMPLICIT NONE

      INTEGER    INDX,NM,I,J,ILAY,IUSE,IFILL,ERRFLAG
      REAL*8     D,AR,AI,UR,UI
      COMPLEX*16 S,M,A,W,SNEW,EIG,RMAT,TMAT,DET,M1,M2,ZQ,ZPHI

      DIMENSION  S(4,4),M(4,4,2),A(4,4),W(4,4)
      DIMENSION  SNEW(4,4),EIG(4),RMAT(2,2),TMAT(2,2),M1(4,4),M2(4,4)
      DIMENSION  D(NM),AR(4),AI(4),UR(4,4),UI(4,4)
      DIMENSION  INDX(4),ZQ(4,NM),ZPHI(4,4,NM)

      DO 321 I=1,4
        DO 320 J=1,4
          S(I,J)=0.D0
 320    CONTINUE
 321  CONTINUE

      DO 40 I=1,4
        S(I,I)=1.D0
 40   CONTINUE

C     Calculate eigenvalues and eigenvectors in first medium

      DO 60 I=1,4
	AR(I)=DREAL(ZQ(I,1))
	AI(I)=DIMAG(ZQ(I,1))
	DO 55 J=1,4
	  UR(J,I)=DREAL(ZPHI(J,I,1))
	  UI(J,I)=DIMAG(ZPHI(J,I,1))
 55     CONTINUE
 60   CONTINUE


C     Sort out solutions into p and s polarisations

      INDX(1)=3
      INDX(2)=1
      INDX(3)=4
      INDX(4)=2
       
C     This bit gives the correct configuration of the forward and
C       backward propagating p and s eigenmodes in the first medium :    
C          forward p  = INDX(1);       forward s  = INDX(2);
C          backward p = INDX(3);       backward s = INDX(4)

      EIG(1)=CDEXP(DCMPLX(-AI(INDX(1)),AR(INDX(1)))*D(1))
      EIG(2)=CDEXP(DCMPLX(-AI(INDX(2)),AR(INDX(2)))*D(1))
      EIG(3)=CDEXP(DCMPLX(-AI(INDX(3)),AR(INDX(3)))*D(1))
      EIG(4)=CDEXP(DCMPLX(-AI(INDX(4)),AR(INDX(4)))*D(1))

      DO 80 I=1,4
	DO 80 J=1,4
	  M(J,I,1)=DCMPLX(UR(J,INDX(I)),UI(J,INDX(I)))
 80   CONTINUE
 
      DO 180 ILAY=1,NM-1
 
	IUSE=MOD(ILAY-1,2)+1
	IFILL=MOD(ILAY+2,2)+1
 
 
C     Calculate eigenvalues and eigenvectors

	DO 110 I=1,4
	  AR(I)=DREAL(ZQ(I,ILAY+1))
	  AI(I)=DIMAG(ZQ(I,ILAY+1))
	  DO 1051 J=1,4
	    UR(J,I)=DREAL(ZPHI(J,I,ILAY+1))
	    UI(J,I)=DIMAG(ZPHI(J,I,ILAY+1))
 1051     CONTINUE
 110    CONTINUE

C     Sort out into p and s polarisations

	IF (ILAY.EQ.NM-1) THEN
          INDX(1)=3
          INDX(2)=1
          INDX(3)=4
          INDX(4)=2
        ELSE
	  CALL SORT(AI,4,INDX)
	  CALL SORT2(AR,AI,4,INDX)
        END IF
 
	DO 120 I=1,4
	  DO 120 J=1,4
	    M(J,I,IFILL)=DCMPLX(UR(J,INDX(I)),UI(J,INDX(I)))
	    M2(J,I)=M(J,I,IFILL)
 120    CONTINUE

	DO 130 I=1,4
	  DO 130 J=1,4
	    M(J,I,IUSE)=M(J,I,IUSE)*EIG(I)
	    M1(J,I)=M(J,I,IUSE)
 130    CONTINUE

C     Calculate transfer matrix A

	CALL INVMAT(M1,A,W,ERRFLAG)
        IF (ERRFLAG.NE.0) RETURN
	CALL CMATX1(W,M2,A,4)

C     Perform S12*A21
 
	W(1,1)=S(1,3)*A(3,1)+S(1,4)*A(4,1)
	W(2,1)=S(1,3)*A(3,2)+S(1,4)*A(4,2)
	W(3,1)=S(2,3)*A(3,1)+S(2,4)*A(4,1)
	W(4,1)=S(2,3)*A(3,2)+S(2,4)*A(4,2)

C     Perform A11-S12*A21

 
	W(1,2)=A(1,1)-W(1,1)
	W(2,2)=A(1,2)-W(2,1)
	W(3,2)=A(2,1)-W(3,1)
	W(4,2)=A(2,2)-W(4,1)

C     Perform (A11-S12*A21)**-1

	DET=1.0D0/(W(1,2)*W(4,2)-W(2,2)*W(3,2))
	W(1,3)=DET*W(4,2)
	W(2,3)=-DET*W(2,2)
	W(3,3)=-DET*W(3,2)
	W(4,3)=DET*W(1,2)

C     Perform S12*A22-A12
 
	W(1,2)=S(1,3)*A(3,3)+S(1,4)*A(4,3)-A(1,3)
	W(2,2)=S(1,3)*A(3,4)+S(1,4)*A(4,4)-A(1,4)
	W(3,2)=S(2,3)*A(3,3)+S(2,4)*A(4,3)-A(2,3)
	W(4,2)=S(2,3)*A(3,4)+S(2,4)*A(4,4)-A(2,4)

C     Perform S11(0,n+1)=((A11-S12*A21)**-1)*S11)
 
	SNEW(1,1)=W(1,3)*S(1,1)+W(2,3)*S(2,1)
	SNEW(1,2)=W(1,3)*S(1,2)+W(2,3)*S(2,2)
	SNEW(2,1)=W(3,3)*S(1,1)+W(4,3)*S(2,1)
	SNEW(2,2)=W(3,3)*S(1,2)+W(4,3)*S(2,2)

C     Perform S12(0,n+1)=((A11-S12*A21)**-1)(S12*A22-A12)
 
	SNEW(1,3)=W(1,3)*W(1,2)+W(2,3)*W(3,2)
	SNEW(1,4)=W(1,3)*W(2,2)+W(2,3)*W(4,2)
	SNEW(2,3)=W(3,3)*W(1,2)+W(4,3)*W(3,2)
	SNEW(2,4)=W(3,3)*W(2,2)+W(4,3)*W(4,2)

C     Perform S22*A21
 
	W(1,1)=S(3,3)*A(3,1)+S(3,4)*A(4,1)
	W(2,1)=S(3,3)*A(3,2)+S(3,4)*A(4,2)
	W(3,1)=S(4,3)*A(3,1)+S(4,4)*A(4,1)
	W(4,1)=S(4,3)*A(3,2)+S(4,4)*A(4,2)
 
C     Perform S22*A21*S11(0,n+1)+S21

	SNEW(3,1)=W(1,1)*SNEW(1,1)+W(2,1)*SNEW(2,1)+S(3,1)
	SNEW(3,2)=W(1,1)*SNEW(1,2)+W(2,1)*SNEW(2,2)+S(3,2)
	SNEW(4,1)=W(3,1)*SNEW(1,1)+W(4,1)*SNEW(2,1)+S(4,1)
	SNEW(4,2)=W(3,1)*SNEW(1,2)+W(4,1)*SNEW(2,2)+S(4,2)

C     Perform S22*A21*S12(0,n+1)+S22*A22
 
	SNEW(3,3)=W(1,1)*SNEW(1,3)+W(2,1)*SNEW(2,3)+
     $  S(3,3)*A(3,3)+S(3,4)*A(4,3)
	SNEW(3,4)=W(1,1)*SNEW(1,4)+W(2,1)*SNEW(2,4)+
     $  S(3,3)*A(3,4)+S(3,4)*A(4,4)
	SNEW(4,3)=W(3,1)*SNEW(1,3)+W(4,1)*SNEW(2,3)+
     $  S(4,3)*A(3,3)+S(4,4)*A(4,3)
	SNEW(4,4)=W(3,1)*SNEW(1,4)+W(4,1)*SNEW(2,4)+
     $  S(4,3)*A(3,4)+S(4,4)*A(4,4)

	EIG(1)=CDEXP(DCMPLX(-AI(INDX(1)),AR(INDX(1)))*D(ILAY+1))
	EIG(2)=CDEXP(DCMPLX(-AI(INDX(2)),AR(INDX(2)))*D(ILAY+1))
	EIG(3)=CDEXP(DCMPLX(-AI(INDX(3)),AR(INDX(3)))*D(ILAY+1))
	EIG(4)=CDEXP(DCMPLX(-AI(INDX(4)),AR(INDX(4)))*D(ILAY+1))

	DO 140 I=1,4
	  DO 140 J=1,4
	    S(J,I)=SNEW(J,I)
 140    CONTINUE
 
 180  CONTINUE
 
      DO 200 I=1,2
	DO 200 J=1,2
	  RMAT(J,I)=S(J+2,I)
          TMAT(J,I)=S(J,I)
 200  CONTINUE
 
      RETURN
      END

*     SMAT ends
