! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE VDTUV1_MOD
CONTAINS
SUBROUTINE VDTUV1(KM,PEPSNM,PIA,YDSP_VOR,YDSP_DIV,YDSP_U,YDSP_V)

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DIM         ,ONLY : R
USE TPM_FIELDS      ,ONLY : F
USE POINTER_MOD


!**** *VDTUV1* - Compute U,V in  spectral space

!     Purpose.
!     --------
!        In Laplace space compute the the winds
!        from vorticity and divergence.

!**   Interface.
!     ----------
!        CALL VDTUV1(...)

!        Explicit arguments :  KM -zonal wavenumber (input-c)
!        --------------------  KFIELD - number of fields (input-c)
!                              PEPSNM - REPSNM for wavenumber KM (input-c)
!                              PVOR(NLEI1,2*KFIELD) - vorticity (input)
!                              PDIV(NLEI1,2*KFIELD) - divergence (input)
!                              PU(NLEI1,2*KFIELD)   - u wind (output)
!                              PV(NLEI1,2*KFIELD)   - v wind (output)
!        Organisation within NLEI1:
!        NLEI1 = NSMAX+4+mod(NSMAX+4+1,2)
!                        overdimensioning
!        1        : n=NSMAX+2
!        2        : n=NSMAX+1
!        3        : n=NSMAX
!        .        :
!        .        :
!        NSMAX+3  : n=0
!        NSMAX+4  : n=-1

!        Implicit arguments :  Eigenvalues of inverse Laplace operator
!        --------------------  from YOMLAP

!     Method.
!     -------

!     Externals.   None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Temperton, 1991, MWR 119 p1303

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-01 From VDTUV1 in IFS CY22R1

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KM
REAL(KIND=JPRB), INTENT(IN)    :: PEPSNM(0:R%NTMAX+2)
REAL (KIND=JPRB), INTENT (INOUT) :: PIA (:,:)
TYPE (PTRS) :: YDSP_VOR (:)
TYPE (PTRS) :: YDSP_DIV (:)
TYPE (PTRS) :: YDSP_U (:)
TYPE (PTRS) :: YDSP_V (:)


!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: II, IJ, IR, J, JN, ISMAX,JI
INTEGER(KIND=JPIM) :: IR_VOR, IR_DIV, IR_U, IR_V
INTEGER(KIND=JPIM) :: II_VOR, II_DIV, II_U, II_V
INTEGER(KIND=JPIM) :: INF

!     LOCAL REAL SCALARS
REAL(KIND=JPRB) :: ZKM
REAL(KIND=JPRB) :: ZN(-1:R%NTMAX+4)
REAL(KIND=JPRB) :: ZLAPIN(-1:R%NSMAX+4)
REAL(KIND=JPRB) :: ZEPSNM(-1:R%NSMAX+4)



!     ------------------------------------------------------------------

!*       1.    COMPUTE U V FROM VORTICITY AND DIVERGENCE.
!              ------------------------------------------

ZKM = KM
ISMAX = R%NSMAX
DO JN=KM-1,ISMAX+2
  IJ = ISMAX+3-JN
  ZN(IJ) = F%RN(JN)
  ZLAPIN(IJ) = F%RLAPIN(JN)
  IF( JN >= 0 ) ZEPSNM(IJ) = PEPSNM(JN)
ENDDO
ZN(0) = F%RN(ISMAX+3)

IF (SIZE (YDSP_VOR) /= SIZE (YDSP_DIV)) STOP 1
IF (SIZE (YDSP_U) /= SIZE (YDSP_V)) STOP 1
IF (SIZE (YDSP_VOR) /= SIZE (YDSP_U)) STOP 1

INF = SIZE (YDSP_VOR)

!*       1.1      U AND V (KM=0) .

IF(KM == 0) THEN
  DO J=1,INF

    IR_VOR = 2*YDSP_VOR (J)%IPTR-1
    IR_DIV = 2*YDSP_DIV (J)%IPTR-1
    IR_U   = 2*YDSP_U   (J)%IPTR-1
    IR_V   = 2*YDSP_V   (J)%IPTR-1

    DO JI=2,ISMAX+3-KM
      PIA(JI,IR_U) = +&
       &ZN(JI+1)*ZEPSNM(JI)*ZLAPIN(JI+1)*PIA(JI+1,IR_VOR)-&
       &ZN(JI-2)*ZEPSNM(JI-1)*ZLAPIN(JI-1)*PIA(JI-1,IR_VOR)
      PIA(JI,IR_V) = -&
       &ZN(JI+1)*ZEPSNM(JI)*ZLAPIN(JI+1)*PIA(JI+1,IR_DIV)+&
       &ZN(JI-2)*ZEPSNM(JI-1)*ZLAPIN(JI-1)*PIA(JI-1,IR_DIV)
    ENDDO
  ENDDO

!*       1.2      U AND V (KM!=0) .

ELSE
  DO J=1,INF

    IR_VOR = 2*YDSP_VOR (J)%IPTR-1; II_VOR = IR_VOR + 1
    IR_DIV = 2*YDSP_DIV (J)%IPTR-1; II_DIV = IR_DIV + 1
    IR_U   = 2*YDSP_U   (J)%IPTR-1; II_U   = IR_U   + 1
    IR_V   = 2*YDSP_V   (J)%IPTR-1; II_V   = IR_V   + 1

    DO JI=2,ISMAX+3-KM
      PIA(JI,IR_U) = -ZKM*ZLAPIN(JI)*PIA(JI,II_DIV)+&
       &ZN(JI+1)*ZEPSNM(JI)*ZLAPIN(JI+1)*PIA(JI+1,IR_VOR)-&
       &ZN(JI-2)*ZEPSNM(JI-1)*ZLAPIN(JI-1)*PIA(JI-1,IR_VOR)
      PIA(JI,II_U) = +ZKM*ZLAPIN(JI)*PIA(JI,IR_DIV)+&
       &ZN(JI+1)*ZEPSNM(JI)*ZLAPIN(JI+1)*PIA(JI+1,II_VOR)-&
       &ZN(JI-2)*ZEPSNM(JI-1)*ZLAPIN(JI-1)*PIA(JI-1,II_VOR)
      PIA(JI,IR_V) = -ZKM*ZLAPIN(JI)*PIA(JI,II_VOR)-&
       &ZN(JI+1)*ZEPSNM(JI)*ZLAPIN(JI+1)*PIA(JI+1,IR_DIV)+&
       &ZN(JI-2)*ZEPSNM(JI-1)*ZLAPIN(JI-1)*PIA(JI-1,IR_DIV)
      PIA(JI,II_V) = +ZKM*ZLAPIN(JI)*PIA(JI,IR_VOR)-&
       &ZN(JI+1)*ZEPSNM(JI)*ZLAPIN(JI+1)*PIA(JI+1,II_DIV)+&
       &ZN(JI-2)*ZEPSNM(JI-1)*ZLAPIN(JI-1)*PIA(JI-1,II_DIV)
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE VDTUV1
END MODULE VDTUV1_MOD

