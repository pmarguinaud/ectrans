! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE SPNSDE1_MOD
CONTAINS
SUBROUTINE SPNSDE1(KM,PEPSNM,PIA,YDSP,YDSP_NS)

USE PARKIND1  ,ONLY : JPIM     ,JPRB

!USE TPM_GEN
USE TPM_DIM         ,ONLY : R
USE TPM_FIELDS      ,ONLY : F
USE POINTER_MOD
!USE TPM_TRANS


!**** *SPNSDE1* - Compute North-South derivative in spectral space

!     Purpose.
!     --------
!        In Laplace space compute the the North-south derivative

!**   Interface.
!     ----------
!        CALL SPNSDE1(...)

!        Explicit arguments :
!        --------------------
!        KM -zonal wavenumber (input-c)
!        PEPSNM - REPSNM for wavenumber KM (input-c)
!        PF  (NLEI1,2*KF_SCALARS) - input field (input)
!        PNSD(NLEI1,2*KF_SCALARS) - N-S derivative (output)

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

!        Implicit arguments :  YOMLAP
!        --------------------

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
!        Original : 00-02-01 From SPNSDE1 in IFS CY22R1

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)    :: KM
REAL(KIND=JPRB),    INTENT(IN)    :: PEPSNM(0:R%NTMAX+2)
REAL(KIND=JPRB),    INTENT(INOUT) :: PIA(:,:)
TYPE (PTRS)                       :: YDSP    (:)
TYPE (PTRS)                       :: YDSP_NS (:)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: IJ, ISKIP, J, JN,JI,ISMAX
INTEGER(KIND=JPIM) :: INF
INTEGER(KIND=JPIM) :: IR, II, IR_NS, II_NS
REAL(KIND=JPRB) :: ZEPSNM(-1:R%NSMAX+4)
REAL(KIND=JPRB) :: ZN(-1:R%NTMAX+4)


!     ------------------------------------------------------------------

!*       1.    COMPUTE NORTH SOUTH DERIVATIVE.
!              -------------------------------


!*       1.1      COMPUTE

ISMAX = R%NSMAX
DO JN=KM-1,ISMAX+2
  IJ = ISMAX+3-JN
  ZN(IJ) = F%RN(JN)
  IF( JN >= 0 ) ZEPSNM(IJ) = PEPSNM(JN)
ENDDO
ZN(0) = F%RN(ISMAX+3)

IF (SIZE (YDSP) /= SIZE (YDSP_NS)) THEN
  WRITE (0, *) __FILE__, ':', __LINE__, SIZE (YDSP)
  WRITE (0, *) __FILE__, ':', __LINE__, SIZE (YDSP_NS)
  STOP 1
ENDIF

INF = SIZE (YDSP)

IF (KM == 0) THEN

  DO J=1,INF

    IR    = 2*YDSP    (J)%IPTR-1
    IR_NS = 2*YDSP_NS (J)%IPTR-1

    DO JI=2,ISMAX+3-KM
      PIA(JI,IR_NS) = -ZN(JI+1)*ZEPSNM(JI)*PIA(JI+1,IR)+ZN(JI-2)*ZEPSNM(JI-1)*PIA(JI-1,IR)
    ENDDO
  ENDDO

ELSE

  DO J=1,INF

    IR    = 2*YDSP    (J)%IPTR-1; II    = IR    + 1
    IR_NS = 2*YDSP_NS (J)%IPTR-1; II_NS = IR_NS + 1

    DO JI=2,ISMAX+3-KM
      PIA(JI,IR_NS) = -ZN(JI+1)*ZEPSNM(JI)*PIA(JI+1,IR)+ZN(JI-2)*ZEPSNM(JI-1)*PIA(JI-1,IR)
      PIA(JI,II_NS) = -ZN(JI+1)*ZEPSNM(JI)*PIA(JI+1,II)+ZN(JI-2)*ZEPSNM(JI-1)*PIA(JI-1,II)
    ENDDO
  ENDDO

ENDIF


!     ------------------------------------------------------------------

END SUBROUTINE SPNSDE1
END MODULE SPNSDE1_MOD
