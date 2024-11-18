! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LTINV1_MOD
CONTAINS
SUBROUTINE LTINV1(KM,KMLOC,KF_OUT_LT,KLEI2,KDIM1,&
 & YDSP,YDSP_VOR,YDSP_DIV,YDSP_U, &
 & YDSP_V,YDSP_S,YDSP_S_NS,FSPGL_PROC)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_FIELDS      ,ONLY : F
USE TPM_DIM         ,ONLY : R
USE TPM_TRANS       ,ONLY : LDIVGP, LVORGP, NF_SC2, NF_SC3A, NF_SC3B, LATLON
USE TPM_FLT
USE TPM_GEOMETRY

!USE PRLE1_MOD
USE PREPSNM_MOD     ,ONLY : PREPSNM
USE PRFI1B_MOD      ,ONLY : PRFI1B
USE PRFI1B1_MOD     ,ONLY : PRFI1B1
USE VDTUV_MOD       ,ONLY : VDTUV
USE VDTUV1_MOD      ,ONLY : VDTUV1
USE SPNSDE_MOD      ,ONLY : SPNSDE
USE SPNSDE1_MOD     ,ONLY : SPNSDE1
USE LEINV_MOD       ,ONLY : LEINV
USE ASRE1B_MOD      ,ONLY : ASRE1B
USE FSPGL_INT_MOD   ,ONLY : FSPGL_INT
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
USE CDMAP_MOD       ,ONLY : CDMAP
USE POINTER_MOD


!**** *LTINV1* - Inverse Legendre transform
!
!     Purpose.
!     --------
!        Tranform from Laplace space to Fourier space, compute U and V
!        and north/south derivatives of state variables.

!**   Interface.
!     ----------
!        *CALL* *LTINV1(...)

!        Explicit arguments :
!        --------------------
!          KM        - zonal wavenumber
!          KMLOC     - local zonal wavenumber
!          PSPVOR    - spectral vorticity
!          PSPDIV    - spectral divergence
!          PSPSCALAR - spectral scalar variables

!        Implicit arguments :  The Laplace arrays of the model.
!        --------------------  The values of the Legendre polynomials
!                              The grid point arrays of the model
!     Method.
!     -------

!     Externals.
!     ----------

!         PREPSNM - prepare REPSNM for wavenumber KM
!         PRFI1B  - prepares the spectral fields
!         VDTUV   - compute u and v from vorticity and divergence
!         SPNSDE  - compute north-south derivatives
!         LEINV   - Inverse Legendre transform
!         ASRE1   - recombination of symmetric/antisymmetric part

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Temperton, 1991, MWR 119 p1303

!     Author.
!     -------
!        Mats Hamrud  *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-01 From LTINV1 in IFS CY22R1
!     ------------------------------------------------------------------

IMPLICIT NONE


INTEGER(KIND=JPIM), INTENT(IN) :: KM
INTEGER(KIND=JPIM), INTENT(IN) :: KMLOC
INTEGER(KIND=JPIM), INTENT(IN) :: KF_OUT_LT
INTEGER(KIND=JPIM), INTENT(IN) :: KLEI2
INTEGER(KIND=JPIM), INTENT(IN) :: KDIM1
TYPE (PTRS)                    :: YDSP      (:)
TYPE (PTRS)                    :: YDSP_VOR  (:)
TYPE (PTRS)                    :: YDSP_DIV  (:)
TYPE (PTRS)                    :: YDSP_U    (:)
TYPE (PTRS)                    :: YDSP_V    (:)
TYPE (PTRS)                    :: YDSP_S    (:)
TYPE (PTRS)                    :: YDSP_S_NS (:)

EXTERNAL  FSPGL_PROC
OPTIONAL  FSPGL_PROC


REAL(KIND=JPRB) :: ZIA(R%NLEI1,KLEI2)
REAL(KIND=JPRB) :: ZEPSNM(0:R%NTMAX+2)
REAL(KIND=JPRB), ALLOCATABLE :: ZSOA1(:,:), ZAOA1(:,:), ZALN(:,:)

INTEGER(KIND=JPIM) :: IFC, ISTA, IIFC, IDGLU, JGL, JFLD
INTEGER(KIND=JPIM) :: IVORL,IVORU,IDIVL,IDIVU,IUL,IUU,IVL,IVU,ISL,ISLO,ISU,IDL,IDU, IGLS
INTEGER(KIND=JPIM) :: IFIRST, ILAST, IDIM1,IDIM3,J3
INTEGER(KIND=JPIM) :: INSDS, INSDE, IUVS, IUVE, IST
INTEGER(KIND=JPIM) :: IF_UV

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.       PERFORM LEGENDRE TRANFORM.
!                 --------------------------

IF (LHOOK) CALL DR_HOOK('LTINV1_MOD',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       1.    PREPARE ZEPSNM.
!              ---------------

CALL PREPSNM(KM,KMLOC,ZEPSNM)

!     ------------------------------------------------------------------


!*       3.    SPECTRAL COMPUTATIONS FOR U,V AND DERIVATIVES.
!              ----------------------------------------------

IF_UV = SIZE (YDSP_U)

IFIRST = 1
ILAST  = 4*IF_UV

CALL PRFI1B1(KM,ZIA,YDSP)

CALL VDTUV1(KM,ZEPSNM,ZIA,YDSP_VOR,YDSP_DIV,YDSP_U,YDSP_V)
IF (SIZE (YDSP_S_NS) > 0) THEN
  CALL SPNSDE1(KM,ZEPSNM,ZIA,YDSP_S,YDSP_S_NS)
ENDIF

!*       4.    INVERSE LEGENDRE TRANSFORM.
!              ---------------------------


ISTA = 1
IFC  = 2*KF_OUT_LT
IF(IF_UV > 0 .AND. .NOT. LVORGP) THEN
  ISTA = ISTA+2*IF_UV
ENDIF
IF(IF_UV > 0 .AND. .NOT. LDIVGP) THEN
  ISTA = ISTA+2*IF_UV
ENDIF

IIFC=IFC
IF(KM == 0)THEN
  IIFC=IFC/2
ENDIF

IF( LATLON.AND.S%LDLL ) THEN

  IDGLU = MIN(R%NDGNH,G%NDGLU(KM))

  IF( (S%LSHIFTLL .AND. KM < 2*IDGLU) .OR.&
   & (.NOT.S%LSHIFTLL .AND. KM < 2*(IDGLU-1)) ) THEN

    ISL = MAX(R%NDGNH-G%NDGLU(KM)+1,1)
    ISLO = S%FA(KMLOC)%ISLD

    ALLOCATE(ZAOA1(KDIM1,R%NLEI3))
    ALLOCATE(ZSOA1(KDIM1,R%NLEI3))
    CALL LEINV(KM,KMLOC,IFC,IIFC,KF_OUT_LT,ISL,IDGLU,ZIA(:,ISTA:ISTA+IFC-1),ZAOA1,ZSOA1)

!*       5.    RECOMBINATION SYMMETRIC/ANTISYMMETRIC PART.
! before (non-linear) mapping !!!!

    ALLOCATE( ZALN(KDIM1,2*R%NDGNH) )
    DO JGL=ISL, R%NDGNH
      IGLS = 2*R%NDGNH+1-JGL
      DO JFLD=1,2*KF_OUT_LT
        ZALN(JFLD, JGL)  =  ZSOA1(JFLD,JGL)+ZAOA1(JFLD,JGL)
        ZALN(JFLD, IGLS) =  ZSOA1(JFLD,JGL)-ZAOA1(JFLD,JGL)
      ENDDO
    ENDDO
    
    IF (SIZE (YDSP_U) > 0 .OR. SIZE (YDSP_S_NS) > 0) THEN
      IST = 1
      IF (LVORGP) IST = IST + 2*SIZE (YDSP_VOR)
      IF (LDIVGP) IST = IST + 2*SIZE (YDSP_DIV)
      IUVS  = IST
      IUVE  = IUVS + 2*(SIZE (YDSP_U) + SIZE (YDSP_V)) - 1
      INSDS = IUVE + 2*(SIZE (YDSP_S)) + 1
      INSDE = INSDS + 2*(SIZE (YDSP_S_NS)) - 1
      
      IGLS = 2*R%NDGNH - ISL + 1
      IF (SIZE (YDSP_U) > 0) THEN
        DO JGL=ISL, IGLS
          DO JFLD=IUVS,IUVE
            ZALN(JFLD, JGL)  =  ZALN(JFLD,JGL)*F%RACTHE(JGL)
          ENDDO
        ENDDO
      ENDIF
      IF (SIZE (YDSP_S_NS) > 0) THEN
        DO JGL=ISL, IGLS
          DO JFLD=INSDS,INSDE
            ZALN(JFLD, JGL)  =  ZALN(JFLD,JGL)*F%RACTHE(JGL)
          ENDDO
        ENDDO
      ENDIF
    ENDIF
    
    DEALLOCATE(ZAOA1)
    DEALLOCATE(ZSOA1)
    
    ! this routine maps to the output latitudes AND fills the FOUBUF
    CALL CDMAP(KM,KMLOC,ISL,ISLO,ZEPSNM(R%NTMAX+1),-1_JPIM,&
     & R%NDGNH,S%NDGNHD,2*KF_OUT_LT,ZALN,ZALN)
    DEALLOCATE(ZALN)
    
  ENDIF

ELSE
  IDGLU = MIN(R%NDGNH,G%NDGLU(KM))
  ISL = MAX(R%NDGNH-G%NDGLU(KM)+1,1)

  ALLOCATE(ZAOA1(KDIM1,R%NLEI3))
  ALLOCATE(ZSOA1(KDIM1,R%NLEI3))
  CALL LEINV(KM,KMLOC,IFC,IIFC,KF_OUT_LT,ISL,IDGLU,ZIA(:,ISTA:ISTA+IFC-1),ZAOA1,ZSOA1)

!     ------------------------------------------------------------------

!*       5.    RECOMBINATION SYMMETRIC/ANTISYMMETRIC PART/FILL FOUBUF
!              --------------------------------------------

  CALL ASRE1B(KF_OUT_LT,KM,KMLOC,ZAOA1,ZSOA1)
  DEALLOCATE(ZAOA1)
  DEALLOCATE(ZSOA1)
  
ENDIF

!     ------------------------------------------------------------------

!     6. OPTIONAL COMPUTATIONS IN FOURIER SPACE

IF(PRESENT(FSPGL_PROC)) THEN
  STOP 1
ENDIF


IF (LHOOK) CALL DR_HOOK('LTINV1_MOD',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------
END SUBROUTINE
END MODULE LTINV1_MOD




