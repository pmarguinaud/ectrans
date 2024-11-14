! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FTINV1_CTL_MOD
CONTAINS
SUBROUTINE FTINV1_CTL(KF_GP,KF_FS,KF_OUT_LT,YDGP)


!**** *FTINV1_CTL - Inverse Fourier transform control

!     Purpose. Control routine for Fourier to Gridpoint transform
!     --------

!**   Interface.
!     ----------
!        CALL FTINV1_CTL(..)

!        Explicit arguments :
!        --------------------
!        PGP     -  gridpoint array
!        KF_UV_G      - global number of spectral u-v fields
!        KF_SCALARS_G - global number of scalar spectral fields
!        KF_UV        - local number of spectral u-v fields
!        KF_SCALARS   - local number of scalar spectral fields
!        KF_SCDERS    - local number of derivatives of scalar spectral fields
!        KF_GP        - total number of output gridpoint fields
!        KF_FS        - total number of fields in fourier space
!        KF_OUT_LT    - total number of fields coming out from inverse LT
!        KVSETUV - "B"  set in spectral/fourier space for
!                   u and v variables
!        KVSETSC - "B" set in spectral/fourier space for
!                  scalar variables
!        KPTRGP - pointer array to fi3elds in gridpoint space

!     Method.
!     -------

!     Externals.  TRLTOG1      - transposition routine
!     ----------  FOURIER_IN  - copy fourier data from Fourier buffer
!                 FTINV       - fourier transform
!                 FSC         - Fourier space computations

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!        R. El Khatib 09-Sep-2020 NSTACK_MEMORY_TR

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_GEN         ,ONLY : NERR   ,NSTACK_MEMORY_TR
USE TPM_TRANS       ,ONLY : FOUBUF, LDIVGP, LSCDERS, LUVDER, LVORGP,LATLON
USE TPM_DISTR       ,ONLY : D, MYPROC, NPROC, MYSETV
USE TPM_FLT         ,ONLY : S
USE FOURIER_IN_MOD  ,ONLY : FOURIER_IN
USE FSC1_MOD        ,ONLY : FSC1
USE FTINV_MOD       ,ONLY : FTINV
USE TRLTOG1_MOD     ,ONLY : TRLTOG1
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
USE POINTER_MOD
!

IMPLICIT NONE

INTEGER(KIND=JPIM) ,INTENT(IN) :: KF_GP
INTEGER(KIND=JPIM) ,INTENT(IN) :: KF_FS
INTEGER(KIND=JPIM) ,INTENT(IN) :: KF_OUT_LT
TYPE (PTRG)                    :: YDGP (:)

INTEGER(KIND=JPIM) :: J
INTEGER(KIND=JPIM) :: IST
INTEGER(KIND=JPIM) :: JGL,IGL
INTEGER(KIND=JPIM) :: IBEG,IEND,IINC

REAL(KIND=JPRB),TARGET  :: ZGTF_STACK(KF_FS*MIN(1,MAX(0,NSTACK_MEMORY_TR)),D%NLENGTF)
REAL(KIND=JPRB),TARGET, ALLOCATABLE :: ZGTF_HEAP(:,:)
REAL(KIND=JPRB),POINTER  :: ZGTF(:,:)

TYPE (PTRS), ALLOCATABLE :: YLGP_U      (:)
TYPE (PTRS), ALLOCATABLE :: YLGP_V      (:)
TYPE (PTRS), ALLOCATABLE :: YLGP_UV     (:)
TYPE (PTRS), ALLOCATABLE :: YLGP_U_EW   (:)
TYPE (PTRS), ALLOCATABLE :: YLGP_V_EW   (:)
TYPE (PTRS), ALLOCATABLE :: YLGP_UV_EW  (:)
TYPE (PTRS), ALLOCATABLE :: YLGP_S      (:)
TYPE (PTRS), ALLOCATABLE :: YLGP_S_EW   (:)
TYPE (PTRS), ALLOCATABLE :: YLGP_S_NS   (:)

INTEGER (KIND=JPIM) :: IF_FS

!     ------------------------------------------------------------------

!    1.  Copy Fourier data to local array

CALL GSTATS(107,0)

IF (NSTACK_MEMORY_TR == 1) THEN
  ZGTF => ZGTF_STACK(:,:)
ELSE
  ALLOCATE(ZGTF_HEAP(KF_FS,D%NLENGTF))
! Now, force the OS to allocate this shared array right now, not when it starts
! to be used which is an OPEN-MP loop, that would cause a threads synchronization lock :
  IF (KF_FS > 0 .AND. D%NLENGTF > 0) THEN
    ZGTF_HEAP(1,1)=HUGE(1._JPRB)
  ENDIF
  ZGTF => ZGTF_HEAP(:,:)
ENDIF

YLGP_U    = GREP (NTYPE_U   )
YLGP_V    = GREP (NTYPE_V   )
YLGP_S    = GREP (NTYPE_S   )
YLGP_U_EW = GREP (NTYPE_U_EW)
YLGP_V_EW = GREP (NTYPE_V_EW)
YLGP_S_EW = GREP (NTYPE_S_EW)
YLGP_S_NS = GREP (NTYPE_S_NS)

YLGP_UV = [YLGP_U, YLGP_V]
YLGP_UV_EW = [YLGP_U_EW, YLGP_V_EW]

IF (MYPROC > NPROC/2)THEN
  IBEG=1
  IEND=D%NDGL_FS
  IINC=1
ELSE
  IBEG=D%NDGL_FS
  IEND=1
  IINC=-1
ENDIF

CALL GSTATS(1639,0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JGL,IGL)
DO JGL=IBEG,IEND,IINC
  IGL = JGL
  CALL FOURIER_IN(ZGTF,KF_OUT_LT,IGL)

!    2.  Fourier space computations

  IF (SIZE (YLGP_U) > 0 .OR. SIZE (YLGP_S_NS) > 0 .OR. (LATLON.AND.S%LDLL) ) THEN
    CALL FSC1(IGL,YLGP_UV,YLGP_S,YLGP_S_NS,YLGP_S_EW,YLGP_UV_EW)
  ENDIF

!    3.  Fourier transform

  IF (KF_FS > 0) THEN
    CALL FTINV(ZGTF,KF_FS,IGL) 
  ENDIF
ENDDO
!$OMP END PARALLEL DO
CALL GSTATS(1639,1)

CALL GSTATS(107,1)

CALL GSTATS(157,0)
CALL TRLTOG1(ZGTF,KF_FS,KF_GP,YDGP)
CALL GSTATS(157,1)

!     ------------------------------------------------------------------

CONTAINS

FUNCTION GREP (KTYPE) RESULT (YLGP)

TYPE (PTRS), ALLOCATABLE :: YLGP (:)
INTEGER :: KTYPE

INTEGER :: ICOUNT, J

ICOUNT = COUNT ((YDGP%ITYPE == KTYPE) .AND. (MYSETV == YDGP%IVSET))

ALLOCATE (YLGP (ICOUNT))

ICOUNT = 0

DO J = 1, SIZE (YDGP)
  IF ((YDGP (J)%ITYPE == KTYPE) .AND. (YDGP (J)%IVSET == MYSETV)) THEN
    ICOUNT = ICOUNT + 1
    YLGP (ICOUNT)%P => ZGTF (YDGP (J)%IRANK_L, :)
    YLGP (ICOUNT)%ITYPE = KTYPE
  ENDIF
ENDDO

END FUNCTION

END SUBROUTINE FTINV1_CTL
END MODULE FTINV1_CTL_MOD
