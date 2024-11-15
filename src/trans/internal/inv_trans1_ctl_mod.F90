! (C) Copyright 2001- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE INV_TRANS1_CTL_MOD
CONTAINS
SUBROUTINE INV_TRANS1_CTL (YDGP, YDSP, FSPGL_PROC)

!**** *INV_TRANS1_CTL* - Control routine for inverse spectral transform.

!     Purpose.
!     --------
!        Control routine for the inverse spectral transform

!**   Interface.
!     ----------
!     CALL INV_TRANS1_CTL(...)

!     Explicit arguments :
!     --------------------
!     KF_UV_G      - global number of spectral u-v fields
!     KF_SCALARS_G - global number of scalar spectral fields
!     KF_GP        - total number of output gridpoint fields
!     KF_FS        - total number of fields in fourier space
!     KF_OUT_LT    - total number of fields coming out from inverse LT
!     KF_UV        - local number of spectral u-v fields
!     KF_SCALARS   - local number of scalar spectral fields
!     KF_SCDERS    - local number of derivatives of scalar spectral fields
!     PSPVOR(:,:)  - spectral vorticity (input)
!     PSPDIV(:,:)  - spectral divergence (input)
!     PSPSCALAR(:,:) - spectral scalarvalued fields (input)
!     KVSETUV(:)  - indicating which 'b-set' in spectral space owns a
!                   vor/div field. Equivalant to NBSETLEV in the IFS.
!                   The length of KVSETUV should be the GLOBAL number
!                   of u/v fields which is the dimension of u and v releated
!                   fields in grid-point space.
!     KVESETSC(:) - indicating which 'b-set' in spectral space owns a
!                   scalar field. As for KVSETUV this argument is required
!                   if the total number of processors is greater than
!                   the number of processors used for distribution in
!                   spectral wave space.
!     FSPGL_PROC  - external procedure to be executed in fourier space
!                   before transposition
!     PGP(:,:,:)  - gridpoint fields (output)

!                  The ordering of the output fields is as follows (all
!                  parts are optional depending on the input switches):

!       vorticity     : KF_UV_G fields
!       divergence    : KF_UV_G fields
!       u             : KF_UV_G fields
!       v             : KF_UV_G fields
!       scalar fields : KF_SCALARS_G fields
!       N-S derivative of scalar fields : KF_SCALARS_G fields
!       E-W derivative of u : KF_UV_G fields
!       E-W derivative of v : KF_UV_G fields
!       E-W derivative of scalar fields : KF_SCALARS_G fields

!     Method.
!     -------

!     Externals.  SHUFFLE     - reshuffle fields for load balancing
!     ----------  FIELD_SPLIT - split fields in NPROMATR packets
!                 LTINV_CTL   - control of Legendre transform
!                 FTINV1_CTL   - control of Fourier transform

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 01-01-03

!     ------------------------------------------------------------------


USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_GEN         ,ONLY : NPROMATR
USE TPM_TRANS       ,ONLY : LDIVGP, LSCDERS, LUVDER, LVORGP
USE TPM_DISTR       ,ONLY : MYPROC, MYSETV

USE SHUFFLE_MOD     ,ONLY : SHUFFLE
USE FIELD_SPLIT_MOD ,ONLY : FIELD_SPLIT
USE LTINV1_CTL_MOD  ,ONLY : LTINV1_CTL
USE FTINV1_CTL_MOD  ,ONLY : FTINV1_CTL
USE POINTER_MOD
!

IMPLICIT NONE

TYPE (PTRG) :: YDGP (:)
TYPE (PTRS) :: YDSP (:)
EXTERNAL  FSPGL_PROC
OPTIONAL  FSPGL_PROC

INTEGER(KIND=JPIM) :: IF_GP
INTEGER(KIND=JPIM) :: IF_FS
INTEGER(KIND=JPIM) :: IF_OUT_LT

IF_OUT_LT = SIZE (YDSP)
IF (.NOT. LVORGP) IF_OUT_LT = IF_OUT_LT - COUNT (YDSP%ITYPE == NTYPE_VOR)
IF (.NOT. LDIVGP) IF_OUT_LT = IF_OUT_LT - COUNT (YDSP%ITYPE == NTYPE_DIV)

IF_GP = SIZE (YDGP)
IF (.NOT. LVORGP) IF_GP = IF_GP - COUNT (YDGP%ITYPE == NTYPE_VOR)
IF (.NOT. LDIVGP) IF_GP = IF_GP - COUNT (YDGP%ITYPE == NTYPE_DIV)

IF_FS = IF_OUT_LT + COUNT ((IAND (YDGP%ITYPE, NTYPE_EW) /= 0) .AND. (YDGP%IVSET == MYSETV))

CALL LTINV1_CTL (IF_OUT_LT, YDGP, YDSP, FSPGL_PROC=FSPGL_PROC)
CALL FTINV1_CTL (IF_GP, IF_FS, IF_OUT_LT, YDGP)

END SUBROUTINE INV_TRANS1_CTL
END MODULE INV_TRANS1_CTL_MOD
