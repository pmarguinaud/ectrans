! (C) Copyright 1995- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE TRLTOG1_MOD

USE PARKIND1, ONLY : JPIM, JPRB

IMPLICIT NONE

PUBLIC TRLTOG1
PRIVATE TRLTOG1_PROLOG, TRLTOG1_COMM, TRLTOG1_CTX


TYPE TRLTOG1_CTX
  INTEGER(KIND=JPIM) :: ISENDCOUNT = -9999
  INTEGER(KIND=JPIM) :: IRECVCOUNT = -9999
  INTEGER(KIND=JPIM) :: INSEND = -9999
  INTEGER(KIND=JPIM) :: INRECV = -9999
  INTEGER(KIND=JPIM), ALLOCATABLE :: ISENDTOT (:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: IRECVTOT (:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: ISEND    (:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: IRECV    (:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: IINDEX(:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: INDOFF(:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: IGPTRSEND(:,:,:)
  INTEGER(KIND=JPIM), ALLOCATABLE :: ISETAL(:), ISETBL(:), ISETWL(:), ISETVL(:)
END TYPE

CONTAINS

SUBROUTINE TRLTOG1(PGLAT,KF_FS,KF_GP,YDGP)

!**** *TRLTOG1 * - head routine for transposition of grid point data from latitudinal
!                 to column structure (this takes place between inverse
!                 FFT and grid point calculations)
!                 TRLTOG1 is the inverse of TRGTOL

!**   Interface.
!     ----------
!        *call* *TRLTOG1(...)

!        Explicit arguments :
!        --------------------
!           PGLAT    -  Latitudinal data ready for direct FFT (input)
!           PGP    -  Blocked grid point data    (output)
!           KVSET    - "v-set" for each field      (input)

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        R. El Khatib *Meteo-France*

!     Modifications.
!     --------------
!        Original  : 18-Aug-2014 from trltog
!        R. El Khatib 09-Sep-2020 NSTACK_MEMORY_TR
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_GEN         ,ONLY : NSTACK_MEMORY_TR
USE TPM_DISTR       ,ONLY : D, NPRTRNS, NPROC
USE TPM_TRANS       ,ONLY : NGPBLKS
USE POINTER_MOD

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KF_FS,KF_GP
REAL(KIND=JPRB),   INTENT(IN)  :: PGLAT(KF_FS,D%NLENGTF)
TYPE(PTRG)                     :: YDGP (:)

TYPE (TRLTOG1_CTX) :: YLCTX

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TRLTOG1',0,ZHOOK_HANDLE)

CALL TRLTOG1_PROLOG(KF_FS,KF_GP,YLCTX,YDGP%IVSET)

CALL TRLTOG1_COMM(PGLAT,KF_FS,KF_GP,YLCTX,YDGP)

IF (LHOOK) CALL DR_HOOK('TRLTOG1',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE TRLTOG1

SUBROUTINE TRLTOG1_PROLOG(KF_FS,KF_GP,YDCTX,KVSET)

!**** *TRLTOG1_PROLOG * - prolog for transposition of grid point data from latitudinal
!                 to column structure (this takes place between inverse
!                 FFT and grid point calculations) : the purpose is essentially 
!                 to compute the size of communication buffers in order to enable
!                 the use of automatic arrays later.
!                 TRLTOG1_PROLOG is the inverse of TRGTOL_PROLOG

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *call* *TRLTOG1_PROLOG(...)

!        Explicit arguments :
!        --------------------
!           KVSET    - "v-set" for each field      (input)

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        R. El Khatib *Meteo-France*

!     Modifications.
!     --------------
!        Original  : 18-Aug-2014 from trltog
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DISTR       ,ONLY : D, MYSETW, NPRTRNS, MYPROC, NPROC
USE TPM_TRANS       ,ONLY : NGPBLKS

USE INIGPTR_MOD     ,ONLY : INIGPTR
USE PE2SET_MOD      ,ONLY : PE2SET
USE TPM_GEN         ,ONLY : NSTACK_MEMORY_TR
!

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)    :: KF_FS,KF_GP
TYPE (TRLTOG1_CTX), INTENT(INOUT) :: YDCTX
INTEGER(KIND=JPIM), INTENT(IN)    :: KVSET(:)

INTEGER(KIND=JPIM) :: IGPTRRECV(NPRTRNS)
INTEGER(KIND=JPIM) :: IFIRSTLAT, IGL, IGLL, ILASTLAT, IPOS, ISETA, ISETB, ISETV
INTEGER(KIND=JPIM) :: ISEND, JFLD, JGL, JL, ISETW, JROC, J
INTEGER(KIND=JPIM) :: INDOFFX,IBUFLENS,IBUFLENR

!     ------------------------------------------------------------------

ALLOCATE (YDCTX%ISENDTOT (NPROC))
ALLOCATE (YDCTX%IRECVTOT (NPROC))
ALLOCATE (YDCTX%ISEND    (NPROC))
ALLOCATE (YDCTX%IRECV    (NPROC))
ALLOCATE (YDCTX%IINDEX(D%NLENGTF))
ALLOCATE (YDCTX%INDOFF(NPROC))
ALLOCATE (YDCTX%IGPTRSEND(2,NGPBLKS,NPRTRNS))
ALLOCATE (YDCTX%ISETAL(NPROC), YDCTX%ISETBL(NPROC), YDCTX%ISETWL(NPROC), YDCTX%ISETVL(NPROC))

!*       0.    Some initializations
!              --------------------

CALL GSTATS(1806,0)

CALL INIGPTR(YDCTX%IGPTRSEND,IGPTRRECV)

INDOFFX  = 0
IBUFLENS = 0
IBUFLENR = 0
YDCTX%INRECV = 0
YDCTX%INSEND = 0

YDCTX%INDOFF = -999

DO JROC=1,NPROC

  CALL PE2SET(JROC,YDCTX%ISETAL(JROC),YDCTX%ISETBL(JROC),YDCTX%ISETWL(JROC),YDCTX%ISETVL(JROC))
  ISEND      = JROC
  ISETA=YDCTX%ISETAL(JROC)
  ISETB=YDCTX%ISETBL(JROC)
  ISETW=YDCTX%ISETWL(JROC)
  ISETV=YDCTX%ISETVL(JROC)
!             count up expected number of fields
  IPOS = 0
  DO JFLD=1,KF_GP
    IF(KVSET(JFLD) == ISETV .OR. KVSET(JFLD) == -1) IPOS = IPOS+1
  ENDDO
  YDCTX%IRECVTOT(JROC) = IGPTRRECV(ISETW)*IPOS
  IF(YDCTX%IRECVTOT(JROC) > 0 .AND. MYPROC /= JROC) THEN
    YDCTX%INRECV = YDCTX%INRECV + 1
    YDCTX%IRECV(YDCTX%INRECV)=JROC
  ENDIF

  IF( JROC /= MYPROC) IBUFLENR = MAX(IBUFLENR,YDCTX%IRECVTOT(JROC))

  IFIRSTLAT = MAX(D%NPTRLS(MYSETW),D%NFRSTLAT(ISETA))
  ILASTLAT  = MIN(D%NPTRLS(MYSETW)+D%NULTPP(MYSETW)-1,D%NLSTLAT(ISETA))

  IPOS = 0
  DO JGL=IFIRSTLAT,ILASTLAT
    IGL  = D%NPTRFRSTLAT(ISETA)+JGL-D%NFRSTLAT(ISETA)
    IPOS = IPOS+D%NONL(IGL,ISETB)
  ENDDO

  YDCTX%ISENDTOT(JROC) = IPOS*KF_FS
  IF( JROC /= MYPROC) THEN
    IBUFLENS = MAX(IBUFLENS,YDCTX%ISENDTOT(JROC))
    IF(YDCTX%ISENDTOT(JROC) > 0) THEN
      YDCTX%INSEND = YDCTX%INSEND+1
      YDCTX%ISEND(YDCTX%INSEND)=JROC
    ENDIF
  ENDIF

  IF(IPOS > 0) THEN
    YDCTX%INDOFF(JROC) = INDOFFX
    INDOFFX = INDOFFX+IPOS
    IPOS = 0
    DO JGL=IFIRSTLAT,ILASTLAT
      IGL  = D%NPTRFRSTLAT(ISETA)+JGL-D%NFRSTLAT(ISETA)
      IGLL = JGL-D%NPTRLS(MYSETW)+1
      DO JL=D%NSTA(IGL,ISETB)+D%NSTAGTF(IGLL),&
       &D%NSTA(IGL,ISETB)+D%NSTAGTF(IGLL)+D%NONL(IGL,ISETB)-1
        IPOS = IPOS+1
        YDCTX%IINDEX(IPOS+YDCTX%INDOFF(JROC)) = JL
      ENDDO
    ENDDO
  ENDIF
ENDDO

YDCTX%ISENDCOUNT=0
YDCTX%IRECVCOUNT=0
DO J=1,NPROC
  YDCTX%ISENDCOUNT=MAX(YDCTX%ISENDCOUNT,YDCTX%ISENDTOT(J))
  YDCTX%IRECVCOUNT=MAX(YDCTX%IRECVCOUNT,YDCTX%IRECVTOT(J))
ENDDO

CALL GSTATS(1806,1)

END SUBROUTINE TRLTOG1_PROLOG

SUBROUTINE TRLTOG1_COMM(PGLAT,KF_FS,KF_GP,YDCTX,YDGP)
 

!**** *trltog * - transposition of grid point data from latitudinal
!                 to column structure. This takes place between inverse
!                 FFT and grid point calculations.
!                 TRLTOG1_COMM is the inverse of TRGTOL

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *call* *trltog(...)

!        Explicit arguments :
!        --------------------
!           PGLAT    -  Latitudinal data ready for direct FFT (input)
!           PGP    -  Blocked grid point data    (output)
!           KVSET    - "v-set" for each field      (input)

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        MPP Group *ECMWF*

!     Modifications.
!     --------------
!        Original  : 95-10-01
!        D.Dent    : 97-08-04 Reorganisation to allow NPRTRV
!                             to differ from NPRGPEW
!        =99-03-29= Mats Hamrud and Deborah Salmond
!                   JUMP in FFT's changed to 1
!                   KINDEX introduced and PCOMBUF not used for same PE
!         01-11-23  Deborah Salmond and John Hague
!                   LIMP_NOOLAP Option for non-overlapping message passing
!                               and buffer packing
!         01-12-18  Peter Towers
!                   Improved vector performance of LTOG_PACK,LTOG_UNPACK
!         03-0-02   G. Radnoti: Call barrier always when nproc>1
!         08-01-01  G.Mozdzynski: cleanup
!         09-01-02  G.Mozdzynski: use non-blocking recv and send
!        R. El Khatib 09-Sep-2020 64 bits addressing for PGLAT
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB    ,JPIA
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE MPL_MODULE  ,ONLY : MPL_RECV, MPL_SEND, MPL_WAIT, JP_NON_BLOCKING_STANDARD, MPL_WAITANY, &
     & JP_BLOCKING_STANDARD, MPL_BARRIER, JP_BLOCKING_BUFFERED

USE TPM_GEN         ,ONLY : NOUT, NTRANS_SYNC_LEVEL, NSTACK_MEMORY_TR
USE TPM_DISTR       ,ONLY : D, MYSETV, MYSETW, MTAGLG,      &
     &                      NPRCIDS, NPRTRNS, MYPROC, NPROC
USE TPM_TRANS       ,ONLY : LDIVGP, LSCDERS, LUVDER, LVORGP, NGPBLKS

USE PE2SET_MOD      ,ONLY : PE2SET
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
USE POINTER_MOD

IMPLICIT NONE


INTEGER(KIND=JPIM), INTENT(IN)            :: KF_FS, KF_GP
TYPE (TRLTOG1_CTX), INTENT(INOUT), TARGET :: YDCTX
REAL(KIND=JPRB),    INTENT(IN)            :: PGLAT(KF_FS,D%NLENGTF)
TYPE(PTRG)                                :: YDGP (:)

! LOCAL VARIABLES

INTEGER(KIND=JPIM) :: IPOSPLUS(YDCTX%INRECV)
INTEGER(KIND=JPIM) :: ISETW(YDCTX%INRECV)
INTEGER(KIND=JPIM) :: JPOS(NGPBLKS,YDCTX%INRECV)
INTEGER(KIND=JPIM) :: IFLDA(KF_GP,YDCTX%INRECV)
INTEGER(KIND=JPIM) :: IREQ_SEND(NPROC)
INTEGER(KIND=JPIM) :: IREQ_RECV(NPROC)

INTEGER(KIND=JPIM) :: IFIRST, IFLD, ILAST, IPOS, ISETA, ISETB, IRECV, ISETV
INTEGER(KIND=JPIM) :: ISEND, ITAG,  JBLK, JFLD, JK, JL, IFLDS, IPROC,JROC, INR, INS
INTEGER(KIND=JPIM) :: II,ILEN,IBUFLENS,IBUFLENR, IFLDT, JI, JJ, J

INTEGER(KIND=JPIA) :: JFLD64

INTEGER(KIND=JPIM) :: JNR
INTEGER(KIND=JPIM) :: IFLDOFF(KF_FS)
INTEGER(KIND=JPIM) :: IRECV_FLD_START,IRECV_FLD_END
INTEGER(KIND=JPIM) :: ISEND_FLD_START,ISEND_FLD_END
INTEGER(KIND=JPIM) :: IGPTROFF(NGPBLKS)

REAL(KIND=JPRB), TARGET :: ZCOMBUFS0(-1:YDCTX%ISENDCOUNT,MERGE (YDCTX%INSEND,0,NSTACK_MEMORY_TR/=0))
REAL(KIND=JPRB), TARGET :: ZCOMBUFR0(-1:YDCTX%IRECVCOUNT,MERGE (YDCTX%INRECV,0,NSTACK_MEMORY_TR/=0))

REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZCOMBUFS(:,:)
REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZCOMBUFR(:,:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE_BAR

!     ------------------------------------------------------------------

!*       0.    Some initializations
!              --------------------

ASSOCIATE (KSENDCOUNT=>YDCTX%ISENDCOUNT, KRECVCOUNT=>YDCTX%IRECVCOUNT, KNSEND=>YDCTX%INSEND, KNRECV=>YDCTX%INRECV, &
& KSENDTOT=>YDCTX%ISENDTOT, KRECVTOT=>YDCTX%IRECVTOT, KSEND=>YDCTX%ISEND, KRECV=>YDCTX%IRECV, KINDEX=>YDCTX%IINDEX, &
& KNDOFF=>YDCTX%INDOFF, KGPTRSEND =>YDCTX%IGPTRSEND, KSETAL=>YDCTX%ISETAL, KSETBL=>YDCTX%ISETBL, KSETWL=>YDCTX%ISETWL, &
& KSETVL=>YDCTX%ISETVL)

IF (NSTACK_MEMORY_TR == 0) THEN
  ALLOCATE(ZCOMBUFS(-1:YDCTX%ISENDCOUNT,YDCTX%INSEND))
  ALLOCATE(ZCOMBUFR(-1:YDCTX%IRECVCOUNT,YDCTX%INRECV))
! Now, force the OS to allocate this shared array right now, not when it starts to be used which is
! an OPEN-MP loop, that would cause a threads synchronization lock :
  IF (YDCTX%INSEND > 0 .AND. YDCTX%ISENDCOUNT >=-1) ZCOMBUFS(-1,1)=HUGE(1._JPRB)
ELSE
  ZCOMBUFS (-1:,1:) => ZCOMBUFS0
  ZCOMBUFR (-1:,1:) => ZCOMBUFR0
ENDIF

ITAG   = MTAGLG

IF (LHOOK) CALL DR_HOOK('TRLTOG1_BAR',0,ZHOOK_HANDLE_BAR)
CALL GSTATS_BARRIER(762)
IF (LHOOK) CALL DR_HOOK('TRLTOG1_BAR',1,ZHOOK_HANDLE_BAR)

CALL GSTATS(805,0)

IF (NTRANS_SYNC_LEVEL <= 0) THEN
   !...Receive loop.........................................................
   DO INR=1,KNRECV
      IRECV=KRECV(INR)
      CALL MPL_RECV(ZCOMBUFR(-1:KRECVTOT(IRECV),INR), &
           & KSOURCE=NPRCIDS(IRECV), &
           & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQ_RECV(INR), &
           & KTAG=ITAG,CDSTRING='TRLTOG1_COMM: NON-BLOCKING IRECV' )
   ENDDO
ENDIF

CALL GSTATS(805,1)

! Copy local contribution
IF( KRECVTOT(MYPROC) > 0 )THEN
  IFLDS = 0
  DO JFLD=1,KF_GP
    IF(YDGP(JFLD)%IVSET == MYSETV .OR. YDGP(JFLD)%IVSET == -1) THEN
      IFLDS = IFLDS+1
      IFLDOFF(IFLDS) = JFLD
    ENDIF
  ENDDO

  IPOS=0
  DO JBLK=1,NGPBLKS
    IGPTROFF(JBLK)=IPOS
    IFIRST = KGPTRSEND(1,JBLK,MYSETW)
    IF(IFIRST > 0) THEN
      ILAST = KGPTRSEND(2,JBLK,MYSETW)
      IPOS=IPOS+ILAST-IFIRST+1
    ENDIF
  ENDDO

  CALL GSTATS(1604,0)
#ifdef __NEC__
! Loops inversion is still better on Aurora machines, according to CHMI. REK.
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(JFLD64,JBLK,JK,IFLD,IPOS,IFIRST,ILAST)
#else
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JFLD64,JBLK,JK,IFLD,IPOS,IFIRST,ILAST)
#endif

  DO JBLK=1,NGPBLKS
    IFIRST = KGPTRSEND(1,JBLK,MYSETW)
    IF(IFIRST > 0) THEN
      ILAST = KGPTRSEND(2,JBLK,MYSETW)
! Address PGLAT over 64 bits because its size may exceed 2 GB for big data and
! small number of tasks.
        DO JFLD64 = 1, IFLDS
          IFLD = IFLDOFF(JFLD64)
          DO JK=IFIRST,ILAST
            IPOS = KNDOFF(MYPROC)+IGPTROFF(JBLK)+JK-IFIRST+1
            YDGP (IFLD)%P (JK, JBLK) = PGLAT(JFLD64,KINDEX(IPOS))
          ENDDO
        ENDDO
    ENDIF
  ENDDO
!$OMP END PARALLEL DO
  CALL GSTATS(1604,1)

ENDIF

!
! loop over the number of processors we need to communicate with.
! NOT MYPROC
!
! Now overlapping buffer packing/unpacking with sends/waits
! Time as if all communications to avoid double accounting

CALL GSTATS(805,0)

!  Pack+send loop.........................................................

ISEND_FLD_START = 1
ISEND_FLD_END   = KF_FS
DO INS=1,KNSEND
  ISEND=KSEND(INS)
  ILEN = KSENDTOT(ISEND)/KF_FS
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JFLD,JL,II)
  DO JL=1,ILEN
    II = KINDEX(KNDOFF(ISEND)+JL)
    DO JFLD=ISEND_FLD_START,ISEND_FLD_END
      ZCOMBUFS((JFLD-ISEND_FLD_START)*ILEN+JL,INS) = PGLAT(JFLD,II)
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
  ZCOMBUFS(-1,INS) = 1
  ZCOMBUFS(0,INS)  = KF_FS
  IF (NTRANS_SYNC_LEVEL <= 1) THEN
     CALL MPL_SEND(ZCOMBUFS(-1:KSENDTOT(ISEND),INS),KDEST=NPRCIDS(ISEND),&
          & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQ_SEND(INS), &
          & KTAG=ITAG,CDSTRING='TRLTOG1_COMM: NON-BLOCKING ISEND')
  ELSE
     CALL MPL_SEND(ZCOMBUFS(-1:KSENDTOT(ISEND),INS),KDEST=NPRCIDS(ISEND),&
          & KMP_TYPE=JP_BLOCKING_BUFFERED, &
          & KTAG=ITAG,CDSTRING='TRLTOG1_COMM: BLOCKING BUFFERED BSEND')
  ENDIF
ENDDO

!  Unpack loop.........................................................

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(INR,IRECV,ISETA,ISETB,ISETV,IFLD,JFLD,IPOS,JBLK,IFIRST,ILAST)
DO INR=1,KNRECV
  IRECV=KRECV(INR)
  
  ISETA=KSETAL(IRECV)
  ISETB=KSETBL(IRECV)
  ISETW(INR)=KSETWL(IRECV)
  ISETV=KSETVL(IRECV)
  IFLD = 0
  DO JFLD=1,KF_GP
    IF(YDGP(JFLD)%IVSET == ISETV .OR. YDGP(JFLD)%IVSET == -1 ) THEN
      IFLD = IFLD+1
      IFLDA(IFLD,INR)=JFLD
    ENDIF
  ENDDO
  IPOS = 0
  IPOSPLUS(INR)=0
  DO JBLK=1,NGPBLKS
    IFIRST = KGPTRSEND(1,JBLK,ISETW(INR))
    IF(IFIRST > 0) THEN
      ILAST = KGPTRSEND(2,JBLK,ISETW(INR))
      JPOS(JBLK,INR)=IPOS
      IPOSPLUS(INR)=IPOSPLUS(INR)+(ILAST-IFIRST+1)
      IPOS=IPOS+(ILAST-IFIRST+1)
    ENDIF
  ENDDO
ENDDO
!$OMP END PARALLEL DO

DO JNR=1,KNRECV
  
  IF (NTRANS_SYNC_LEVEL <= 0) THEN
     CALL MPL_WAITANY(KREQUEST=IREQ_RECV(1:KNRECV),KINDEX=INR,&
          & CDSTRING='TRLTOG1_COMM: WAIT FOR ANY RECEIVES')
  ELSE
     INR = JNR
     IRECV=KRECV(INR)
     CALL MPL_RECV(ZCOMBUFR(-1:KRECVTOT(IRECV),INR), &
          & KSOURCE=NPRCIDS(IRECV), &
          & KMP_TYPE=JP_BLOCKING_STANDARD, &
          & KTAG=ITAG,CDSTRING='TRLTOG1_COMM: BLOCKING RECV' )
  ENDIF

  IPOS=IPOSPLUS(INR)
  IRECV_FLD_START = ZCOMBUFR(-1,INR)
  IRECV_FLD_END   = ZCOMBUFR(0,INR)
  
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IFLDT,IFIRST,ILAST,JK,JJ,JI,JBLK)
  DO JJ=IRECV_FLD_START,IRECV_FLD_END
    IFLDT=IFLDA(JJ,INR)
    DO JBLK=1,NGPBLKS
      IFIRST = KGPTRSEND(1,JBLK,ISETW(INR))
      IF(IFIRST > 0) THEN
        ILAST = KGPTRSEND(2,JBLK,ISETW(INR))
        DO JK=IFIRST,ILAST
          JI=(JJ-IRECV_FLD_START)*IPOS+JPOS(JBLK,INR)+JK-IFIRST+1
          YDGP (IFLDT)%P (JK,JBLK) = ZCOMBUFR(JI,INR)
        ENDDO
      ENDIF
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
ENDDO

IF (NTRANS_SYNC_LEVEL <= 1) THEN
   IF(KNSEND > 0) THEN
      CALL MPL_WAIT(KREQUEST=IREQ_SEND(1:KNSEND),CDSTRING='TRLTOG1_COMM: WAIT FOR ISENDS')
   ENDIF
ENDIF

IF (NTRANS_SYNC_LEVEL >= 1) THEN
   CALL MPL_BARRIER(CDSTRING='TRLTOG1_COMM: BARRIER AT END')
ENDIF

CALL GSTATS(805,1)

CALL GSTATS_BARRIER2(762)

IF (NSTACK_MEMORY_TR == 0) THEN
  DEALLOCATE(ZCOMBUFR)
  DEALLOCATE(ZCOMBUFS)
ENDIF

END ASSOCIATE

END SUBROUTINE TRLTOG1_COMM

END MODULE TRLTOG1_MOD
