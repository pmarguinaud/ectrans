SUBROUTINE UNPAIR(N, A)
! unpair.pl: one block per case
INTEGER :: N, J
REAL :: A(N)
! pair, ACC first: DATA moves to the line of TARGET DATA
#ifdef ACCGPU
!$ACC DATA PRESENT(A) &
!$ACC& ASYNC(1)
#endif
#ifdef OMPGPU
!$OMP TARGET DATA &
!$OMP& MAP(PRESENT,ALLOC:A)
#endif
! pair, OMP first: the kernel moves to the OMP line, comment of the ACC half kept
#ifdef OMPGPU
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SHARED(A,N) MAP(TO:N)
#endif
#ifdef ACCGPU
! acc comment
!$ACC PARALLEL LOOP FIRSTPRIVATE(N) PRESENT(A) ASYNC(1)
#endif
DO J=1,N
  A(J) = 0.
ENDDO
! acc-only
#ifdef ACCGPU
!$ACC WAIT(1)
#endif
! omp-only, optional END
#ifdef OMPGPU
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
#endif
! pair, END DATA + WAIT vs END TARGET DATA
#ifdef OMPGPU
!$OMP END TARGET DATA
#endif
#ifdef ACCGPU
!$ACC WAIT(1)
!$ACC END DATA
#endif
! code pair: not removable
#ifdef ACCGPU
CALL ACC_SPECIFIC()
#endif
#ifdef OMPGPU
CALL OMP_SPECIFIC()
#endif
! OMP directive without counterpart: not removable
#ifdef ACCGPU
!$ACC UPDATE HOST(A)
#endif
#ifdef OMPGPU
!$OMP TARGET UPDATE FROM(A)
!$OMP TARGET EXIT DATA MAP(DELETE:A)
#endif
END SUBROUTINE
