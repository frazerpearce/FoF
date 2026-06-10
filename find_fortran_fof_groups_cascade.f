      PROGRAM FOF
C     Timed FoF cascade over nested IC files.
C     Build with: gfortran -std=legacy -O3 -march=native -flto
C
C     Reads files:
C        data/ics_0004096.txt ... data/ics_1048576.txt
C
C     Writes:
C        outdir/fortran_fof_timings.txt
C
C     The timed region is the FoF analysis only.  Position-file reads and
C     timing-file writes are deliberately outside the measured section.

      PARAMETER (NMAX=1048576, NHASH=4194301, MINKEEP=10)
      PARAMETER (NTRIAL=3)
      DOUBLE PRECISION BOX, B
      PARAMETER (BOX=1.0D0, B=0.2D0)

      INTEGER PARENT(NMAX), RANK(NMAX), ROOT(NMAX), GSIZE(NMAX)
      INTEGER GMEMBER(NMAX), ORDER(NMAX), MEMBERS(NMAX)
      INTEGER GROUPID(NMAX), GFIRST(NMAX), WORK(NMAX)
      INTEGER GHEAD(NMAX), GTAIL(NMAX), GNEXT(NMAX)
      INTEGER HKEY(NHASH), HHEAD(NHASH), NEXT(NMAX)
      INTEGER USED(NMAX)
      INTEGER IXCELL(NMAX), IYCELL(NMAX), IZCELL(NMAX)
      DOUBLE PRECISION X(NMAX), Y(NMAX), Z(NMAX)
      DOUBLE PRECISION T0, T1, ELAPSED
      CHARACTER*80 FNAME

      CALL SYSTEM('mkdir -p outdir')

      OPEN(40, FILE='outdir/fortran_fof_timings.txt',
     &     STATUS='UNKNOWN', IOSTAT=IOS)
      IF (IOS .NE. 0) THEN
         PRINT *, 'Could not open outdir/fortran_fof_timings.txt'
         STOP 1
      ENDIF

      WRITE(40, '(A)')
     & '# file n trial seconds raw kept'

      DO 200 IPOW = 12, 20
         NEXPECTED = 2**IPOW
         WRITE(FNAME, 1000) NEXPECTED
 1000    FORMAT('data/ics_', I7.7, '.txt')

         CALL LOADPOSITIONS(FNAME, NMAX, X, Y, Z, N, IOS)
         IF (IOS .NE. 0) THEN
            PRINT *, 'Could not read ', FNAME
            CLOSE(40)
            STOP 1
         ENDIF
         IF (N .NE. NEXPECTED) THEN
            PRINT *, 'Warning: ', FNAME, ' has ', N,
     &               ' particles, expected ', NEXPECTED
         ENDIF

         DO 190 ITRIAL = 1, NTRIAL
            CALL CPU_TIME(T0)
            CALL ANALYSEFOF(N, X, Y, Z, PARENT, RANK, ROOT, GSIZE,
     &           GMEMBER, ORDER, MEMBERS, GROUPID, GFIRST, WORK,
     &           GHEAD, GTAIL, GNEXT, HKEY, HHEAD, NEXT, USED,
     &           IXCELL, IYCELL, IZCELL, NRAW, NKEEP)
            CALL CPU_TIME(T1)
            ELAPSED = T1 - T0

            WRITE(40, 1010) FNAME(1:21), N, ITRIAL, ELAPSED, NRAW, NKEEP
            WRITE(*, 1020) FNAME(1:21), N, ITRIAL, ELAPSED, NRAW, NKEEP
 190     CONTINUE
 200  CONTINUE

 1010 FORMAT(A21,1X,I0,1X,I0,1X,ES16.8E2,1X,I0,1X,I0)
 1020 FORMAT(A21,' N=',I0,' trial=',I0,' time=',ES12.5E2,
     &       ' raw=',I0,' kept=',I0)

      CLOSE(40)
      END

      SUBROUTINE LOADPOSITIONS(FNAME, NMAX, X, Y, Z, N, IOSOUT)
      CHARACTER*(*) FNAME
      INTEGER NMAX, N, IOSOUT
      DOUBLE PRECISION X(*), Y(*), Z(*)

      IOSOUT = 0
      OPEN(10, FILE=FNAME, STATUS='OLD', IOSTAT=IOS)
      IF (IOS .NE. 0) THEN
         IOSOUT = IOS
         RETURN
      ENDIF

      N = 0
 10   CONTINUE
      IF (N .GE. NMAX) GOTO 20
      READ(10, *, END=20, IOSTAT=IOS) X(N+1), Y(N+1), Z(N+1)
      IF (IOS .NE. 0) THEN
         IOSOUT = IOS
         CLOSE(10)
         RETURN
      ENDIF
      N = N + 1
      GOTO 10
 20   CONTINUE
      CLOSE(10)
      RETURN
      END

      SUBROUTINE ANALYSEFOF(N, X, Y, Z, PARENT, RANK, ROOT, GSIZE,
     &     GMEMBER, ORDER, MEMBERS, GROUPID, GFIRST, WORK, GHEAD,
     &     GTAIL, GNEXT, HKEY, HHEAD, NEXT, USED, IXCELL, IYCELL,
     &     IZCELL, NGROUPS, NKEEP)

      PARAMETER (NHASH=4194301, MINKEEP=10)
      DOUBLE PRECISION BOX, B, LL, LL2, DCELL
      PARAMETER (BOX=1.0D0, B=0.2D0)
      INTEGER PARENT(*), RANK(*), ROOT(*), GSIZE(*)
      INTEGER GMEMBER(*), ORDER(*), MEMBERS(*), GROUPID(*), GFIRST(*)
      INTEGER WORK(*), GHEAD(*), GTAIL(*), GNEXT(*)
      INTEGER HKEY(*), HHEAD(*), NEXT(*), USED(*)
      INTEGER IXCELL(*), IYCELL(*), IZCELL(*)
      DOUBLE PRECISION X(*), Y(*), Z(*), NDBLE
      DOUBLE PRECISION DX, DY, DZ, R2

      EXTERNAL FINDROOT
      INTEGER FINDROOT

      IF (N .LE. 0) THEN
         NGROUPS = 0
         NKEEP = 0
         RETURN
      ENDIF

      NDBLE = DBLE(N)
      LL = B * (BOX**3 / NDBLE)**(1.0D0 / 3.0D0)
      LL2 = LL * LL
      NCELL = MAX(1, INT(BOX / LL))
      DCELL = DBLE(NCELL)

      DO 30 I = 1, N
         PARENT(I) = I
         RANK(I) = 0
         ROOT(I) = 0
         GSIZE(I) = 0
         GMEMBER(I) = 0
         GROUPID(I) = 0
         GFIRST(I) = 0
         GHEAD(I) = 0
         GTAIL(I) = 0
         GNEXT(I) = 0
         ORDER(I) = I
 30   CONTINUE

C     Hash occupied mesh cells and chain their particles.
      DO 40 I = 1, NHASH
         HKEY(I) = -1
         HHEAD(I) = 0
 40   CONTINUE
      NUSED = 0
      DO 50 I = 1, N
         IXCELL(I) = MIN(NCELL-1, INT(X(I)*DCELL))
         IYCELL(I) = MIN(NCELL-1, INT(Y(I)*DCELL))
         IZCELL(I) = MIN(NCELL-1, INT(Z(I)*DCELL))
         ICELL = IXCELL(I) + NCELL*(IYCELL(I)+NCELL*IZCELL(I))
         IHASH = 1 + MOD(ICELL, NHASH)
 42      IF (HKEY(IHASH) .EQ. ICELL) GOTO 46
         IF (HKEY(IHASH) .EQ. -1) THEN
            HKEY(IHASH) = ICELL
            NUSED = NUSED + 1
            USED(NUSED) = IHASH
            GOTO 46
         ENDIF
         IHASH = IHASH + 1
         IF (IHASH .GT. NHASH) IHASH = 1
         GOTO 42
 46      NEXT(I) = HHEAD(IHASH)
         HHEAD(IHASH) = I
 50   CONTINUE

C     Search each occupied cell and its 13 forward neighbours.
      DO 59 IU = 1, NUSED
         IHOLD = USED(IU)
         ICELL = HKEY(IHOLD)
         IZ = ICELL / (NCELL*NCELL)
         IREM = ICELL - IZ*NCELL*NCELL
         IY = IREM / NCELL
         IX = IREM - IY*NCELL
         DO 58 KZ = 0, 1
            JZ = MOD(IZ + KZ + NCELL, NCELL)
            DO 57 KY = -1, 1
               IF (KZ .EQ. 0 .AND. KY .LT. 0) GOTO 57
               JY = MOD(IY + KY + NCELL, NCELL)
               DO 56 KX = -1, 1
                  IF (KZ .EQ. 0 .AND. KY .EQ. 0 .AND.
     &                KX .LT. 0) GOTO 56
                  JX = MOD(IX + KX + NCELL, NCELL)
                  JCELL = JX + NCELL*(JY + NCELL*JZ)
                  IHASH = 1 + MOD(JCELL, NHASH)
 51               IF (HKEY(IHASH) .EQ. JCELL) GOTO 47
                  IF (HKEY(IHASH) .EQ. -1) GOTO 56
                  IHASH = IHASH + 1
                  IF (IHASH .GT. NHASH) IHASH = 1
                  GOTO 51
 47               I = HHEAD(IHOLD)
 52               IF (I .EQ. 0) GOTO 56
                  IF (KX .EQ. 0 .AND. KY .EQ. 0 .AND.
     &                KZ .EQ. 0) THEN
                     J = NEXT(I)
                  ELSE
                     J = HHEAD(IHASH)
                  ENDIF
 55               IF (J .EQ. 0) GOTO 54
            DX = DABS(X(I) - X(J))
            IF (DX .GT. 0.5D0 * BOX) DX = BOX - DX
            IF (DX .GT. LL) GOTO 53
            DY = DABS(Y(I) - Y(J))
            IF (DY .GT. 0.5D0 * BOX) DY = BOX - DY
            IF (DY .GT. LL) GOTO 53
            DZ = DABS(Z(I) - Z(J))
            IF (DZ .GT. 0.5D0 * BOX) DZ = BOX - DZ
            IF (DZ .GT. LL) GOTO 53
            R2 = DX*DX + DY*DY + DZ*DZ
            IF (R2 .LE. LL2) CALL UNIONROOT(I, J, PARENT, RANK)
 53               J = NEXT(J)
                  GOTO 55
 54               I = NEXT(I)
                  GOTO 52
 56            CONTINUE
 57         CONTINUE
 58      CONTINUE
 59   CONTINUE

      NGROUPS = 0
      DO 70 I = 1, N
         ROOT(I) = FINDROOT(I, PARENT)
         IGID = GROUPID(ROOT(I))
         IF (IGID .EQ. 0) THEN
            NGROUPS = NGROUPS + 1
            IGID = NGROUPS
            GROUPID(ROOT(I)) = IGID
            GMEMBER(NGROUPS) = ROOT(I)
            GFIRST(NGROUPS) = I
            GHEAD(NGROUPS) = I
            GTAIL(NGROUPS) = I
            ORDER(NGROUPS) = NGROUPS
         ELSE
            GNEXT(GTAIL(IGID)) = I
            GTAIL(IGID) = I
         ENDIF
         GSIZE(IGID) = GSIZE(IGID) + 1
 70   CONTINUE

C     Sort groups by descending size, then ascending first member.
      CALL SORTGROUPS(NGROUPS, ORDER, WORK, GSIZE, GFIRST)

      NKEEP = 0
      DO 110 I = 1, NGROUPS
         IF (GSIZE(ORDER(I)) .GE. MINKEEP) THEN
            NKEEP = NKEEP + 1
         ENDIF
 110  CONTINUE

      RETURN
      END

      INTEGER FUNCTION FINDROOT(I, PARENT)
      INTEGER PARENT(*)
      ICUR = I
 10   CONTINUE
      IF (PARENT(ICUR) .NE. ICUR) THEN
         INEXT = PARENT(ICUR)
         PARENT(ICUR) = PARENT(INEXT)
         ICUR = PARENT(ICUR)
         GOTO 10
      ENDIF
      FINDROOT = ICUR
      RETURN
      END

      SUBROUTINE UNIONROOT(I, J, PARENT, RANK)
      INTEGER PARENT(*), RANK(*)
      EXTERNAL FINDROOT
      INTEGER FINDROOT
      IROOT = FINDROOT(I, PARENT)
      JROOT = FINDROOT(J, PARENT)
      IF (IROOT .EQ. JROOT) RETURN
      IF (RANK(IROOT) .LT. RANK(JROOT)) THEN
         PARENT(IROOT) = JROOT
      ELSEIF (RANK(IROOT) .GT. RANK(JROOT)) THEN
         PARENT(JROOT) = IROOT
      ELSE
         PARENT(JROOT) = IROOT
         RANK(IROOT) = RANK(IROOT) + 1
      ENDIF
      RETURN
      END

      SUBROUTINE SORTGROUPS(NG, ORDER, WORK, GSIZE, GFIRST)
      INTEGER ORDER(*), WORK(*), GSIZE(*), GFIRST(*)
      LOGICAL BETTER
      EXTERNAL BETTER

      IWIDTH = 1
 10   IF (IWIDTH .GE. NG) GOTO 90
      ILEFT = 1
 20   IF (ILEFT .GT. NG) GOTO 70
      IMID = MIN(ILEFT + IWIDTH - 1, NG)
      IRIGHT = MIN(ILEFT + 2*IWIDTH - 1, NG)
      IP = ILEFT
      IQ = IMID + 1
      IR = ILEFT
 30   IF (IP .GT. IMID .OR. IQ .GT. IRIGHT) GOTO 40
      IF (BETTER(ORDER(IP), ORDER(IQ), GSIZE, GFIRST)) THEN
         WORK(IR) = ORDER(IP)
         IP = IP + 1
      ELSE
         WORK(IR) = ORDER(IQ)
         IQ = IQ + 1
      ENDIF
      IR = IR + 1
      GOTO 30
 40   IF (IP .GT. IMID) GOTO 50
      WORK(IR) = ORDER(IP)
      IP = IP + 1
      IR = IR + 1
      GOTO 40
 50   IF (IQ .GT. IRIGHT) GOTO 60
      WORK(IR) = ORDER(IQ)
      IQ = IQ + 1
      IR = IR + 1
      GOTO 50
 60   ILEFT = IRIGHT + 1
      GOTO 20
 70   DO 80 I = 1, NG
         ORDER(I) = WORK(I)
 80   CONTINUE
      IWIDTH = 2 * IWIDTH
      GOTO 10
 90   CONTINUE
      RETURN
      END

      LOGICAL FUNCTION BETTER(IA, IB, GSIZE, GFIRST)
      INTEGER GSIZE(*), GFIRST(*)
      BETTER = .FALSE.
      IF (GSIZE(IA) .GT. GSIZE(IB)) BETTER = .TRUE.
      IF (GSIZE(IA) .EQ. GSIZE(IB) .AND.
     &    GFIRST(IA) .LE. GFIRST(IB)) BETTER = .TRUE.
      RETURN
      END
