
! atomic compare: compare two structures
! find the d, given dmax as a upper limit
! when used in filtering, this will be accurate
! this is very efficient, if dmax is not large
SUBROUTINE at_comp(n, x, y, nele, eles, dmax, d, msel)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n, nele
  REAL(8), INTENT(IN) :: x(n, 3), y(n, 3), dmax
  INTEGER, INTENT(IN) :: eles(n) ! must start from 1
  REAL(8), INTENT(OUT) :: d
  INTEGER, INTENT(OUT) :: msel(n)
  
  INTEGER :: nea(nele), ndi(nele, n), ndu(nele, n)
  INTEGER :: sel(n), est(n)
  INTEGER :: i, j, ai, si, ei
  REAL(8) :: dl(n), dmin, dh, dy(n, n)
  
  nea = 0
  ndi = 0
  DO i = 1, n
    nea(eles(i)) = nea(eles(i)) + 1
    ndi(eles(i), nea(eles(i))) = i
  END DO
  DO i = 1, n
    DO j = 1, i - 1
      dy(j, i) = SQRT(SUM((y(j, :) - y(i, :)) ** 2))
    END DO
  END DO
  ndu = ndi
  ai = 1
  sel(1) = 0
  dmin = dmax
  msel = 0
  DO WHILE (.TRUE.)
    ei = eles(ai)
    si = 0
    DO i = 1, nea(ei)
      IF (ndu(ei, i) > sel(ai)) THEN
        si = ndu(ei, i)
        est(ai) = i
        ndu(ei, i) = 0
        EXIT
      END IF
    END DO
    IF (si == 0) THEN
      IF (ai /= 1) THEN
        ai = ai - 1
        ei = eles(ai)
        ndu(ei, est(ai)) = ndi(ei, est(ai))
        CYCLE
      ELSE
        EXIT
      END IF
    END IF
    sel(ai) = si
    IF (ai == 1) THEN
      dl(ai) = 0
    ELSE
      dl(ai) = dl(ai - 1)
    END IF
    DO i = 1, ai - 1
      dh = ABS(SQRT(SUM((x(sel(ai), :) - x(sel(i), :)) ** 2)) - dy(i, ai))
      IF (dh > dl(ai)) dl(ai) = dh
    END DO
    IF (dl(ai) > dmin .OR. ai == n) THEN
      IF (ai == n .AND. dl(ai) < dmin) THEN
        dmin = dl(ai)
        msel = sel
      END IF
      ndu(ei, est(ai)) = ndi(ei, est(ai))
    ELSE
      ai = ai + 1
      IF (ai <= n) sel(ai) = 0
    END IF
  END DO
  d = dmin
END SUBROUTINE

! return all possible atomic matches below dmax
SUBROUTINE at_comp_list(n, x, y, nele, eles, dmax, nmsel, msel)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n, nele
  REAL(8), INTENT(IN) :: x(n, 3), y(n, 3), dmax
  INTEGER, INTENT(IN) :: eles(n) ! must start from 1
  INTEGER, INTENT(INOUT) :: nmsel
  INTEGER, INTENT(OUT) :: msel(nmsel, n)
  
  INTEGER :: nea(nele), ndi(nele, n), ndu(nele, n)
  INTEGER :: sel(n), est(n)
  INTEGER :: i, j, ai, si, ei, nmsel_max
  REAL(8) :: dl(n), dh, dy(n, n)

  nmsel_max = nmsel
  nmsel = 1
  
  nea = 0
  ndi = 0
  DO i = 1, n
    nea(eles(i)) = nea(eles(i)) + 1
    ndi(eles(i), nea(eles(i))) = i
  END DO
  DO i = 1, n
    DO j = 1, i - 1
      dy(j, i) = SQRT(SUM((y(j, :) - y(i, :)) ** 2))
    END DO
  END DO
  ndu = ndi
  ai = 1
  sel(1) = 0
  msel = 0
  DO WHILE (.TRUE.)
    ei = eles(ai)
    si = 0
    DO i = 1, nea(ei)
      IF (ndu(ei, i) > sel(ai)) THEN
        si = ndu(ei, i)
        est(ai) = i
        ndu(ei, i) = 0
        EXIT
      END IF
    END DO
    IF (si == 0) THEN
      IF (ai /= 1) THEN
        ai = ai - 1
        ei = eles(ai)
        ndu(ei, est(ai)) = ndi(ei, est(ai))
        CYCLE
      ELSE
        EXIT
      END IF
    END IF
    sel(ai) = si
    IF (ai == 1) THEN
      dl(ai) = 0
    ELSE
      dl(ai) = dl(ai - 1)
    END IF
    DO i = 1, ai - 1
      dh = ABS(SQRT(SUM((x(sel(ai), :) - x(sel(i), :)) ** 2)) - dy(i, ai))
      IF (dh > dl(ai)) dl(ai) = dh
    END DO
    IF (dl(ai) > dmax .OR. ai == n) THEN
      IF (ai == n .AND. dl(ai) < dmax) THEN
        IF (nmsel <= nmsel_max) THEN
          msel(nmsel, :) = sel
        END IF
        nmsel = nmsel + 1
      END IF
      ndu(ei, est(ai)) = ndi(ei, est(ai))
    ELSE
      ai = ai + 1
      IF (ai <= n) sel(ai) = 0
    END IF
  END DO
END SUBROUTINE
