
! atomic compare: compare two structures
! find the d, given dmax as a upper limit
! will directly compare coordinates
! for surface/cluster case
SUBROUTINE co_comp(n, x, y, nele, eles, dmax, d, msel)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n, nele
  REAL(8), INTENT(IN) :: x(n, 3), y(n, 3), dmax
  INTEGER, INTENT(IN) :: eles(n) ! must start from 1
  REAL(8), INTENT(OUT) :: d
  INTEGER, INTENT(OUT) :: msel(n) ! output parameter
  
  INTEGER :: nea(nele), ndi(nele, n), ndu(nele, n)
  INTEGER :: sel(n), est(n)
  INTEGER :: i, ai, si, ei
  REAL(8) :: dl(n), dmin, dh
  
  nea = 0 ! how many atoms in each atoms
  ndi = 0 ! atom indices for each atoms
  DO i = 1, n
    nea(eles(i)) = nea(eles(i)) + 1
    ndi(eles(i), nea(eles(i))) = i
  END DO
  ndu = ndi
  ai = 1
  sel(1) = 0
  dmin = dmax
  msel = 0
  DO WHILE (.TRUE.)
    ei = eles(ai) ! ai: cur atom index
    si = 0 ! selected atom
    DO i = 1, nea(ei) ! enum all atoms of this type
      IF (ndu(ei, i) > sel(ai)) THEN ! from the smallest, next time select bigger
        si = ndu(ei, i) ! find the atom index
        est(ai) = i
        ndu(ei, i) = 0 ! cannot be selected
        EXIT
      END IF
    END DO
    IF (si == 0) THEN ! cannot select any
      IF (ai /= 1) THEN
        ai = ai - 1
        ei = eles(ai)
        ndu(ei, est(ai)) = ndi(ei, est(ai))
        CYCLE
      ELSE
        EXIT
      END IF
    END IF
    sel(ai) = si ! sel(ai) cur selected atom index at cur pos
    IF (ai == 1) THEN
      dl(ai) = 0
    ELSE
      dl(ai) = dl(ai - 1)
    END IF
    dh = SQRT(SUM((x(sel(ai), :) - y(ai, :)) ** 2))
    IF (dh > dl(ai)) dl(ai) = dh
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
