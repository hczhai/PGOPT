
! Find the difference and atomic matching of
! two structures, by direct comparing coordinates
! using Kuhn-Munkres Algorithm
! the metric is the mean of all difference

! km algorithm with custom weights
SUBROUTINE kmc_comp(n, ww, nele, eles, d, msel)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, nele
    REAL(8), INTENT(IN) :: ww(n, n)
    INTEGER, INTENT(IN) :: eles(n) ! must start from 1
    REAL(8), INTENT(OUT) :: d
    INTEGER, INTENT(OUT) :: msel(n) ! output parameter

    REAL(8), ALLOCATABLE :: w(:, :)
    INTEGER, ALLOCATABLE :: wi(:), wx(:)
    INTEGER :: i, k, m
    DO k = 1, nele
        m = COUNT(eles == k)
        ALLOCATE ( wi(m), w(m, m), wx(m) )
        m = 0
        DO i = 1, n
            IF (eles(i) /= k) CYCLE
            m = m + 1
            wi(m) = i
        END DO
        w = ww(wi, wi)
        CALL kuhn_munkres(m, w, wx)
        msel(wi) = wi(wx)
        DEALLOCATE ( wi, w, wx )
    END DO
    d = 0d0
    DO i = 1, n
        d = d + ww(msel(i), i) / n
    END DO
END SUBROUTINE

SUBROUTINE weights(n, x, y, nele, eles, cell, ww)
    IMPLICIT NONE
    REAL(8), PARAMETER :: eps = 1D-10
    INTEGER, INTENT(IN) :: n, nele
    REAL(8), INTENT(IN) :: x(n, 3), y(n, 3), cell(3, 3)
    INTEGER, INTENT(IN) :: eles(n) ! must start from 1
    REAL(8), INTENT(OUT) :: ww(n, n)

    INTEGER :: i, j, k, m
    REAL(8) :: d
    DO i = 1, n
        DO j = 1, n
            IF (eles(i) /= eles(j)) CYCLE
            ww(i, j) = SQRT(SUM((x(i, :) - y(j, :)) ** 2))
            IF (ABS(cell(1, 1)) > eps) THEN
                DO k = -1, 1
                    DO m = -1, 1
                        IF (k == 0 .AND. m == 0) CYCLE
                        d = SQRT(SUM((x(i, :) - y(j, :) - &
                            k * cell(1, :) - m * cell(2, :)) ** 2))
                        IF (d < ww(i, j)) ww(i, j) = d
                    END DO
                END DO
            END IF
        END DO
    END DO
END SUBROUTINE

SUBROUTINE km_comp(n, x, y, nele, eles, cell, d, msel)
    IMPLICIT NONE
    REAL(8), PARAMETER :: eps = 1D-10
    INTEGER, INTENT(IN) :: n, nele
    REAL(8), INTENT(IN) :: x(n, 3), y(n, 3), cell(3, 3)
    INTEGER, INTENT(IN) :: eles(n) ! must start from 1
    REAL(8), INTENT(OUT) :: d
    INTEGER, INTENT(OUT) :: msel(n) ! output parameter

    REAL(8), ALLOCATABLE :: w(:, :)
    INTEGER, ALLOCATABLE :: wi(:), wx(:)
    REAL(8) :: ww(n, n)
    INTEGER :: i, j, k, m
    DO i = 1, n
        DO j = 1, n
            IF (eles(i) /= eles(j)) CYCLE
            ww(i, j) = SQRT(SUM((x(i, :) - y(j, :)) ** 2))
            IF (ABS(cell(1, 1)) > eps) THEN
                DO k = -1, 1
                    DO m = -1, 1
                        IF (k == 0 .AND. m == 0) CYCLE
                        d = SQRT(SUM((x(i, :) - y(j, :) - &
                            k * cell(1, :) - m * cell(2, :)) ** 2))
                        IF (d < ww(i, j)) ww(i, j) = d
                    END DO
                END DO
            END IF
        END DO
    END DO
    DO k = 1, nele
        m = COUNT(eles == k)
        ALLOCATE ( wi(m), w(m, m), wx(m) )
        m = 0
        DO i = 1, n
            IF (eles(i) /= k) CYCLE
            m = m + 1
            wi(m) = i
        END DO
        w = ww(wi, wi)
        CALL kuhn_munkres(m, w, wx)
        msel(wi) = wi(wx)
        DEALLOCATE ( wi, w, wx )
    END DO
    d = 0d0
    DO i = 1, n
        d = d + ww(msel(i), i) / n
    END DO
END SUBROUTINE

! km algorithm with custom weights, output best k solutions
SUBROUTINE kmc_comp_best(n, k, ww, nele, eles, d, msel)
    IMPLICIT NONE
    REAL(8), PARAMETER :: inf = 1D5, eps = 1D-10
    INTEGER, INTENT(IN) :: n, nele
    REAL(8), INTENT(IN) :: ww(n, n)
    INTEGER, INTENT(IN) :: eles(n), k ! must start from 1
    REAL(8), INTENT(OUT) :: d(k)
    INTEGER, INTENT(OUT) :: msel(n, k) ! output parameter

    REAL(8) :: w(n, n, k), wc(n, n)
    REAL(8) :: xd
    INTEGER :: p, mp, u, pi, v, i, xmsel(n)

    w(:, :, 1) = ww
    CALL kmc_comp(n, w(:, :, 1), nele, eles, d(1), msel(:, 1))
    mp = 1
    DO p = 1, k - 1
        DO u = 1, n
            wc = w(:, :, p)
            wc(msel(u, p), u) = inf
            CALL kmc_comp(n, wc, nele, eles, xd, xmsel)
            v = mp + 1
            pi = 0
            DO v = p + 1, mp
                IF (d(v) > xd) EXIT
                IF (ABS(d(v) - xd) < eps .AND. ALL(xmsel == msel(:, v))) THEN
                    pi = 1
                    EXIT
                END IF
            END DO
            IF (pi == 1) CYCLE
            DO i = MIN(k, mp + 1), v + 1, -1
                w(:, :, i) = w(:, :, i - 1)
                d(i) = d(i - 1)
                msel(:, i) = msel(:, i - 1)
            END DO
            IF (v /= k + 1) THEN
                w(:, :, v) = wc
                d(v) = xd
                msel(:, v) = xmsel
                IF (mp < k) mp = mp + 1
            END IF
        END DO
    END DO
END SUBROUTINE

RECURSIVE FUNCTION match(n, s, t, lx, ly, slack, w, x, u) RESULT(l)
    IMPLICIT NONE
    REAL(8), PARAMETER :: eps = 1E-10
    INTEGER, INTENT(IN) :: n, u
    REAL(8), INTENT(IN) :: lx(n), ly(n), w(n, n)
    REAL(8), INTENT(INOUT) :: slack(n)
    LOGICAL, INTENT(INOUT) :: s(n), t(n)
    INTEGER, INTENT(INOUT) :: x(n)
    LOGICAL :: l
    INTEGER :: v
    REAL(8) :: k
    s(u) = .TRUE.
    DO v = 1, n
        IF (t(v)) CYCLE
        k = lx(u) + ly(v) + w(u, v)
        IF (ABS(k) < eps) THEN
            t(v) = .TRUE.
            IF (x(v) == 0 .OR. match(n, s, t, lx, ly, slack, w, x, x(v))) THEN
                x(v) = u
                l = .TRUE.
                RETURN
            END IF
        ELSE
            IF (k < slack(v)) slack(v) = k
        END IF
    END DO
    l = .FALSE.
END FUNCTION

! https://blog.sengxian.com/algorithms/km
! Kuhn-Munkres Algorithm O(n^3)
SUBROUTINE kuhn_munkres(n, w, x)
    IMPLICIT NONE
    INTERFACE
        RECURSIVE FUNCTION match(n, s, t, lx, ly, slack, w, x, u) RESULT(l)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: n, u
            REAL(8), INTENT(IN) :: lx(n), ly(n), w(n, n)
            REAL(8), INTENT(INOUT) :: slack(n)
            LOGICAL, INTENT(INOUT) :: s(n), t(n)
            INTEGER, INTENT(INOUT) :: x(n)
            LOGICAL :: l
        END FUNCTION
    END INTERFACE
    REAL(8), PARAMETER :: inf = 1E10
    INTEGER, INTENT(IN) :: n
    REAL(8), INTENT(IN) :: w(n, n)
    INTEGER, INTENT(OUT) :: x(n)

    REAL(8) :: lx(n), ly(n), slack(n), d
    LOGICAL :: s(n), t(n)
    INTEGER :: i, j

    DO i = 1, n
        lx(i) = -MINVAL(w(i, :))
    END DO
    ly = 0
    x = 0
    DO i = 1, n
        slack = inf
        DO
            s = .FALSE.
            t = .FALSE.
            IF (match(n, s, t, lx, ly, slack, w, x, i)) THEN
                EXIT
            ELSE
                d = MINVAL(slack, .NOT. t)
                DO j = 1, n
                    IF (s(j)) lx(j) = lx(j) - d
                    IF (t(j)) ly(j) = ly(j) + d
                END DO
            END IF
        END DO
    END DO
END SUBROUTINE