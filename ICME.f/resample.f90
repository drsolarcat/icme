
MODULE resample
  USE cubic_spline_gcv

  IMPLICIT NONE
CONTAINS
  FUNCTION resample_by_length(series_in, n_in, n_out) RESULT(series_out)
    INTEGER :: n_in, n_out
    REAL, DIMENSION(n_in) :: series_in
    REAL, DIMENSION(n_out) :: series_out

    REAL :: dt_in, dt_out

    dt_in=1.0
    dt_out = (n_in-1)*dt_in/(n_out-1)

    CALL resample_generic(series_in, series_out, n_in, n_out, dt_in, dt_out)
  END FUNCTION resample_by_length

  FUNCTION resample_by_timestep(series_in, dt_in, dt_out) RESULT(series_out)
    REAL :: dt_in, dt_out
    REAL, DIMENSION(:) :: series_in
    REAL, DIMENSION(:), ALLOCATABLE :: series_out

    INTEGER :: n_in, n_out

    n_in = SIZE(series_in)
    n_out = (n_in-1)*dt_in/dt_out+1

    ALLOCATE(series_out(n_out))

    CALL resample_generic(series_in, series_out, n_in, n_out, dt_in, dt_out)
  END FUNCTION resample_by_timestep

  SUBROUTINE resample_generic(series_in, series_out, n_in, n_out, &
    dt_in, dt_out)
    REAL :: dt_in, dt_out
    INTEGER :: n_in, n_out
    REAL, DIMENSION(n_in) :: series_in
    REAL, DIMENSION(n_out) :: series_out

    REAL, DIMENSION(n_in) :: t_in
    INTEGER :: i, k
    REAL :: d
    INTEGER :: job, ier
    REAL(dp) :: y(n_in), c(n_in-1, 3), df(n_in), var, se(n_in), wk(7*(n_in+2))

    df = 1.0D0
    var = -1.0D0
    job = 1

    DO i = 1, n_in
      t_in(i) = (i-1)*dt_in
    END DO

    CALL cubgcv(t_in, series_in, df, n_in, y, c, n_in-1, var, job, se, wk, ier)

    DO i = 1, n_out
      k = MINLOC(ABS(t_in-i*dt_out), DIM=1)
      IF (k*dt_in > i*dt_out) k = k-1
      d = i*dt_out - k*dt_in
      series_out(i) = ((c(k,3)*d+c(k,2))*d+c(k,1))*d+y(k)
    END DO
  END SUBROUTINE resample_generic
END MODULE resample

