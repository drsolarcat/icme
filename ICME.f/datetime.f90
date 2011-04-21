
MODULE datetime
  IMPLICIT NONE
CONTAINS
  FUNCTION datetime2unixtime(year, month, day, hour, minute, second) &
    RESULT(unixtime)
    INTEGER, INTENT(IN) :: year, month, day, hour, minute, second
    INTEGER :: unixtime, leap_years, doy

    doy=dayofmonth2dayofyear(year, month, day)
    leap_years=FLOOR((year-1-1972)/4.0)+1

    unixtime=second+minute*60+hour*3600+&
      (doy+leap_years*366+(year-1971)*365)*86400

  END FUNCTION datetime2unixtime

  FUNCTION dayofmonth2dayofyear(year, month, day) RESULT(day_of_year)
    INTEGER, INTENT(IN) :: year, month, day
    INTEGER :: day_of_year
    INTEGER :: m, is_leap

    is_leap=0

    IF (MOD(year, 4) == 0 .AND. .NOT. (MOD(year, 100) == 0 .AND. &
      MOD(year, 400) /= 0)) THEN
      is_leap=1
    END IF

    day_of_year=day

    DO m=1, 12
      IF (m == month) EXIT
      SELECT CASE(m)
        CASE (1, 3, 5, 7, 8, 10, 12)
          day_of_year=day_of_year+31
        CASE (4, 6, 9, 11)
          day_of_year=day_of_year+30
        CASE (2)
          day_of_year=day_of_year+28+is_leap
      END SELECT
    END DO
  END FUNCTION dayofmonth2dayofyear
END MODULE datetime

