
PROGRAM data_read
  USE data

  IMPLICIT NONE

  CHARACTER(LEN=128) :: data_input_path, data_output_path, buffer
  INTEGER :: begin_year, begin_month, begin_day, begin_hour, &
    begin_minute, begin_second, end_year, end_month, end_day, end_hour, &
    end_minute, end_second
  TYPE(data_file_type) :: data_input_file
  TYPE(data_rows_type) :: data_output
  INTEGER :: data_status
  INTEGER, PARAMETER :: data_unit = 10
  INTEGER :: i

  IF (IARGC() == 14) THEN
    CALL GETARG(1, buffer)
    READ(buffer, '(A)') data_input_path
    CALL GETARG(2, buffer)
    READ(buffer, '(A)') data_output_path
    CALL GETARG(3, buffer)
    READ(buffer,*) begin_year
    CALL GETARG(4, buffer)
    READ(buffer,*) begin_month
    CALL GETARG(5, buffer)
    READ(buffer,*) begin_day
    CALL GETARG(6, buffer)
    READ(buffer,*) begin_hour
    CALL GETARG(7, buffer)
    READ(buffer,*) begin_minute
    CALL GETARG(8, buffer)
    READ(buffer,*) begin_second
    CALL GETARG(9, buffer)
    READ(buffer,*) end_year
    CALL GETARG(10, buffer)
    READ(buffer,*) end_month
    CALL GETARG(11, buffer)
    READ(buffer,*) end_day
    CALL GETARG(12, buffer)
    READ(buffer,*) end_hour
    CALL GETARG(13, buffer)
    READ(buffer,*) end_minute
    CALL GETARG(14, buffer)
    READ(buffer,*) end_second
  ELSE
    STOP
  END IF

  data_input_file%path = data_input_path

  CALL read_data(data_input_file)
  CALL filter_data_by_date(data_input_file%data_struct, data_output, &
    begin_year, begin_month, begin_day, begin_hour, &
    begin_minute, begin_second, end_year, end_month, end_day, end_hour, &
    end_minute, end_second)

  100 FORMAT (1X, T1, 6(I4), 4(F8.2), F7.1, 3(F8.1), F9.6, F8.3, F11.1, &
      F7.1, F10.5)

  OPEN(UNIT=data_unit, FILE=data_output_path, ACTION='WRITE', &
    FORM='FORMATTED', IOSTAT=data_status)
  DO i = 1, data_output%length
    WRITE(data_unit, 100) &
      data_output%rows(i)%year, &
      data_output%rows(i)%month, &
      data_output%rows(i)%day, &
      data_output%rows(i)%hour, &
      data_output%rows(i)%minute, &
      data_output%rows(i)%second, &
      data_output%rows(i)%b, &
      data_output%rows(i)%bx, &
      data_output%rows(i)%by, &
      data_output%rows(i)%bz, &
      data_output%rows(i)%vp, &
      data_output%rows(i)%vx, &
      data_output%rows(i)%vy, &
      data_output%rows(i)%vz, &
      data_output%rows(i)%pth, &
      data_output%rows(i)%np, &
      data_output%rows(i)%tp, &
      data_output%rows(i)%vth, &
      data_output%rows(i)%beta
  END DO
  CLOSE(data_unit)

END PROGRAM data_read

