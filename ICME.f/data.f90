
MODULE data
!
!  Purpose:
!    To read data from files and apply to it some necessary operations like
!    filtering and conversion to/from different unit systems.
!
!  Author:
!    Alexey Isavnin
!
  USE types
  USE datetime

  IMPLICIT NONE
CONTAINS

  SUBROUTINE filter_data_by_unixtime(data_input, data_output, &
    begin_unixtime, end_unixtime)
  !
  !  Purpose:
  !    To filter data by unix timestamps.
  !
    TYPE(data_rows_type), INTENT(INOUT) :: data_input
    TYPE(data_rows_type), OPTIONAL, INTENT(OUT) :: data_output
    INTEGER, INTENT(IN) :: begin_unixtime, end_unixtime
    INTEGER :: unixtime
    INTEGER :: i
    INTEGER, PARAMETER :: tmp_length=5000
    TYPE(data_row_type), DIMENSION(tmp_length) :: tmp_data_rows
    TYPE(data_rows_type) :: tmp_data

    DO i=1, data_input%length
      unixtime = datetime2unixtime(&
        data_input%rows(i)%year, data_input%rows(i)%month, &
        data_input%rows(i)%day, data_input%rows(i)%hour, &
        data_input%rows(i)%minute, data_input%rows(i)%second)

      IF (i == 1 .AND. unixtime > begin_unixtime) EXIT
      IF (unixtime >= begin_unixtime .AND. unixtime <= end_unixtime) THEN
        tmp_data%length=tmp_data%length+1
        tmp_data_rows(tmp_data%length)=data_input%rows(i)
        CYCLE
      END IF
      IF (unixtime > end_unixtime) EXIT
    END DO
    ALLOCATE(tmp_data%rows(tmp_data%length))
    tmp_data%rows=tmp_data_rows(1:tmp_data%length)
    IF (PRESENT(data_output)) THEN
      data_output=tmp_data
    ELSE
      data_input=tmp_data
    END IF
  END SUBROUTINE filter_data_by_unixtime

  SUBROUTINE filter_data_by_date(data_input, data_output, begin_year, &
    begin_month, begin_day, begin_hour, begin_minute, begin_second, &
    end_year, end_month, end_day, end_hour, end_minute, end_second)
    TYPE(data_rows_type), INTENT(INOUT) :: data_input
    TYPE(data_rows_type), INTENT(OUT) :: data_output
    INTEGER, INTENT(IN) :: begin_year, begin_month, begin_day, begin_hour, &
      begin_minute, begin_second, end_year, end_month, end_day, end_hour, &
      end_minute, end_second

    CALL filter_data_by_unixtime(data_input, data_output, &
      datetime2unixtime(begin_year, begin_month, begin_day, begin_hour, &
      begin_minute, begin_second), &
      datetime2unixtime(end_year, end_month, end_day, end_hour, &
      end_minute, end_second))
  END SUBROUTINE filter_data_by_date

  SUBROUTINE read_data(data_file_struct)
  !
  !  Purpose:
  !    To read data from ascii file.
  !
    TYPE(data_file_type), INTENT(INOUT) :: data_file_struct
    INTEGER, PARAMETER :: tmp_length=2000000, data_unit=20
    TYPE(data_row_type), DIMENSION(tmp_length) :: tmp_data_rows
    INTEGER :: data_status, i

    200 FORMAT (1X, T1, 6(I4), 4(F8.2), F7.1, 3(F8.1), F9.6, F8.3, F11.1, &
      F7.1, F10.5)

    OPEN(UNIT=data_unit, FILE=data_file_struct%path, STATUS='OLD', &
      ACTION='READ', FORM='FORMATTED', IOSTAT=data_status)
    IF (data_status == 0) THEN
      READ(data_unit, *)
      READ(data_unit, *)
      i=1
      DO
        READ(data_unit, 200, IOSTAT=data_status) &
          tmp_data_rows(i)%year, &
          tmp_data_rows(i)%month, &
          tmp_data_rows(i)%day, &
          tmp_data_rows(i)%hour, &
          tmp_data_rows(i)%minute, &
          tmp_data_rows(i)%second, &
          tmp_data_rows(i)%b, &
          tmp_data_rows(i)%bx, &
          tmp_data_rows(i)%by, &
          tmp_data_rows(i)%bz, &
          tmp_data_rows(i)%vp, &
          tmp_data_rows(i)%vx, &
          tmp_data_rows(i)%vy, &
          tmp_data_rows(i)%vz, &
          tmp_data_rows(i)%pth, &
          tmp_data_rows(i)%np, &
          tmp_data_rows(i)%tp, &
          tmp_data_rows(i)%vth, &
          tmp_data_rows(i)%beta
        IF (data_status == -1) EXIT
        IF (data_status > 0) CYCLE
        i=i+1
      END DO
      CLOSE(data_unit)
      data_file_struct%data_struct%length = i-1
      ALLOCATE(data_file_struct%data_struct%rows( &
        data_file_struct%data_struct%length))
      data_file_struct%data_struct%rows = &
        tmp_data_rows(1:data_file_struct%data_struct%length)
    END IF
  END SUBROUTINE read_data

  SUBROUTINE convert_data_to_si(data_struct)
  !
  !  Purpose:
  !    To convert data from default units used in the data files to SI units.
  !
    TYPE(data_rows_type), INTENT(INOUT) :: data_struct
    data_struct%rows(:)%b=data_struct%rows(:)%b*1e-9     ! nT -> T
    data_struct%rows(:)%bx=data_struct%rows(:)%bx*1e-9   ! nT -> T
    data_struct%rows(:)%by=data_struct%rows(:)%by*1e-9   ! nT -> T
    data_struct%rows(:)%bz=data_struct%rows(:)%bz*1e-9   ! nT -> T
    data_struct%rows(:)%vp=data_struct%rows(:)%vp*1e3    ! km/s -> m/s
    data_struct%rows(:)%vx=data_struct%rows(:)%vx*1e3    ! km/s -> m/s
    data_struct%rows(:)%vy=data_struct%rows(:)%vy*1e3    ! km/s -> m/s
    data_struct%rows(:)%vz=data_struct%rows(:)%vz*1e3    ! km/s -> m/s
    data_struct%rows(:)%pth=data_struct%rows(:)%pth*1e-9 ! nPa -> Pa
    data_struct%rows(:)%np=data_struct%rows(:)%np*1e6    ! cm^-3 -> m^-3
    data_struct%rows(:)%vth=data_struct%rows(:)%vth*1e3  ! km/s -> m/s
  END SUBROUTINE convert_data_to_si
END MODULE data

