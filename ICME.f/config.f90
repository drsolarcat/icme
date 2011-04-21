
MODULE config
  USE types

  IMPLICIT NONE
CONTAINS
  SUBROUTINE read_config(config_file_struct)
    TYPE(config_file_type), INTENT(INOUT) :: config_file_struct
    INTEGER, PARAMETER :: tmp_length=1000, config_unit=10
    TYPE(config_row_type), DIMENSION(tmp_length) :: tmp_config_rows
    INTEGER :: config_status, i

    100 FORMAT (1X, T1, I1, A3, 8I3, A6, 2(2X, I4, 1X, I2, 1X, I2, 1X, I2, 1X, &
      I2), I5, 2I4, 2I5, I4, F6.2, F8.3, F7.3, I3)

    OPEN(UNIT=config_unit, FILE=config_file_struct%path, STATUS='OLD', &
      ACTION='READ', FORM='FORMATTED', IOSTAT=config_status)
    IF (config_status == 0) THEN
      READ(config_unit, *)
      i = 1
      DO
        READ(config_unit, 100, IOSTAT=config_status) &
          tmp_config_rows(i)%flag, &
          tmp_config_rows(i)%quality, &
          tmp_config_rows(i)%to_plot, &
          tmp_config_rows(i)%to_compute_gsr, &
          tmp_config_rows(i)%to_compute_mva, &
          tmp_config_rows(i)%to_compute_cylinder_model, &
          tmp_config_rows(i)%to_compute_torus_model, &
          tmp_config_rows(i)%to_compute_hidalgo_model, &
          tmp_config_rows(i)%to_compute_owens_model, &
          tmp_config_rows(i)%to_save, &
          tmp_config_rows(i)%spacecraft, &
          tmp_config_rows(i)%begin_year, &
          tmp_config_rows(i)%begin_month, &
          tmp_config_rows(i)%begin_day, &
          tmp_config_rows(i)%begin_hour, &
          tmp_config_rows(i)%begin_minute, &
          tmp_config_rows(i)%end_year, &
          tmp_config_rows(i)%end_month, &
          tmp_config_rows(i)%end_day, &
          tmp_config_rows(i)%end_hour, &
          tmp_config_rows(i)%end_minute, &
          tmp_config_rows(i)%sampling_rate, &
          tmp_config_rows(i)%min_theta, &
          tmp_config_rows(i)%max_theta, &
          tmp_config_rows(i)%min_phi, &
          tmp_config_rows(i)%max_phi, &
          tmp_config_rows(i)%num_xpoints, &
          tmp_config_rows(i)%step_ratio, &
          tmp_config_rows(i)%min_y, &
          tmp_config_rows(i)%max_y, &
          tmp_config_rows(i)%order
        IF (config_status == -1 .OR. i+1 > tmp_length) EXIT
        IF (config_status > 0) CYCLE
        i=i+1
      END DO
      CLOSE(config_unit)
      config_file_struct%config_struct%length=i-1
      ALLOCATE(config_file_struct%config_struct%rows( &
        config_file_struct%config_struct%length))
      config_file_struct%config_struct%rows= &
        tmp_config_rows(1:config_file_struct%config_struct%length)
    END IF
  END SUBROUTINE read_config
END MODULE config

