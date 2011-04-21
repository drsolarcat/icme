
PROGRAM icme

  USE types     ! library of derived data types
  USE fgsl      ! GSL library
  USE config    ! config file routines
  USE data      ! data files routines
  USE datetime  ! datetime routines
  USE mva       ! Maximum variance analysis
  USE gsr       ! Grad-Shafranov analysis
  USE cm        ! Cylinder model

  IMPLICIT NONE

  TYPE(config_file_type) :: config_file

  ! quality of events to be analyzed
  CHARACTER(len=1), PARAMETER :: quality='g'

  TYPE(data_file_type) :: data_file
  TYPE(event_type) :: event

  ! iterators
  INTEGER :: i, j, k, l, m

  !!!!!!!!!!!!!!! START OF THE MAIN PROGRAM !!!!!!!!!!!!!!!

  ! path to config file
  config_file%path='../data/config_f'

  ! read the config file
  CALL read_config(config_file)

  ! start looping through events
  DO i=1, config_file%config_struct%length
    ! select only events flagged with 1 and of needed quality
    IF (config_file%config_struct%rows(i)%flag == 1 .AND. &
      config_file%config_struct%rows(i)%quality == quality) THEN

      ! a piece of code to get the path to the data file
      WRITE (data_file%path, '(I3, A4)') &
        config_file%config_struct%rows(i)%sampling_rate, '.dat'
      data_file%path=ADJUSTL(data_file%path)
      SELECT CASE(TRIM(ADJUSTL(config_file%config_struct%rows(i)%spacecraft)))
        CASE ('WIND')
          data_file%path='../data/wind_'//data_file%path
        CASE ('ACE')
          data_file%path='../data/ace_'//data_file%path
        CASE ('STA')
          data_file%path='../data/stereo_a_'//data_file%path
        CASE ('STB')
          data_file%path='../data/stereo_b_'//data_file%path
        CASE ('MHDA')
          data_file%path='../data/mhd_a_'//data_file%path
        CASE ('MHDB')
          data_file%path='../data/mhd_b_'//data_file%path
      END SELECT
      data_file%path=TRIM(data_file%path)

      ! output the path to the data file
      WRITE(*,*) '---------------------'
      WRITE(*,*) 'Path to the data file'
      WRITE(*,*) data_file%path

      ! read all the data from the data file
      CALL read_data(data_file)

      ! output the length of the data arrays
      WRITE(*,*) '------------------'
      WRITE(*,*) 'Data arrays length'
      WRITE(*,*) data_file%data_struct%length



      ! compute unix timestamps for begin and end of the event
      config_file%config_struct%rows(i)%begin_unixtime = &
        datetime2unixtime(config_file%config_struct%rows(i)%begin_year, &
        config_file%config_struct%rows(i)%begin_month, &
        config_file%config_struct%rows(i)%begin_day, &
        config_file%config_struct%rows(i)%begin_hour, &
        config_file%config_struct%rows(i)%begin_minute, 0)
      config_file%config_struct%rows(i)%end_unixtime = &
        datetime2unixtime(config_file%config_struct%rows(i)%end_year, &
        config_file%config_struct%rows(i)%end_month, &
        config_file%config_struct%rows(i)%end_day, &
        config_file%config_struct%rows(i)%end_hour, &
        config_file%config_struct%rows(i)%end_minute, 0)

      event%config_struct=config_file%config_struct%rows(i)

      ! output event date and spacecraft
      WRITE(*,*) '-------------------------'
      WRITE(*,*) 'Event date and spacecraft'
      WRITE(*,"(1X,I4.4,'-',I2.2,'-',I2.2,' ',A)") event%config_struct%begin_year, &
        event%config_struct%begin_month, event%config_struct%begin_day, &
        event%config_struct%spacecraft

      ! filter data based on the begin/end times
      CALL filter_data_by_unixtime(data_file%data_struct, event%data_struct, &
        event%config_struct%begin_unixtime, event%config_struct%end_unixtime)

!      CALL convert_data_to_si(event%data_struct)

      event%data_struct%sampling_rate=event%config_struct%sampling_rate

      WRITE(*,*) '-------------------'
      WRITE(*,*) 'Event arrays length'
      WRITE(*,*) event%data_struct%length

      ! conduct MVA analysis
      IF (event%config_struct%to_compute_mva == 1) THEN
        CALL mva_compute_mvab(event)
        WRITE(*,*) '---------'
        WRITE(*,*) 'MVAB axes'
        WRITE(*,*) event%mvab%x
        WRITE(*,*) event%mvab%y
        WRITE(*,*) event%mvab%z
        WRITE(*,*) 'criterion: ', event%mvab%criterion
        CALL mva_compute_mvub(event)
        WRITE(*,*) '---------'
        WRITE(*,*) 'MVUB axes'
        WRITE(*,*) event%mvub%x
        WRITE(*,*) event%mvub%y
        WRITE(*,*) event%mvub%z
        WRITE(*,*) 'criterion: ', event%mvub%criterion
      END IF

      IF (event%config_struct%to_compute_gsr == 1) THEN
        CALL compute_gsr()
      END IF

      IF (event%config_struct%to_compute_cylinder_model == 1) THEN
        CALL compute_cylinder_model()
      END IF

      IF (event%config_struct%to_compute_torus_model == 1) THEN
        !CALL compute_torus_model()
      END IF

      IF (event%config_struct%to_compute_hidalgo_model == 1) THEN
        !CALL compute_hidalgo_model()
      END IF

      IF (event%config_struct%to_compute_owens_model == 1) THEN
        !CALL compute_owens_model()
      END IF
    END IF
  END DO
CONTAINS
  SUBROUTINE compute_gsr()
    CALL gsr_compute_dht(event)
    WRITE(*,*) '---------'
    WRITE(*,*) 'dHT speed'
    WRITE(*,*) event%dht%vht
    WRITE(*,*) 'correleation: ', event%dht%cc

    CALL mva_compute_pmvab(event)
    WRITE(*,*) '----------'
    WRITE(*,*) 'PMVAB axes'
    WRITE(*,*) event%pmvab%x
    WRITE(*,*) event%pmvab%y
    WRITE(*,*) event%pmvab%z

    CALL gsr_compute_axes(event)
    !CALL compute_axes(event_bx, event_by, event_bz, event_pth, event_vhtx, &
    !  event_vhty, event_vhtz, config_sampling_rate(i), gsr_x, gsr_y, gsr_z, gsr_theta_opt, &
    !  gsr_phi_opt)
    !CALL compute_gsr(event_bx, event_by, event_bz)
  END SUBROUTINE

  SUBROUTINE compute_cylinder_model()
    CALL cm_fit(event, 240.0)
  END SUBROUTINE compute_cylinder_model
END PROGRAM icme

