
MODULE cm
  USE types
  USE constants
  USE levenberg_marquardt
  USE resample

  IMPLICIT NONE

  TYPE(cm_rows_type), SAVE :: real_data
CONTAINS
  FUNCTION cm_model(u0, b0, r0_tmp, theta_a, phi_a, p, s, n, dt) RESULT(cm_struct)
    !
    !  Purpose:
    !    To calculate model parameters of the flux rope. This is a model only,
    !    which has to be fitted with real data.
    !

    ! INPUT
    ! bulk flow velocity of the solar wind, or the speed of the MC at the center
    REAL :: u0 ! km/s
    ! the intensity of the magnetic field at the cylinder axis
    REAL :: b0 ! nT
    ! the radius of the MC cylinder at time t=0
    REAL :: r0_tmp
    REAL :: r0 ! AU
    ! the latitude and longitude angles of the cylinder axis
    REAL :: theta_a, phi_a ! degrees
    ! the impact parameter
    REAL :: p ! fracture of r0 [-1,1]
    ! the sign of the magnetic chirality of the MC
    INTEGER :: s ! -1 or 1
    ! number of points to model
    INTEGER :: n
    ! sampling rate of modelled measurements
    REAL :: dt ! s

    ! OUTPUT
    TYPE(cm_rows_type) :: cm_struct

    ! LOCAL
    ! the expansion rate
    REAL :: e ! 1/s
    REAL :: sin_phi, cos_phi, sin_theta, cos_theta
    REAL, DIMENSION(n) :: t, x1, y1, z1, alpha, r, ro, v_ro, b_phi, b_zeta, &
      bx1, by1, bz1, vx1, vy1, vz1, sin_beta, cos_beta
    INTEGER :: i

    ! convert r0 from AU to km
    r0=r0_tmp*ASTRONOMICAL_UNIT/1000.0 ! km

    ALLOCATE(cm_struct%rows(n))

    sin_phi=SQRT(SIN(theta_a*DEG2RAD)**2+COS(theta_a*DEG2RAD)**2*SIN(phi_a*DEG2RAD)**2)
    cos_phi=COS(theta_a*DEG2RAD)*COS(phi_a*DEG2RAD)
    sin_theta=SIN(theta_a*DEG2RAD)/sin_phi
    cos_theta=COS(theta_a*DEG2RAD)*SIN(phi_a*DEG2RAD)/sin_phi

    ! construct time vector
    DO i=1,n
      t(i)=(i-1)*dt ! s
    END DO

    ! new coordinate system, x1-y1 plane contains the flux rope axis
    x1=-r0*SQRT(1-p**2)/sin_phi+u0*t ! km
    y1=0 ! km
    z1=r0*p ! km

    ! calculate expansion rate
    e=(SQRT((x1(1)+u0*t(n))**2*sin_phi**2+y1(1)**2+z1(1)**2)/r0-1)/t(n) ! 1/s

    ! the radial distance from the cylinder axis to the spacecraft
    ro=SQRT(x1**2*sin_phi**2+y1**2+z1**2) ! km

    ! elevation angle of the spacecraft position from the x1-y1 plain
    sin_beta=-x1*sin_phi/ro
    cos_beta=z1/ro

    ! flux rope radius
    r=r0*(1+e*t) ! km

    ! expansion speed
    v_ro=e*ro/(1+e*t) ! km/s

    alpha=2.405/r

    ! compute magnetic field components
    DO i=1,n
!      b_phi(i)=s*b0*BESJ1(alpha(i)*ro(i))/(1+e*t(i))**2 ! nT
!      b_zeta(i)=b0*BESJ0(alpha(i)*ro(i))/(1+e*t(i))**2 ! nT
      b_phi(i)=s*b0*BESSEL_J1(alpha(i)*ro(i))/(1+e*t(i))**2 ! nT
      b_zeta(i)=b0*BESSEL_J0(alpha(i)*ro(i))/(1+e*t(i))**2 ! nT
    END DO

    ! expansion speed in new coordinates
    vx1=-v_ro*cos_beta*sin_phi ! km/s
    vy1=v_ro*cos_beta*cos_phi ! km/s
    vz1=v_ro*sin_beta ! km/s

    ! magnetic field in new coordinates
    bx1=b_zeta*cos_phi+b_phi*sin_beta*sin_phi ! nT
    by1=b_zeta*sin_phi-b_phi*sin_beta*cos_phi ! nT
    bz1=b_phi*cos_beta ! nT

    ! rotate the axes back into initial coordinates
    !   velocity
    cm_struct%rows(:)%vx=-v_ro*cos_beta*sin_phi-u0 ! km/s
    cm_struct%rows(:)%vy=vy1*cos_theta-vz1*sin_theta ! km/s
    cm_struct%rows(:)%vz=vy1*sin_theta+vz1*cos_theta ! km/s
    !   magnetic field
    cm_struct%rows(:)%bx=bx1 ! nT
    cm_struct%rows(:)%by=by1*cos_theta-bz1*sin_theta ! nT
    cm_struct%rows(:)%bz=by1*sin_theta+bz1*cos_theta ! nT
  END FUNCTION cm_model

  SUBROUTINE cm_f(m, n, p, fvec, iflag)
    INTEGER, INTENT(IN) :: m, n
    REAL, INTENT(IN) :: p(:)
    REAL, INTENT(IN OUT) :: fvec(:)
    INTEGER, INTENT(IN OUT) :: iflag

    TYPE(cm_rows_type) :: model_data

    model_data=cm_model(p(1), p(2), p(3), p(4), p(5), p(6), 1, &
      real_data%length, real_data%sampling_rate)

    fvec = (model_data%rows(:)%bx-real_data%rows(:)%bx)**2+ &
      (model_data%rows(:)%by-real_data%rows(:)%by)**2+ &
      (model_data%rows(:)%bz-real_data%rows(:)%bz)**2+ &
      0.1*(SQRT(model_data%rows(:)%vx**2+model_data%rows(:)%vy**2+model_data%rows(:)%vz**2)- &
      SQRT(real_data%rows(:)%vx**2+real_data%rows(:)%vy**2+real_data%rows(:)%vz**2))**2
  END SUBROUTINE cm_f

  SUBROUTINE cm_fit(event, sampling_rate)
    TYPE(event_type) :: event
    REAL :: sampling_rate

    REAL :: object=0.0, object_final=1.0E7, simp, stopcr, rnd
    REAL, DIMENSION(6) :: params_lower, params_upper, params, params_final
    INTEGER, DIMENSION(6) :: iwa
    INTEGER :: info
    REAL, DIMENSION(:), ALLOCATABLE :: fvec
    INTEGER :: i, k

    real_data%sampling_rate=sampling_rate
    real_data%length=(event%data_struct%length-1)*event%data_struct%sampling_rate/sampling_rate+1

    ALLOCATE(real_data%rows(real_data%length), fvec(real_data%length))

    real_data%rows(:)%bx=resample_by_timestep(event%data_struct%rows(:)%bx, &
      event%data_struct%sampling_rate, sampling_rate)
    real_data%rows(:)%by=resample_by_timestep(event%data_struct%rows(:)%by, &
      event%data_struct%sampling_rate, sampling_rate)
    real_data%rows(:)%bz=resample_by_timestep(event%data_struct%rows(:)%bz, &
      event%data_struct%sampling_rate, sampling_rate)
    real_data%rows(:)%vx=resample_by_timestep(event%data_struct%rows(:)%vx, &
      event%data_struct%sampling_rate, sampling_rate)
    real_data%rows(:)%vy=resample_by_timestep(event%data_struct%rows(:)%vy, &
      event%data_struct%sampling_rate, sampling_rate)
    real_data%rows(:)%vz=resample_by_timestep(event%data_struct%rows(:)%vz, &
      event%data_struct%sampling_rate, sampling_rate)

    params_lower = (/ ABS(SUM(real_data%rows(:)%vx)/real_data%length)-50, &
      0.1, 0.01, 0.0, 0.0, -1.0 /)
    params_upper = (/ ABS(SUM(real_data%rows(:)%vx)/real_data%length)+50, &
      30.0, 0.2, 180.0, 360.0, 1.0 /)
    params = (/ ABS(SUM(real_data%rows(:)%vx)/real_data%length), &
      15.0, 0.01, 45.0, 45.0, 0.1 /)

    CALL RANDOM_SEED()
    DO i = 1, 1000
      DO k = 1, 6
        CALL RANDOM_NUMBER(rnd)
        params(k) = params_lower(k)+(params_upper(k)-params_lower(k))*rnd
      END DO
      CALL lmdif1(cm_f, real_data%length, 6, params, fvec, 1.0E-04, info, iwa)
      object = SUM(fvec)/real_data%length
      IF (ISNAN(object)) CYCLE
      IF (object < object_final .AND. &
        params(1) <= params_upper(1) .AND. params(1) >= params_lower(1) .AND. &
        params(2) <= params_upper(2) .AND. params(2) >= params_lower(2) .AND. &
!        params(3) <= params_upper(3) .AND. params(3) >= params_lower(3) .AND. &
!        params(4) <= params_upper(4) .AND. params(4) >= params_lower(4) .AND. &
!        params(5) <= params_upper(5) .AND. params(5) >= params_lower(5) .AND. &
        params(6) <= params_upper(6) .AND. params(6) >= params_lower(6)) THEN
        WRITE(*,*) object
        object_final = object
        params_final = params
      END IF
    END DO
    object = object_final
    params = params_final

!    DO
!      IF (params(4) >= 0.0 .AND. params(4) <= 90.0) EXIT
!      params(4) = params(4)
!    END DO

!    MODULO(params(4), 360)

    WRITE(*, 900) object, params
    900 FORMAT(' Success !'/' Objective function = ', f12.6/ ' at: ', 6f12.6)
!    WRITE(*, 910) var
!    910 FORMAT(' Elements of var = ', 6f12.6)
  END SUBROUTINE cm_fit
END MODULE cm

