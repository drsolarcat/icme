
MODULE gsr
  USE fgsl
  USE extblas
  USE types
  USE mva
  USE constants

  IMPLICIT NONE
CONTAINS
  SUBROUTINE gsr_compute_dht(event)
    TYPE(event_type), INTENT(INOUT) :: event
    INTEGER :: n, i
    REAL, DIMENSION(event%data_struct%length) :: vx, vy, vz, bx, by, bz
    REAL, DIMENSION(3*event%data_struct%length) :: c1, c2
    REAL :: vhtx, vhty, vhtz

    n=event%data_struct%length
    vx=event%data_struct%rows(:)%vx
    vy=event%data_struct%rows(:)%vy
    vz=event%data_struct%rows(:)%vz
    bx=event%data_struct%rows(:)%bx
    by=event%data_struct%rows(:)%by
    bz=event%data_struct%rows(:)%bz

    CALL gsr_loop_dht( &
      MINVAL(vx)-10.0, 10.0, MAXVAL(vx)+10.0, &
      MINVAL(vy)-10.0, 10.0, MAXVAL(vy)+10.0, &
      MINVAL(vz)-10.0, 10.0, MAXVAL(vz)+10.0, &
      vx, vy, vz, bx, by, bz, n, vhtx, vhty, vhtz)
    CALL gsr_loop_dht(vhtx-10.0, 1.0, vhtx+10.0, vhty-10.0, 1.0, vhty+10.0, &
      vhtz-10.0, 1.0, vhtz+10.0, vx, vy, vz, bx, by, bz, n, vhtx, vhty, vhtz)
    CALL gsr_loop_dht(vhtx-1.0, 0.1, vhtx+1.0, vhty-1.0, 0.1, vhty+1.0, &
      vhtz-1.0, 0.1, vhtz+1.0, vx, vy, vz, bx, by, bz, n, vhtx, vhty, vhtz)

    event%dht%vht(1)=vhtx
    event%dht%vht(2)=vhty
    event%dht%vht(3)=vhtz

    CALL extblas_mvcross_expand( &
      -event%data_struct%rows(:)%vx, &
      -event%data_struct%rows(:)%vy, &
      -event%data_struct%rows(:)%vz, &
      event%data_struct%rows(:)%bx, &
      event%data_struct%rows(:)%by, &
      event%data_struct%rows(:)%bz, &
      c1(1:n), c1(n+1:2*n), c1(2*n+1:3*n))
    CALL extblas_mvcross_expand( &
      -event%dht%vht(1)*extblas_ones_vector(event%data_struct%length), &
      -event%dht%vht(2)*extblas_ones_vector(event%data_struct%length), &
      -event%dht%vht(3)*extblas_ones_vector(event%data_struct%length), &
      event%data_struct%rows(:)%bx, &
      event%data_struct%rows(:)%by, &
      event%data_struct%rows(:)%bz, &
      c2(1:n), c2(n+1:2*n), c2(2*n+1:3*n))

    event%dht%cc=ABS(fgsl_stats_correlation(c1, 1, c2, 1, 3*n))
  END SUBROUTINE gsr_compute_dht

  SUBROUTINE gsr_loop_dht(begin_vfx, delta_vfx, end_vfx, begin_vfy, delta_vfy, &
    end_vfy, begin_vfz, delta_vfz, end_vfz, vx, vy, vz, bx, by, bz, &
    n, vhtx, vhty, vhtz)
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: begin_vfx, delta_vfx, end_vfx, &
      begin_vfy, delta_vfy, end_vfy, begin_vfz, delta_vfz, end_vfz
    REAL, DIMENSION(n), INTENT(IN) :: vx, vy, vz, bx, by, bz
    REAL, INTENT(OUT) :: vhtx, vhty, vhtz
    REAL :: vfx, vfy, vfz, d, dmin
    REAL, DIMENSION(n) :: ex, ey, ez

    dmin=0.0

    vfx=begin_vfx
    DO WHILE (vfx <= end_vfx)
      vfy=begin_vfy
      DO WHILE (vfy <= end_vfy)
        vfz=begin_vfz
        DO WHILE (vfz <= end_vfz)
          CALL extblas_mvcross_expand( &
            vx-vfx*extblas_ones_vector(n), vy-vfy*extblas_ones_vector(n), &
            vz-vfz*extblas_ones_vector(n), bx, by, bz, ex, ey, ez)
          d = SUM(ex**2+ey**2+ez**2)/n
          IF (dmin == 0.0 .OR. d < dmin) THEN
            dmin = d
            vhtx = vfx
            vhty = vfy
            vhtz = vfz
          END IF
          vfz=vfz+delta_vfz
        END DO
        vfy=vfy+delta_vfy
      END DO
      vfx=vfx+delta_vfx
    END DO
  END SUBROUTINE gsr_loop_dht

  SUBROUTINE gsr_compute_axes(event)
    TYPE(event_type), INTENT(INOUT) :: event
    INTEGER :: file_status,i
    INTEGER, DIMENSION(2) :: max_angles
    REAL, DIMENSION(3) :: x,y,z
    TYPE(curve_type) :: curve

    CALL mva_compute_pmvab(event)
    ALLOCATE(event%residue(10,37),event%branch_length(10,37))
    CALL gsr_optimize_axis(event, 0, 10, 90, 0, 10, 360)

    max_angles=MAXLOC(event%residue**(-1)*event%branch_length)
    WRITE(*,*) max_angles

    z=extblas_rotate_vector(event%pmvab%z, event%pmvab%y, max_angles(1)*PI/180)
    z=extblas_rotate_vector(z, event%pmvab%z, max_angles(2)*PI/180)
    z=z/extblas_norm(z)
    x=-(event%dht%vht-DOT_PRODUCT(event%dht%vht, z)*z)
    x=x/extblas_norm(x)
    y=extblas_vcross(z, x)
    y=y/extblas_norm(y)

    curve=gsr_compute_curve(event, x, y, z)

    OPEN(UNIT=50, FILE='./curve.dat', &
      ACTION='WRITE', FORM='FORMATTED', IOSTAT=file_status)
    DO i=1,curve%length
      WRITE(50,'(1X,F15.5,F15.5)') curve%x(i), curve%y(i)
    END DO

    CLOSE(50)

    OPEN(UNIT=30, FILE='./residue.dat', &
      ACTION='WRITE', FORM='FORMATTED', IOSTAT=file_status)
    OPEN(UNIT=40, FILE='./length.dat', &
      ACTION='WRITE', FORM='FORMATTED', IOSTAT=file_status)
    DO i=1,10
      WRITE(30,'(1X,37F15.5)') event%residue(i,:)
      WRITE(40,'(1X,37I15)') event%branch_length(i,:)
    END DO

    CLOSE(30)
    CLOSE(40)


!    IF (event%config_struct%min_theta /= 0 .OR. &
!      event%config_struct%max_theta /=0 .OR. &
!      event%config_struct%min_phi /= 0 .OR. &
!      event%config_struct%max_phi /= 0) THEN
!      CALL gsr_optimize_axis(event, &
!        event%config_struct%min_theta, 1, event%config_struct%max_theta, &
!        event%config_struct%min_phi, 1, event%config_struct%max_phi)
!    END IF
!    CALL gsr_optimize_axis(event, &
!      event%gsr%theta-5, 1, event%gsr%theta+5, &
!      event%gsr%phi-5, 1, event%gsr%phi+5)
  END SUBROUTINE gsr_compute_axes

  SUBROUTINE gsr_optimize_axis(event, theta_min, theta_delta, theta_max, &
    phi_min, phi_delta, phi_max)
    TYPE(event_type), INTENT(INOUT) :: event
    INTEGER, INTENT(IN) :: theta_min, theta_delta, theta_max, &
      phi_min, phi_delta, phi_max
    INTEGER :: theta, phi, i,k
    REAL, DIMENSION(3) :: x, y, z
    TYPE(curve_type) :: curve

    i=1
    theta=theta_min
    DO WHILE (theta <= theta_max)
      z=extblas_rotate_vector(event%pmvab%z, event%pmvab%y, theta*PI/180)
      k=1
      phi=phi_min
      DO WHILE (phi <= phi_max)
        z=extblas_rotate_vector(z, event%pmvab%z, phi*PI/180)
        z=z/extblas_norm(z)
        x=-(event%dht%vht-DOT_PRODUCT(event%dht%vht, z)*z)
        x=x/extblas_norm(x)
        y=extblas_vcross(z, x)
        y=y/extblas_norm(y)
        curve=gsr_compute_curve(event, x, y, z)
        CALL gsr_retreive_branches_simple(curve)
        CALL gsr_compute_residue(curve)
        event%residue(i,k)=curve%residue
        event%branch_length(i,k)=curve%branch_length
        phi=phi+phi_delta
        k=k+1
      END DO
      theta=theta+theta_delta
      i=i+1
    END DO
  END SUBROUTINE gsr_optimize_axis

  FUNCTION gsr_compute_curve(event, x, y, z) RESULT(curve)
    TYPE(event_type), INTENT(IN) :: event
    REAL, DIMENSION(3), INTENT(IN) :: x, y, z
    REAL, DIMENSION(event%data_struct%length) :: bx,by,bz
    REAL :: dx
    INTEGER :: i
    TYPE(fgsl_function) :: f_obj
    TYPE(curve_type) :: curve

    CALL extblas_mvprojection( &
      event%data_struct%rows(:)%bx, &
      event%data_struct%rows(:)%by, &
      event%data_struct%rows(:)%bz, x, y, z, bx, by, bz)

    dx=-DOT_PRODUCT(event%dht%vht, x)*event%config_struct%sampling_rate;

    ALLOCATE(curve%x(event%data_struct%length), &
      curve%y(event%data_struct%length))
    curve%length=event%data_struct%length-1

    DO i=1,curve%length
      curve%x(i)=simp(i+1,dx,by(1:i+1))
      curve%y(i)=(event%data_struct%rows(i+1)%pth+bz(i+1)**2/2/ &
        VACUUM_PERMEABILITY)*1e9
    END DO
  END FUNCTION gsr_compute_curve

  SUBROUTINE gsr_compute_residue(curve)
    TYPE(curve_type), INTENT(INOUT) :: curve
    TYPE(fgsl_interp_accel) :: f_accelerator_in,f_accelerator_out
    TYPE(fgsl_spline) :: f_spline_in,f_spline_out
    REAL :: x_left,x_right,x_delta,x
    REAL, DIMENSION(curve%length) :: y_in,y_out
    INTEGER, DIMENSION(curve%branch_in%length) :: p_in
    INTEGER, DIMENSION(curve%branch_out%length) :: p_out
    INTEGER :: f_status,i

    IF (curve%is_branched .AND. curve%branch_in%length > 100 .AND. &
      curve%branch_out%length > 100) THEN
      f_accelerator_in=fgsl_interp_accel_alloc()
      f_accelerator_out=fgsl_interp_accel_alloc()
      f_spline_in=fgsl_spline_alloc(FGSL_INTERP_CSPLINE,curve%branch_in%length)
      f_spline_out=fgsl_spline_alloc(FGSL_INTERP_CSPLINE, &
        curve%branch_out%length)

      CALL fgsl_sort_index(p_in,curve%branch_in%x,1,curve%branch_in%length)
      CALL fgsl_sort_index(p_out,curve%branch_out%x,1,curve%branch_out%length)
      f_status=fgsl_permute(p_in,curve%branch_in%x,1,curve%branch_in%length)
      f_status=fgsl_permute(p_in,curve%branch_in%y,1,curve%branch_in%length)
      f_status=fgsl_permute(p_out,curve%branch_out%x,1,curve%branch_out%length)
      f_status=fgsl_permute(p_out,curve%branch_out%y,1,curve%branch_out%length)

      f_status=fgsl_spline_init(f_spline_in,curve%branch_in%x, &
        curve%branch_in%y,curve%branch_in%length)
      f_status=fgsl_spline_init(f_spline_out,curve%branch_out%x, &
        curve%branch_out%y,curve%branch_out%length)

      x_left=MAX(MINVAL(curve%branch_in%x),MINVAL(curve%branch_out%x))
      x_right=MIN(MAXVAL(curve%branch_in%x),MAXVAL(curve%branch_out%x))
      x_delta=(x_right-x_left)/(curve%length-1)

      DO i=1,curve%length
        x=x_left+(i-1)*x_delta
        y_in(i)=fgsl_spline_eval(f_spline_in,x,f_accelerator_in)
        y_out(i)=fgsl_spline_eval(f_spline_out,x,f_accelerator_out)
      END DO

      curve%residue=SQRT(SUM((y_in-y_out)**2))/ &
        ABS(MIN(MAXVAL(y_in),MAXVAL(y_out))-MAX(MINVAL(y_in),MINVAL(y_out)))
      curve%branch_length=MIN(curve%branch_in%length,curve%branch_out%length)

      CALL fgsl_spline_free(f_spline_in)
      CALL fgsl_spline_free(f_spline_out)
      CALL fgsl_interp_accel_free(f_accelerator_in)
      CALL fgsl_interp_accel_free(f_accelerator_out)
    ELSE
      curve%residue=1e7
      curve%branch_length=0
    END IF
  END SUBROUTINE gsr_compute_residue

  SUBROUTINE gsr_retreive_branches_simple(curve)
    TYPE(curve_type) :: curve
    REAL, DIMENSION(curve%length) :: x,y
    REAL :: min_left_x,min_right_x
    INTEGER :: max_index,min_left_index,min_right_index
    LOGICAL :: is_branched

    x=curve%x
    y=curve%y

    IF (ABS(MINVAL(x)) > ABS(MAXVAL(x))) x=-x

    max_index=MAXLOC(x,1)
    min_left_index=1
    min_right_index=curve%length
    is_branched=.TRUE.

    IF (max_index == 1 .OR. max_index == curve%length) THEN
      is_branched=.FALSE.
    ELSE
      min_left_index=MINLOC(x(1:max_index),1)
      min_right_index=MINLOC(x(max_index:curve%length),1)+max_index-1
      min_left_x=x(min_left_index)
      min_right_x=x(min_right_index)
      IF (min_left_x >= min_right_x) THEN
        min_right_index=MINLOC(ABS(x(max_index:curve%length)-min_left_x),1)+ &
          max_index-1
      ELSE
        min_left_index=MINLOC(ABS(x(1:max_index)-min_right_x),1)
      END IF
    END IF

    curve%is_branched=is_branched

    curve%branch_in%left_index=min_left_index
    curve%branch_in%right_index=max_index
    curve%branch_out%left_index=max_index
    curve%branch_out%right_index=min_right_index

    curve%branch_in%length=max_index-min_left_index+1
    curve%branch_out%length=min_right_index-max_index+1
    ALLOCATE( &
      curve%branch_in%x(curve%branch_in%length), &
      curve%branch_in%y(curve%branch_in%length), &
      curve%branch_out%x(curve%branch_out%length), &
      curve%branch_out%y(curve%branch_out%length))

    curve%branch_in%x= &
      curve%x(curve%branch_in%left_index:curve%branch_in%right_index)
    curve%branch_in%y= &
      curve%y(curve%branch_in%left_index:curve%branch_in%right_index)
    curve%branch_out%x= &
      curve%x(curve%branch_out%left_index:curve%branch_out%right_index)
    curve%branch_out%y= &
      curve%y(curve%branch_out%left_index:curve%branch_out%right_index)
  END SUBROUTINE gsr_retreive_branches_simple

!  SUBROUTINE gsr_retreive_branches(curve)
!    TYPE(curve_type) :: curve
!    REAL, DIMENSION(curve%length) :: x, y
!    INTEGER :: max_index_x, min_left_index, min_right_index, max_index
!    LOGICAL :: is_branched

!    x=curve%x
!    y=curve%y
!    delta=0.05*curve%length

!    IF (ABS(MINVAL(curve%x)) > ABS(MAXVAL(curve%x))) THEN
!      x=-x
!    END IF

!    max_index_x=fgsl_stats_max_index(x,1,curve%length)

!    is_branched=.TRUE.
!    min_left_index=1
!    min_right_index=curve%length
!    max_index=max_index_x

!    IF (max_index_x == 1 .OR. max_index_x == curve%length) THEN
!      is_branched=.FALSE.
!    ELSE
!      min_left_index_x=fgsl_stats_min_index(x(1:max_index_x),1,max_index_x)
!      min_right_index_x=fgsl_stats_min_index(x(max_index_x:curve%length),1, &
!        curve%length-max_index_x+1)+max_index_x-1
!      min_left_x=x(min_left_index_x)
!      min_right_x=x(min_right_index_x)
!      IF (min_left_x == MAX(min_left_x,min_right_x)) THEN
!        x_right=x(max_index_x:curve%length)
!        !x_right=x_right(x_right > min_left_x) !!!
!        min_index_x_right=fgsl_stats_min_index(x_right, 1, curve%length-max_index_x+1)
!        min_right_index_x=min_index_x_right+max_index_x-1
!      ELSE
!        x_left=x(1:max_index_x)
!        !x_left=x_left(x_left < min_right_x) !!!
!        max_index_x_left=fgsl_stats_max_index(x_left, 1, max_index_x)
!        min_left_index_x=max_index_x_left+1
!      END IF
!      IF (min_left_index_x >= min_right_index_x) is_branched=.FALSE.
!      index_y_left=MAX(min_left_index_x,max_index_x-2*delta)
!      index_y_right=MIN(min_right_index_x,max_index_x+2*delta)
!      max_index_y=fgsl_stats_max_index(y(index_y_left:index_y_right))
!      max_index_y=max_index_x-2*delta+max_index_y
!      IF (max_index_y == index_y_left .OR. max_index_y == index_y_right .OR. &
!        (max_index_y >= min_right_index_x .OR. &
!         max_index_y <= min_left_index_x)) THEN
!        is_branched=.FALSE.
!      ELSE
!        CALL extblas_find_local_minima(y(min_left_index_x:max_index_y),2, &
!          min_left_y,min_index_y_left)
!        min_left_index_y=0
!        IF (.TRUE.)
!          k=min_left_y_length
!          DO WHILE(k >= 1)
!            index0=min_index_y_left(k)+min_left_index_x-1
!            index1=MAX(1,index0-delta)
!            index2=MIN(curve%length, index0+delta)
!            f_status=fgsl_fit_linear(index1:index0,1,y(index1:index0),1, &
!              index0-index1+1,left_c0,left_c1,&
!              left_cov00,left_cov01,left_cov11,left_sumsq)
!            f_status_fgsl_fit_linear(index0:index2,1,y(index0:index2),1, &
!              index2-index0+1,right_c0,right_c1, &
!              right_cov00,right_cov01,right_cov11,right_sumsq)
!            IF (left_c1 < 0 .AND. right_c1 > 0) THEN
!              min_left_y=min_left_y(k)
!              min_left_index_y=index0
!              EXIT
!            END IF
!            k=k-1
!          END DO
!        END IF
!        IF (min_left_index_y == 0) THEN
!          min_left_y=y(min_left_index_x)
!          min_left_index_y=min_left_index_x
!        END IF
!        CALL extblas_find_local_minima(y(max_index_y:min_right_index_x),2, &
!          min_right_y,min_index_y_right)
!        min_right_index_y=0
!        IF (.TRUE.) THEN
!          k=1
!          DO WHILE(k <= min_right_y_length)
!            index0=min_index_y_right(k)+max_index_y-1
!            index1=MAX(1,index0-delta)
!            index2=MIN(curve%length,index0+delta)
!            f_status=fgsl_fit_linear(index1:index0,1,y(index1:index0),1, &
!              index0-index1+1,left_c0,left_c1,&
!              left_cov00,left_cov01,left_cov11,left_sumsq)
!            f_status_fgsl_fit_linear(index0:index2,1,y(index0:index2),1, &
!              index2-index0+1,right_c0,right_c1, &
!              right_cov00,right_cov01,right_cov11,right_sumsq)
!            IF (left_c1 < 0 .AND. right_c1 > 0) THEN
!              min_right_y=min_right_y(k)
!              min_right_index_y=index0
!              EXIT
!            END IF
!            k=k+1
!          END DO
!        END IF
!        IF (min_right_index_y == 0) THEN
!          min_right_y=y(min_right_index_x)
!          min_right_index_y=min_right_index_x
!        END IF
!        IF (y(min_left_index_y) >= y(min_right_index_y)) THEN
!          min_index_y_right=fgsl_stats_min_index(ABS(y(max_index_y:min_right_index_y)-min_left_y))
!      END IF
!    END IF

!  END SUBROUTINE gsr_retreive_branches

  SUBROUTINE gsr_save_curve(curve)
    TYPE(curve_type) :: curve
  END SUBROUTINE gsr_save_curve

END MODULE gsr

