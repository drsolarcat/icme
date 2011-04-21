
MODULE mva
  USE fgsl
  USE types
  USE extblas

  IMPLICIT NONE
CONTAINS
  SUBROUTINE mva_compute_mvab(event)
    TYPE(event_type), INTENT(INOUT) :: event
    REAL, DIMENSION(3,3), TARGET :: evec
    REAL, DIMENSION(3), TARGET :: eval

    CALL mva_compute_mvab_matrix(event)
    CALL mva_compute_var_axes(event%mvab%matrix, evec, eval)

    event%mvab%x=evec(1,:)
    event%mvab%y=evec(2,:)
    event%mvab%z=evec(3,:)

    event%mvab%lambda_x=eval(1)
    event%mvab%lambda_y=eval(2)
    event%mvab%lambda_z=eval(3)

    event%mvab%criterion=event%mvab%lambda_y/event%mvab%lambda_z
  END SUBROUTINE mva_compute_mvab

  SUBROUTINE mva_compute_mvub(event)
    TYPE(event_type), INTENT(INOUT) :: event
    REAL, DIMENSION(3,3), TARGET :: evec
    REAL, DIMENSION(3), TARGET :: eval

    CALL mva_compute_mvub_matrix(event)
    CALL mva_compute_var_axes(event%mvub%matrix, evec, eval)

    event%mvub%x=evec(1,:)
    event%mvub%y=evec(2,:)
    event%mvub%z=evec(3,:)

    event%mvub%lambda_x=eval(1)
    event%mvub%lambda_y=eval(2)
    event%mvub%lambda_z=eval(3)

    event%mvub%criterion=event%mvub%lambda_y/event%mvub%lambda_z
  END SUBROUTINE mva_compute_mvub

  SUBROUTINE mva_compute_pmvab(event)
    TYPE(event_type), INTENT(INOUT) :: event
    REAL, DIMENSION(3,3) :: var_matrix, p_matrix, pmp_matrix, evec
    REAL, DIMENSION(3) :: eval, vht

    vht=event%dht%vht/extblas_norm(event%dht%vht)

    p_matrix(1,1)=1-vht(1)**2
    p_matrix(1,2)=-vht(1)*vht(2)
    p_matrix(1,3)=-vht(1)*vht(3)
    p_matrix(2,1)=-vht(2)*vht(1)
    p_matrix(2,2)=1-vht(2)**2
    p_matrix(2,3)=-vht(2)*vht(3)
    p_matrix(3,1)=-vht(3)*vht(1)
    p_matrix(3,2)=-vht(3)*vht(2)
    p_matrix(3,3)=1-vht(3)**2

    event%pmvab%p_matrix=p_matrix

    var_matrix=mva_compute_var_matrix( &
      event%data_struct%rows(:)%bx, &
      event%data_struct%rows(:)%by, &
      event%data_struct%rows(:)%bz, &
      event%data_struct%length)

    event%pmvab%matrix=var_matrix

    pmp_matrix=extblas_matrix_multiply(p_matrix, &
      extblas_matrix_multiply(var_matrix, p_matrix, 3, 3), 3, 3)

    CALL mva_compute_var_axes(pmp_matrix, evec, eval)

    event%pmvab%pmp_matrix=pmp_matrix

    event%pmvab%x=-event%dht%vht/extblas_norm(event%dht%vht)
    event%pmvab%y=evec(1,:)
    event%pmvab%z=extblas_vcross(event%pmvab%x, event%pmvab%y)
    event%pmvab%z=event%pmvab%z/extblas_norm(event%pmvab%z)
  END SUBROUTINE mva_compute_pmvab

  SUBROUTINE mva_compute_mvab_matrix(event)
    TYPE(event_type), INTENT(INOUT) :: event

    event%mvab%matrix=mva_compute_var_matrix( &
      event%data_struct%rows(:)%bx, &
      event%data_struct%rows(:)%by, &
      event%data_struct%rows(:)%bz, &
      event%data_struct%length)
  END SUBROUTINE mva_compute_mvab_matrix

  SUBROUTINE mva_compute_mvub_matrix(event)
    TYPE(event_type), INTENT(INOUT) :: event

    event%mvub%matrix=mva_compute_var_matrix( &
      event%data_struct%rows(:)%bx/event%data_struct%rows(:)%b, &
      event%data_struct%rows(:)%by/event%data_struct%rows(:)%b, &
      event%data_struct%rows(:)%bz/event%data_struct%rows(:)%b, &
      event%data_struct%length)
  END SUBROUTINE mva_compute_mvub_matrix

  FUNCTION mva_compute_var_matrix(data_x, data_y, data_z, data_length) &
    RESULT(var_matrix)
    INTEGER, INTENT(IN) :: data_length
    REAL, DIMENSION(data_length), INTENT(IN) :: data_x, data_y, data_z
    REAL, DIMENSION(3,3) :: var_matrix

    var_matrix(1,1)=fgsl_stats_mean(data_x*data_x, 1, data_length)- &
      fgsl_stats_mean(data_x, 1, data_length)* &
      fgsl_stats_mean(data_x, 1, data_length)
    var_matrix(1,2)=fgsl_stats_mean(data_x*data_y, 1, data_length)- &
      fgsl_stats_mean(data_x, 1, data_length)* &
      fgsl_stats_mean(data_y, 1, data_length)
    var_matrix(1,3)=fgsl_stats_mean(data_x*data_z, 1, data_length)- &
      fgsl_stats_mean(data_x, 1, data_length)* &
      fgsl_stats_mean(data_z, 1, data_length)
    var_matrix(2,1)=fgsl_stats_mean(data_y*data_x, 1, data_length)- &
      fgsl_stats_mean(data_y, 1, data_length)* &
      fgsl_stats_mean(data_x, 1, data_length)
    var_matrix(2,2)=fgsl_stats_mean(data_y*data_y, 1, data_length)- &
      fgsl_stats_mean(data_y, 1, data_length)* &
      fgsl_stats_mean(data_y, 1, data_length)
    var_matrix(2,3)=fgsl_stats_mean(data_y*data_z, 1, data_length)- &
      fgsl_stats_mean(data_y, 1, data_length)* &
      fgsl_stats_mean(data_z, 1, data_length)
    var_matrix(3,1)=fgsl_stats_mean(data_z*data_x, 1, data_length)- &
      fgsl_stats_mean(data_z, 1, data_length)* &
      fgsl_stats_mean(data_x, 1, data_length)
    var_matrix(3,2)=fgsl_stats_mean(data_z*data_y, 1, data_length)- &
      fgsl_stats_mean(data_z, 1, data_length)* &
      fgsl_stats_mean(data_y, 1, data_length)
    var_matrix(3,3)=fgsl_stats_mean(data_z*data_z, 1, data_length)- &
      fgsl_stats_mean(data_z, 1, data_length)* &
      fgsl_stats_mean(data_z, 1, data_length)
  END FUNCTION mva_compute_var_matrix

  SUBROUTINE mva_compute_var_axes(var_matrix, evec, eval)
    REAL, DIMENSION(3,3), TARGET, INTENT(IN) :: var_matrix
    REAL, DIMENSION(3,3), TARGET, INTENT(OUT) :: evec
    REAL, DIMENSION(3), TARGET, INTENT(OUT) :: eval
    INTEGER :: f_status
    TYPE(fgsl_matrix) :: f_var_matrix, f_evec
    TYPE(fgsl_vector) :: f_eval
    TYPE(fgsl_eigen_symmv_workspace) :: f_workspace

    f_var_matrix=fgsl_matrix_init(1.0)
    f_evec=fgsl_matrix_init(1.0)
    f_eval=fgsl_vector_init(1.0)

    f_status=fgsl_matrix_align(var_matrix, 3, 3, 3, f_var_matrix)
    f_status=fgsl_matrix_align(evec, 3, 3, 3, f_evec)
    f_status=fgsl_vector_align(eval, 3, f_eval, 3, 0, 1)

    f_workspace=fgsl_eigen_symmv_alloc(3)

    f_status=fgsl_eigen_symmv(f_var_matrix, f_eval, f_evec, f_workspace)
    f_status=fgsl_eigen_symmv_sort(f_eval, f_evec, FGSL_EIGEN_SORT_ABS_DESC)

    CALL fgsl_matrix_free(f_var_matrix)
    CALL fgsl_matrix_free(f_evec)
    CALL fgsl_vector_free(f_eval)
    CALL fgsl_eigen_symmv_free(f_workspace)
  END SUBROUTINE mva_compute_var_axes
END MODULE mva

