
MODULE types
  IMPLICIT NONE

  ! abstract types
  TYPE :: abstract_file_type
    CHARACTER(LEN=128) :: path
  END TYPE abstract_file_type

  TYPE :: abstract_data_type
    INTEGER :: length=0
    REAL :: sampling_rate=0
  END TYPE abstract_data_type

  TYPE :: abstract_curve_type
    REAL, DIMENSION(:), ALLOCATABLE :: x, y
    INTEGER :: length
  END TYPE abstract_curve_type

  ! config types
  TYPE :: config_row_type
    INTEGER :: flag, to_plot, to_compute_gsr, to_compute_mva, &
      to_compute_cylinder_model, to_compute_torus_model, &
      to_compute_hidalgo_model, to_compute_owens_model, to_save, begin_year, &
      begin_month, begin_day, begin_hour, begin_minute, end_year, end_month, &
      end_day, end_hour, end_minute, sampling_rate, min_theta, max_theta, &
      min_phi, max_phi, num_xpoints, order, begin_unixtime, end_unixtime
    CHARACTER(len=1) :: quality
    CHARACTER(len=4) :: spacecraft
    REAL :: step_ratio, min_y, max_y
  END TYPE config_row_type

  TYPE, EXTENDS(abstract_data_type) :: config_rows_type
    TYPE(config_row_type), DIMENSION(:), ALLOCATABLE :: rows
  END TYPE config_rows_type

  TYPE, EXTENDS(abstract_file_type) :: config_file_type
    TYPE(config_rows_type) :: config_struct
  END TYPE config_file_type

  ! data types
  TYPE :: data_row_type
    INTEGER :: year, month, day, hour, minute, second, data_unixtime
    REAL :: B, Bx, By, Bz, Vp, Vx, Vy, Vz, Pth, Np, Tp, Vth, beta
  END TYPE data_row_type

  TYPE, EXTENDS(abstract_data_type) :: data_rows_type
    TYPE(data_row_type), DIMENSION(:), ALLOCATABLE :: rows
  END TYPE data_rows_type

  TYPE, EXTENDS(abstract_file_type) :: data_file_type
    TYPE(data_rows_type) :: data_struct
  END TYPE data_file_type

  ! mva types
  TYPE :: mva_type
    REAL, DIMENSION(3, 3) :: matrix
    REAL, DIMENSION(3) :: x, y, z
    REAL :: lambda_x, lambda_y, lambda_z, criterion
  END TYPE mva_type

  TYPE, EXTENDS(mva_type) :: pmva_type
    REAL, DIMENSION(3, 3) :: p_matrix, pmp_matrix
  END TYPE pmva_type

  ! gsr and other types
  TYPE :: gsr_type
    REAL, DIMENSION(3) :: x, y, z
    INTEGER :: theta, phi
  END TYPE gsr_type

  TYPE :: dht_type
    REAL, DIMENSION(3) :: vht
    REAL :: cc
  END TYPE dht_type

  TYPE, EXTENDS(abstract_curve_type) :: branch_type
    INTEGER :: left_index,right_index
  END TYPE branch_type

  TYPE, EXTENDS(abstract_curve_type) :: curve_type
    LOGICAL :: is_branched
    TYPE(branch_type) :: branch_in
    TYPE(branch_type) :: branch_out
    REAL :: residue
    INTEGER :: branch_length
  END TYPE curve_type

  TYPE :: event_type
    TYPE(data_rows_type) :: data_struct
    TYPE(config_row_type) :: config_struct
    TYPE(dht_type) :: dht
    TYPE(gsr_type) :: gsr
    TYPE(mva_type) :: mvab
    TYPE(mva_type) :: mvub
    TYPE(pmva_type) :: pmvab
    REAL, DIMENSION(:,:), ALLOCATABLE :: residue
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: branch_length
  END TYPE event_type

  ! cm types
  TYPE cm_row_type
    REAL :: bx,by,bz,vx,vy,vz
  END TYPE cm_row_type

  TYPE, EXTENDS(abstract_data_type) :: cm_rows_type
    TYPE(cm_row_type), DIMENSION(:), ALLOCATABLE :: rows
  END TYPE cm_rows_type
END MODULE types

