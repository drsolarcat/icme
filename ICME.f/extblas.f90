
MODULE extblas
  IMPLICIT NONE
CONTAINS

  ! Cross product
  ! vector cross product
  FUNCTION extblas_vcross(a,b) RESULT(c)
    REAL, DIMENSION(3), INTENT(IN) :: a,b
    REAL, DIMENSION(3) :: c

    c=RESHAPE(extblas_mcross(RESHAPE(a,(/1,3/)),RESHAPE(b,(/1,3/))), (/3/))
  END FUNCTION extblas_vcross

  ! multiple vectors (matrix) cross product
  FUNCTION extblas_mcross(a, b) RESULT(c)
    REAL, DIMENSION(:,:), INTENT(IN) :: a, b
    REAL, DIMENSION(:,:), ALLOCATABLE :: c

    IF (ANY(SHAPE(a) /= SHAPE(b))) THEN
      WRITE(*,*) 'Error: passed arrays have different shapes'
      STOP
    END IF

    IF (SIZE(a,2) /= 3) THEN
      WRITE(*,*) 'Error: passed vectors are not 3d'
      STOP
    END IF

    ALLOCATE(c(SIZE(a,1), 3))

    c(:,1)=a(:,2)*b(:,3)-a(:,3)*b(:,2)
    c(:,2)=a(:,3)*b(:,1)-a(:,1)*b(:,3)
    c(:,3)=a(:,1)*b(:,2)-a(:,2)*b(:,1)
  END FUNCTION extblas_mcross

  ! expand multiple vectors (matrix) into component vectors
  SUBROUTINE extblas_expand(a, ax, ay, az)
    REAL, DIMENSION(:,:), INTENT(IN) :: a
    REAL, DIMENSION(SIZE(a,1)), INTENT(OUT) :: ax, ay, az

    ax=a(:,1)
    ay=a(:,2)
    az=a(:,3)
  END SUBROUTINE extblas_expand

  ! multiple vectors (matrix) cross product with following expansion into
  ! components vectors
  SUBROUTINE extblas_mcross_expand(a,b,cx,cy,cz)
    REAL, DIMENSION(:,:), INTENT(IN) :: a,b
    REAL, DIMENSION(SIZE(a,1)), INTENT(OUT) :: cx,cy,cz

    CALL extblas_expand(extblas_mcross(a,b),cx,cy,cz)
  END SUBROUTINE extblas_mcross_expand

  ! multiple vectors cross product
  SUBROUTINE extblas_mvcross_expand(ax,ay,az,bx,by,bz,cx,cy,cz)
    REAL, DIMENSION(:), INTENT(IN) :: ax,ay,az,bx,by,bz
    REAL, DIMENSION(SIZE(ax,1)), INTENT(OUT) :: cx,cy,cz
    INTEGER :: n

    n=SIZE(ax,1)
    CALL extblas_mcross_expand(RESHAPE((/ax,ay,az/),(/n,3/)), &
      RESHAPE((/bx,by,bz/),(/n,3/)), cx, cy, cz)
  END SUBROUTINE extblas_mvcross_expand

  ! other
  FUNCTION extblas_matrix_multiply(a, b, m, n) RESULT(c)
    INTEGER, INTENT(IN) :: m, n
    REAL, DIMENSION(m,n), INTENT(IN) :: a
    REAL, DIMENSION(n,m), INTENT(IN) :: b
    REAL, DIMENSION(m,m) :: c
    INTEGER :: i, k

    DO i=1,m
      DO k=1,m
        c(i,k)=SUM(a(i,:)*b(:,k))
      END DO
    END DO
  END FUNCTION extblas_matrix_multiply

  FUNCTION extblas_norm(v) RESULT(norm)
    REAL, DIMENSION(3), INTENT(IN) :: v
    REAL :: norm

    norm=SQRT(v(1)**2+v(2)**2+v(3)**2)
  END FUNCTION extblas_norm

  FUNCTION extblas_normalize(v) RESULT(nv)
    REAL, DIMENSION(3), INTENT(IN) :: v
    REAL, DIMENSION(3) :: nv

    nv=v/extblas_norm(v)
  END FUNCTION extblas_normalize

  FUNCTION extblas_rotate_vector(a, b, phi) RESULT(c)
    REAL, DIMENSION(3), INTENT(IN) :: a, b
    REAL, INTENT(IN) :: phi
    REAL, DIMENSION(3) :: c, b_abs

    b_abs=b/extblas_norm(b)
    c=a*COS(phi)+DOT_PRODUCT(a, b_abs)*(1-COS(phi))*b_abs+ &
      extblas_vcross(b_abs, a)*SIN(phi)
  END FUNCTION extblas_rotate_vector

  FUNCTION extblas_vprojection(a, x, y, z) RESULT(b)
    REAL, DIMENSION(3), INTENT(IN) :: a, x, y, z
    REAL, DIMENSION(3) :: b

    b(1)=DOT_PRODUCT(a, x)
    b(2)=DOT_PRODUCT(a, y)
    b(3)=DOT_PRODUCT(a, z)
  END FUNCTION extblas_vprojection

  SUBROUTINE extblas_mvprojection(ax, ay, az, x, y, z, bx, by, bz)
    REAL, DIMENSION(:), INTENT(IN) :: ax, ay, az
    REAL, DIMENSION(3), INTENT(IN) :: x, y, z
    REAL, DIMENSION(SIZE(ax,1)), INTENT(OUT) :: bx, by, bz
    INTEGER :: i

    DO i=1,SIZE(ax,1)
      bx(i)=DOT_PRODUCT((/ax(i),ay(i),az(i)/),x)
      by(i)=DOT_PRODUCT((/ax(i),ay(i),az(i)/),y)
      bz(i)=DOT_PRODUCT((/ax(i),ay(i),az(i)/),z)
    END DO
  END SUBROUTINE extblas_mvprojection

  FUNCTION simp(n,h,f) RESULT(s)
    INTEGER, INTENT (IN) :: n
    REAL, INTENT (IN) :: h
    REAL, INTENT (IN), DIMENSION (n) :: f
    REAL :: s0,s1,s2,s
    INTEGER :: i

    s=0.0
    s0=0.0
    s1=0.0
    s2=0.0

    DO i=2,n-1,2
      s1=s1+f(i-1)
      s0=s0+f(i)
      s2=s2+f(i+1)
    END DO

    s=h*(s1+4.0*s0+s2)/3.0

    IF (MOD(n,2) == 0) s=s+h*(5.0*f(n)+8.0*f(n-1)-f(n-2))/12.0
  END FUNCTION simp

  FUNCTION extblas_ones_vector(n) RESULT(ones)
    INTEGER :: n
    REAL, DIMENSION(n) :: ones

    ones = 1.0
  END FUNCTION extblas_ones_vector
END MODULE extblas

