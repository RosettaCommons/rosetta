module vector
  implicit none
  contains

  ! Computes the vector dot product of the two given vectors.

  function dot_product (U, V)
    double precision, dimension(3), intent(in) :: U, V
    double precision :: dot_product

    dot_product = U(1)*V(1) + U(2)*V(2) + U(3)*V(3)

  end function dot_product

  ! Computes the vector cross product of the two given vectors.

  function cross_product (U, V)
    double precision, dimension(3), intent(in) :: U, V
    double precision, dimension(3) :: cross_product

    cross_product = (/                                            &
      U(2)*V(3) - U(3)*V(2),                                      &
      U(3)*V(1) - U(1)*V(3),                                      &
      U(1)*V(2) - U(2)*V(1) /)

  end function cross_product

  ! Returns a normalized copy of the given input vector.

  function normalize(V)
    double precision, dimension(3), intent(in) :: V
    double precision, dimension(3) :: normalize
    double precision :: magnitude

    magnitude = sqrt(V(1)*V(1) + V(2)*V(2) + V(3)*V(3))
    normalize = (/                                                &
      V(1) / magnitude,                                           &
      V(2) / magnitude,                                           &
      V(3) / magnitude /)

  end function normalize
      
  ! Computes the determinant of the given 3x3 matrix.

  function det_3 (M1, M2, M3)
    double precision, dimension(3), intent(in) :: M1, M2, M3
    double precision :: det_3

    det_3 =                                                       &
          + M1(1) * M2(2) * M3(3)                                 &
          - M1(1) * M3(2) * M2(3)                                 &
          - M2(1) * M1(2) * M3(3)                                 &
          + M2(1) * M3(2) * M1(3)                                 &
          + M3(1) * M1(2) * M2(3)                                 &
          - M3(1) * M2(2) * M1(3)

  end function det_3
end module
