program build_jacobian
  double precision :: r_n(3, 3), r_ca(3, 3), r_c(3, 3)
  double precision :: jacobian

  ! Coordinates come from the first three residues of the "marked 
  ! loop" structure.

  r_n(:,1) =  (/  97.259, 39.712, 56.443 /)
  r_n(:,2) =  (/  99.345, 41.629, 56.722 /)
  r_n(:,3) =  (/ 100.357, 42.410, 59.363 /)

  r_ca(:,1) = (/  98.448, 39.501, 55.622 /)
  r_ca(:,2) = (/ 100.335, 42.688, 56.897 /)
  r_ca(:,3) = (/ 100.816, 42.188, 60.731 /)

  r_c(:,1) =  (/  99.429, 40.636, 55.791 /)
  r_c(:,2) =  (/ 101.057, 42.540, 58.216 /)
  r_c(:,3) =  (/  99.665, 42.264, 61.706 /)

  call find_jacobian(r_n, r_ca, r_c, jacobian)

  print *, 'Jacobian Determinant:', jacobian

end program build_jacobian

! Finds the determinant of the jacobian given nine position vectors.  
! This subroutine was adapted from the algortihms described by Dodd, 
! Boone, & Theodorou (1987) and Wu & Deem (1999).  The need for random 
! rotation has been eliminated by combining the 4th and 5th rows with a 
! correction factor to produce a 4x4 form of the jacobian.

subroutine find_jacobian(r_n, r_ca, r_c, jacobian)
  use vector
  implicit none

  double precision, intent(in) :: r_n(3,3), r_ca(3,3), r_c(3,3)
  double precision, intent(out) :: jacobian
  double precision :: axis(3,6), R1(3,6), R2(3,6), j(4,4), det
  integer :: i

  R1(:,1) =  r_n(:,1);    R2(:,1) = r_ca(:,1)
  R1(:,2) = r_ca(:,1);    R2(:,2) =  r_c(:,1)
  R1(:,3) =  r_n(:,2);    R2(:,3) = r_ca(:,2)
  R1(:,4) = r_ca(:,2);    R2(:,4) =  r_c(:,2)
  R1(:,5) =  r_n(:,3);    R2(:,5) = r_ca(:,3)
  R1(:,6) = r_ca(:,3);    R2(:,6) =  r_c(:,3)

  do i = 1,6
    axis(:,i) = normalize(R2(:,i) - R1(:,i))
  enddo

  do i = 1, 4
    ! jacobian elements (j11 - j14), (j21 - j24), (j31 - j34)
    j(1:3,i) = cross_product(axis(:,i), R1(:,6) - R1(:,i))

    ! jacobian elements (j41 - j45)
    j(4,  i) = dot_product(axis(:,i),                             &
                        cross_product(axis(:,5), axis(:,6)))
  end do

  ! compute determinant
  det = - j(4,1) * det_3(j(1:3,2), j(1:3,3), j(1:3,4))            &
        + j(4,2) * det_3(j(1:3,1), j(1:3,3), j(1:3,4))            &
        - j(4,3) * det_3(j(1:3,1), j(1:3,2), j(1:3,4))            &
        + j(4,4) * det_3(j(1:3,1), j(1:3,2), j(1:3,3))

  if (det == 0.d0) then
    print *, 'Null determinant, resetting to 10^-100'
    det = 10**(-100)
  end if

  jacobian = 1.0d0 / abs(det)
end subroutine find_jacobian

