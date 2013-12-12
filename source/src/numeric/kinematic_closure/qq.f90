! Copyright (c) 2011-2013 by Evangelos A. Coutsias and Michael J. Wester
! Department of Mathematics and Statistics
! University of New Mexico, Albuquerque, New Mexico, USA
! Written by Evangelos A. Coutsias

! -----------------------------------------------------------------------------
!==============================================================================
MODULE dixon
!==============================================================================
USE globals
USE geom
!-----------------------------------------------------------------------------
CONTAINS
!subroutines
! Dixon_Resultant_Real_Roots ! triaxial_coefficients ! RECONSTRUCT ! triangle ! charpoly
!functions
! point_value2 ! point_value4 ! point_value6 ! point_value8 ! point_value16
! point_product ! point_product4x2 ! point_product4x4 ! point_product8x6
! poly_product ! poly_product2x2 ! poly_product4x2  ! poly_product4x4 ! poly_product4sq
!                poly_product6x6 ! poly_product12x4 ! poly_product8x8 ! poly_product8sq
!-----------------------------------------------------------------------------
SUBROUTINE Dixon_Resultant_Real_Roots( A, B, C, D, n_roots, u, method, flag)
   implicit none
!  real (DP), external :: p
!  real (DP), intent(in) :: alpha(3), delta(0:2), eta(3), theta(3), xi(0:2)
   integer, intent(in) :: method
   integer, intent(inout) :: flag
   integer,   intent(out) :: n_roots
   real (DP), intent(out) :: u(3, 16)
!  integer                   n_soln
!  real (DP)                u2(3, 16)
!-----------------------------------------------------------------------------
   integer i, j, n, info, k(16), ktmp
   real (DP) A(0:2, 0:2), B(0:2, 0:2), C(0:2, 0:2), D(0:2, 0:2)
   real (DP) R(0:2, 8, 8), gA(16, 16), gB(16, 16)
   real (DP) e_alpha_r(16), e_alpha_i(16), e_beta(16), e(16), etmp, v(16, 16)
   real (DP) denom
   real start_time, stop_time
!-----------------------------------------------------------------------------
!  do i = 0, 2
!     do j = 0, 2
!        A(i,j) = p(i,1,2, alpha,delta,eta,theta,xi                      &
!                        )*p(0,j,3, alpha,delta,eta,theta,xi             &
!                                 ) - p(i,0,2, alpha,delta,eta,theta,xi  &
!                                            )*p(1,j,3, alpha,delta,eta, &
!                                                       theta,xi)
!        B(i,j) = p(i,2,2, alpha,delta,eta,theta,xi                      &
!                        )*p(0,j,3, alpha,delta,eta,theta,xi             &
!                                 ) - p(i,0,2, alpha,delta,eta,theta,xi  &
!                                            )*p(2,j,3, alpha,delta,eta, &
!                                                       theta,xi)
!        C(i,j) = p(i,2,2, alpha,delta,eta,theta,xi                      &
!                        )*p(1,j,3, alpha,delta,eta,theta,xi             &
!                                 ) - p(i,1,2, alpha,delta,eta,theta,xi  &
!                                            )*p(2,j,3, alpha,delta,eta, &
!                                                       theta,xi)
!        D(i,j) = p(j,i,1, alpha,delta,eta,theta,xi)
!     end do
!  end do
!-----------------------------------------------------------------------------
if (method == 1) then
   do i = 0, 2
      R(i, :, :) = reshape(source =                                           &
         (/ 0.0_DP, A(0,i), A(1,i), A(2,i), 0.0_DP, B(0,i), B(1,i), B(2,i),   &
            A(0,i), A(1,i), A(2,i), 0.0_DP, B(0,i), B(1,i), B(2,i), 0.0_DP,   &
            0.0_DP, B(0,i), B(1,i), B(2,i), 0.0_DP, C(0,i), C(1,i), C(2,i),   &
            B(0,i), B(1,i), B(2,i), 0.0_DP, C(0,i), C(1,i), C(2,i), 0.0_DP,   &
            0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, D(0,i), D(1,i), D(2,i),   &
            0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, D(0,i), D(1,i), D(2,i), 0.0_DP,   &
            0.0_DP, D(0,i), D(1,i), D(2,i), 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP,   &
            D(0,i), D(1,i), D(2,i), 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP /),&
         shape = (/ 8, 8 /), order = (/ 2, 1 /))
   end do

   ! R_16 := expand(LinearAlgebra:-Determinant(R[0] + R[1]*u_3 + R[2]*u_3^2))
   ! e := sort([fsolve(R_16, u_3)])

   gA = 0.0_DP
   do i = 1, 8
      gA(i, 8+i) = 1.0_DP
   end do
   gA(9:16, 1:8 ) = -R(0, :, :)
   gA(9:16, 9:16) = -R(1, :, :)

   gB = 0.0_DP
   do i = 1, 8
      gB(i, i) = 1.0_DP
   end do
   gB(9:16, 9:16) =  R(2, :, :)

   ! e is the vector of eigenvalues and v is the matrix of eigenvectors.
   call cpu_time(start_time)
!  call la_ggev(gA, gB, e_alpha_r, e_alpha_i, e_beta, VR = v, INFO = info)
   call rgg(16, 16, gA, gB, e_alpha_r, e_alpha_i, e_beta, 1, v, info)
   call cpu_time(stop_time)
!  print '(/"      Eigenvecs elapsed time: ",f7.3," seconds")', &
!        stop_time - start_time

   if (info /= 0) then
      print '("Dixon_Resultant_Real_Roots: Eigenvector calculation failure,", &
&             " info = ",i3)', &
            info
   end if

   ! Select only the real eigenvalues, enforce that they contain no imaginary
   ! parts, and then sort them in order of increasing value.
   ! n is the number of real eigenvalues discovered.
   k = 0  ! Integer array; this line fills with zero.
   n = 0  ! Integer index variable.

   ! k holds indices of real eigenvalues (and eigenvectors).
   ! e holds the real eigenvalues themselves.
   
   do i = 1, 16
      if (abs(e_alpha_i(i)) < epsilon(0.0_DP)) then
         n = n + 1
         k(n) = i
         e(n) = e_alpha_r(i)/e_beta(i)
      end if
   end do

   ! Shitty n**2 sort of eigenvalues.
   do i = 1, n - 1
      do j = i + 1, n
         if (e(i) > e(j)) then
            etmp = e(i)
            e(i) = e(j)
            e(j) = etmp

            ktmp = k(i)
            k(i) = k(j)
            k(j) = ktmp
         end if
      end do
   end do

   ! Sort the eigenvectors into the same order as the eigenvalues and also
   ! enforce that they contain no imaginary parts.
   ! e is v(9, :) / v(1, :).
   n_roots = n
   do i = 1, n
      if (abs(v(1,k(i))) .gt. .00000001) then
         denom = 1/v(1,k(i))
         u(1, i) = v(2, k(i)) * denom
         u(2, i) = v(5, k(i)) * denom
      else ! still need to guard!!!
         denom = 1/v(3,k(i))
         u(1, i) = v(4, k(i)) * denom
         u(2, i) = v(7, k(i)) * denom
      endif

      ! u3 is just the eigenvalues; I already knew that.
      u(3, i) = e(i)
   end do
!      print*, 'using eig'
!      do i = 1, n_roots
!         print*, i, u(1,i), u(2,i), u(3,i)
!      enddo
else
       call charpoly(A, B, C, D, n_roots, u, flag)
!      print*, 'using charpoly'
!      do i = 1, n_roots
!         print*, i, u2(1,i), u2(2,i), u2(3,i)
!      enddo
   endif

end SUBROUTINE Dixon_Resultant_Real_Roots
!==============================================================================
