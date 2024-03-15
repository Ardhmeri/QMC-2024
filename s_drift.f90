subroutine drift(a, Rtot, N, nn, Nuc, Bi)
  implicit none
  double precision, intent(in)   :: a
  integer, intent(in)            :: N, nn
  double precision, intent(in)   :: Rtot(3*N), Nuc(3*nn)
  double precision, intent(out)  :: Bi(3*N)  ! The drift vector
  double precision               :: r(3), b(3)
  double precision               :: distance_e_nuc, Den
  integer                        :: i, j

  do i = 1, 3*N, 3
     b = 0.d0  ! Reset drift vector for each electron
     Den = 0.d0  ! Reset denominator for normalization
     do j = 1, 3*nn, 3
        r(1:3) = Rtot(i:i+2) - Nuc(j:j+2)  ! Calculate distance vector between electron and nucleus
        distance_e_nuc = dsqrt(sum(r**2))  ! Compute distance between electron and nucleus
        if (distance_e_nuc > 0.d0) then
           b = b - a * r / distance_e_nuc * dexp(-a * distance_e_nuc)  !Update drift vector
           Den = Den + dexp(-a * distance_e_nuc)  ! Update denominator for normalization
        end if
     end do
     if (Den > 0.d0) then  ! Check if denominator is not zero
        b = b / Den  ! Normalize drift vector
     end if
     Bi(i:i+2) = b  ! Assign the calculated drift vector to output array
  end do
end subroutine drift
