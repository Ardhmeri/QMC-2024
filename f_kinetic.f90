double precision function kinetic(a, Rtot, N, nn, Nuc)
  implicit none
  double precision, intent(in)  :: a
  integer, intent(in)           :: N, nn
  double precision, intent(in)  :: Rtot(3*N), Nuc(3*nn)
  double precision              :: distance_e_nuc
  double precision              :: r(3), Num, Den, Arg
  integer                       :: i, j

  kinetic = 0.d0  ! Initialize local kinetic energy
  Num = 0.d0     ! Initialize numerator for the kinetic energy calculation
  Den = 0.d0     ! Initialize denominator for the kinetic energy calculation
  Arg = 0.d0     ! Initialize argument for the kinetic energy calculation

  ! Compute local kinetic energy
  do i = 1, 3*N, 3
    do j = 1, 3*nn, 3
      r(1:3) = Rtot(i:i+2) - Nuc(j:j+2)   ! Calculate distance vector between electron and nucleus
      distance_e_nuc = dsqrt(sum(r**2))   ! Compute distance between electron and nucleus
      
      if (distance_e_nuc > 0.d0) then
        Num = Num + (a**2 - (2*a / distance_e_nuc)) * dexp(-a * distance_e_nuc)
! Update numerator
        Den = Den + dexp(-a * distance_e_nuc)
! Update denominator
      else
        kinetic = 1.d30  ! Set to a large positive value to indicate divergence
        return           ! Exit the function
      end if
    end do
    
    ! Calculate the argument for the local kinetic energy at each electron
    Arg = Num / Den
    
    ! Accumulate the local kinetic energy
    kinetic = kinetic + Arg
    
    ! Reset variables for the next electron
    Num = 0.d0
    Den = 0.d0
    Arg = 0.d0
  end do

  ! Multiply by -0.5 to get the final local kinetic energy
  kinetic = -0.5d0 * kinetic

end function kinetic
