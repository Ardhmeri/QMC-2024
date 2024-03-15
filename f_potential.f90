double precision function potential(Rtot, N, nn, Z, Nuc)
  implicit none
  double precision, intent(in)  :: Rtot(3*N)
  double precision, intent(in)  :: Nuc(3*nn)
  integer, intent(in)           :: Z, N, nn
  double precision              :: distance_e_nuc, distance_ee,distance_nn
  double precision              :: r(3)
  integer                       :: i, j, k

  potential = 0.d0

  ! Calculate potential between electrons and nuclei
  do i = 1, 3*N, 3
    do j = 1, 3*nn, 3
      r(1) = Rtot(i) - Nuc(j)
      r(2) = Rtot(i+1) - Nuc(j+1)
      r(3) = Rtot(i+2) - Nuc(j+2)
      distance_e_nuc = dsqrt(sum(r**2))

      if (distance_e_nuc > 0.d0) then
        potential = potential - (Z / distance_e_nuc)
      else
        potential = 1.d30  ! Set to a large positive value to indicate divergence
        return            ! Exit the function
      end if
    end do
  end do

  ! Calculate potential between electrons
  if (N == 2) then
    r(1) = Rtot(1) - Rtot(4)
    r(2) = Rtot(2) - Rtot(5)
    r(3) = Rtot(3) - Rtot(6)
    distance_ee = dsqrt(sum(r**2))
    
    if (distance_ee > 0.d0) then
      potential = potential + 1 / distance_ee
    else
      potential = 1.d30  ! Set to a large positive value to indicate divergence
      return            ! Exit the function
    end if
  end if

  ! Calculate potential between nuclei
  if (nn > 1) then
    do i = 1, 3*nn, 3
      do j = i + 3, 3*nn, 3
        r(1) = Nuc(i) - Nuc(j)
        r(2) = Nuc(i+1) - Nuc(j+1)
        r(3) = Nuc(i+2) - Nuc(j+2)
        distance_nn = dsqrt(sum(r**2))
        
        if (distance_nn > 0.d0) then
          potential = potential + Z*Z / distance_nn
        else
          potential = 1.d30  ! Set to a large positive value to indicatedivergence
          return            ! Exit the function
        end if
      end do
    end do
  end if

end function potential
