subroutine random_gauss(z, k)
  implicit none
  integer, intent(in)            :: k
  double precision, intent(out)  :: z(k)
  double precision               :: u(k+1)
  double precision, parameter    :: two_pi = 2.0 * acos(-1.0)
  integer                        :: i

  call random_number(u)
  
  do i = 1, k, 2
    z(i) = dsqrt(-2.0 * dlog(u(i)))
    z(i+1) = z(i) * dsin(two_pi * u(i+1))
    z(i) = z(i) * dcos(two_pi * u(i+1))
  end do

  if (mod(k, 2) == 1) then
    z(k) = dsqrt(-2.0 * dlog(u(k)))
    z(k) = z(k) * dcos(two_pi * u(k+1))
  end if

end subroutine random_gauss
