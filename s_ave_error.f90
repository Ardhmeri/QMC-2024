Subroutine ave_error(x, n, average, error)
  implicit none
  integer, intent(in)            :: n
  real*8, intent(in)               :: x(n)
  real*8, intent(out)              :: average, error
  real*8                         :: variance

  ! Check if the array has at least one element
  if (n < 1) then
    stop 'Error: Array dimension is less than 1'
  ! If array has only one element, average is that element and error is 0
  else if (n == 1) then
    average = x(1)
    error = 0.0
  ! Calculate average and error for array with more than one element
  else
    average = sum(x(:)) / dble(n)
    variance = sum((x(:) - average)**2) / dble(n-1)
    error = dsqrt(variance / dble(n))
  end if

end subroutine ave_error
