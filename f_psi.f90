double precision function psi(a, Rtot, N, nn, Nuc)
  implicit none
  double precision, intent(in)  :: a
  integer, intent(in)           :: N, nn
  double precision, intent(in)  :: Rtot(3*N), Nuc(3*nn)
  integer                       :: i, j
  double precision              :: r(3), phi, Norm, distance_e_nuc
  
  psi = 1.d0
  phi = 0.d0
  Norm = 0.d0

  ! Compute wave function
  do i = 1, 3*N, 3
    phi = 0.d0  ! Reset phi for each electron
    do j = 1, 3*nn, 3
      r(1:3) = Rtot(i:i+2) - Nuc(j:j+2)  ! Calculate the distance vector betweenelectron and nucleus
      distance_e_nuc = dsqrt(sum(r**2))  ! Compute distance between electron and nucleus
      
      ! Calculate contribution to phi from electron-nucleus interaction
      phi = phi + dexp(-a * distance_e_nuc)
    end do
    
    psi = psi * phi  ! Multiply phi for each electron to get the total wavefunction
  end do

  ! Normalize the wave function
  Norm = 1.d0 / sqrt(real(N))  ! Calculate normalization factor
  psi = psi * Norm             ! Apply normalization factor to wave function
  
end function psi
