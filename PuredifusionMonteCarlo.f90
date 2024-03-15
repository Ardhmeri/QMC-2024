Subroutine QuantumMonteCarloStep(a, N, nn, Z, Nuc, dt, nmax, energy,accep, tau, E_ref)

implicit none
    
! Input parameters
    real*8, intent(in) :: a, dt, tau, Nuc(3*nn)
    integer, intent(in) :: N, nn, Z
    integer*8, intent(in) :: nmax
    real*8, intent(in) :: E_ref
    
! Output parameters
    real*8, intent(out) :: energy, accep
    
! Local variables
    real*8 :: sq_dt, chi(3*N), B2_old, B2_new, prod, u, q
    real*8 :: psi_old, psi_new, argexpo
    real*8 :: R_old(3*N), R_new(3*N)
    real*8 :: B_old(3*N), B_new(3*N)
    real*8 :: local_energy, weight, normalization_factor, current_tau,e_loc , psi
    integer :: step_index, acceptance_count, i

! Compute square root of time step
    sq_dt = sqrt(dt)
    
! Initialize variables
    energy = 0.d0
    acceptance_count = 0
    normalization_factor = 0.d0
    weight = 1.d0
    current_tau = 0.d0
    call random_gauss(R_old, 3*N)

! Compute the drift term at the initial position
    call drift(a, R_old, N, nn, Nuc, B_old)

! Compute the square of the drift term
    B2_old = 0.d0
    do i = 1, 3*N
        B2_old = B2_old + B_old(i)*B_old(i)
    end do

! Compute the wave function at the initial position
    psi_old = psi(a, R_old, N, nn, Nuc)
    
! Perform Monte Carlo steps
    do step_index = 1, nmax
       
! Compute the local energy at the current position
        local_energy = e_loc(a, R_old, N, nn, Nuc, Z)
        
! Update the weight
        weight = weight * exp(-dt * (local_energy - E_ref))
        
! Accumulate the weight for normalization
        normalization_factor = normalization_factor + weight
        
! Accumulate the energy
        energy = energy + weight * local_energy
        
! Update the current tau value
        current_tau = current_tau + dt
        
! Reset weight and current tau value when tau is reached
        if (current_tau > tau) then
            weight = 1.d0
            current_tau = 0.d0
        end if
        
! Generate random displacement
        call random_gauss(chi, 3*N)
        
! Update the position using drift and diffusion
        R_new(:) = R_old(:) + dt * B_old(:) + chi(:) * sq_dt
        call drift(a, R_new, N, nn, Nuc, B_new)

! Compute the square of the new drift term
        B2_new = 0.d0
        do i = 1, 3*N
            B2_new = B2_new + B_new(i) * B_new(i)
        end do

! Compute the wave function at the new position
        psi_new = psi(a, R_new, N, nn, Nuc)

! Calculate the Metropolis acceptance ratio
        prod = 0.d0
        do i = 1, 3*N
            prod = prod + (B_new(i) + B_old(i)) * (R_new(i) - R_old(i))
        end do
        argexpo = 0.5d0 * (B2_new - B2_old) * dt + prod
        q = psi_new / psi_old
        q = exp(-argexpo) * q**2

! Determine whether to accept or reject the move based on the
! Metropolis criterion
        call random_number(u)
        if (u <= q) then

! Accept the move
            acceptance_count = acceptance_count + 1
            R_old(:) = R_new(:)
            B_old(:) = B_new(:)
            B2_old = B2_new
            psi_old = psi_new
        end if
    end do

! Compute the average energy and acceptance ratio
    energy = energy / normalization_factor
    accep = real(acceptance_count) / real(nmax)
    
end subroutine QuantumMonteCarloStep
