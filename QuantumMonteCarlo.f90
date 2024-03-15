program QuantumMonteCarlo

implicit none

! Declare variables

 integer :: molecule_choice 
 real*8 :: a 
 integer :: num_elec, i , k 
 integer :: num_nuc
 integer :: atomic_number
 integer*8 :: num_iterations
 integer :: num_runs
 real*8 :: time_step
 real*8 :: propagation_length
 real*8 :: energy_reference
 integer :: run_index
 real*8, allocatable :: energies(:), acceptance_ratios(:)
 real*8 :: average_energy, energy_error
 real*8 :: average_acceptance, acceptance_error
 real*8, allocatable :: nucleus_positions(:)

! Display molecule options

 write(*,*) 'Select the molecule type:'
 write(*,*) '1: Hydrogen (H) atom'
 write(*,*) '2: Helium (He) atom'
 write(*,*) '3: Hydrogen ion (H2+)'
 write(*,*) '4: Hydrogen molecule (H2)'
 write(*,*) '5: Hydrogen ion (H3+)'
 read(*,*) molecule_choice 

! Determine filename based on molecule choice

select case (molecule_choice)
case (1)
   open(unit=1, file='H_atom.dat', status='old', action='read')
case (2)
    open(unit=1, file='He_atom.dat', status='old', action='read')
case (3)
    open(unit=1, file='H2+_ion.dat', status='old', action='read')
case (4)
    open(unit=1, file='H2_molecule.dat', status='old',action='read')
case (5)
    open(unit=1, file='H3+_ion.dat', status='old', action='read')
case default
    write(*,*) 'Invalid molecule choice.'
     stop
end select

! Read parameters from file

 read(1,*) a
 read(1,*) num_elec
 read(1,*) num_nuc
 read(1,*) atomic_number
 read(1,*) num_iterations
 read(1,*) num_runs
 read(1,*) time_step
 read(1,*) propagation_length
 read(1,*) energy_reference
 
 allocate(nucleus_positions(3*num_nuc))

 do i=1,3*num_nuc
        read(1,*) nucleus_positions(i)
 enddo
    
close(1)


! Allocate memory for storing results

allocate(energies(num_runs))
allocate(acceptance_ratios(num_runs))


! Calculate average energy and acceptance ratio along with error estimates
    
!call ave_error(energies,num_runs,average_energy, energy_error)
do k=1, num_runs
     call QuantumMonteCarloStep(a,num_elec,num_nuc,atomic_number,nucleus_positions,time_step,num_iterations,energies(k),acceptance_ratios(k),propagation_length,energy_reference)
enddo 

! Output results
call ave_error(energies,num_runs,average_energy, energy_error)   

write(*,*) 'Average energy: ', average_energy, ' +/- ', energy_error

call ave_error(acceptance_ratios,num_runs,average_energy, energy_error)

write(*,*) 'Average acceptance ratio: ',average_energy , ' +/- ',energy_error

end program QuantumMonteCarlo
