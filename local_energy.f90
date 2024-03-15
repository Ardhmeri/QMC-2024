double precision function e_loc(a,Rtot,N,nn,Nuc,Z)
  implicit none !this function compute the local energy in R
  double precision, intent(in)  :: a
  integer, intent(in)           :: N,nn,Z
  double precision, intent(in)  :: Rtot(3*N),Nuc(3*nn)
  double precision, external    :: kinetic   ! we use the kinetic local energy
  double precision, external    :: potential ! and the potential calculated
                                             ! before
  e_loc = kinetic(a,Rtot,N,nn,Nuc)+potential(Rtot,N,nn,Z,Nuc)

end function e_loc
