program main
  implicit none
  real(8) De, a, mass, lambda, alpha, z, hg
  real(8) Nn, Cn
  real(8) xe, x_min, x_max
  real(8) h
  real(8) my_sum
  real(8), allocatable :: x(:), my_psi(:)
  integer(4) i, n, nSteps
  real(8), parameter :: hbar=1.054571800130d-34 ! J-S
  real(8), parameter :: Ese=1.602176620898d-19 ! electron charge C
  real(8), parameter :: Lse=1.0d-09 ! Nanometer to meter factor
  real(8), intrinsic :: gamma
  real(8), external :: fact
  open(1,file='in.dat', action='read')
  open(2,file='psi.dat', action='write')
  read(1,*) De  ! in the unit of eV
  read(1,*) a   ! in the unit of per nm
  read(1,*) xe  ! in the unit of nm
  read(1,*) x_min  ! in the unit of nm
  read(1,*) x_max  ! in the unit of nm
  read(1,*) nSteps ! segemnts between x_min and x_max
  read(1,*) mass   ! mass in kg
  read(1,*) n      ! nth eigenfunction, starts from 0 (zero)
  close(1)
  allocate(x(0:nSteps), my_psi(0:nSteps))

  h=abs(x_max-x_min)/dble(nSteps-1)

  lambda=dsqrt(2.0d0*mass*De*Ese)*a*Lse/hbar
  write(*,*) lambda
  alpha=2.0d0*(lambda-dble(n))-1.0d0

  Nn=dsqrt(fact(n)*alpha/gamma(2.0d0*lambda-dble(n)))
  write(*,11) n, Nn
  11 format(i5, e24.8)
  Cn=gamma(alpha+dble(n+1.0d0))/(gamma(alpha+1.0d0)*gamma(dble(n+1.0d0)))

  do i=0, nSteps
    x(i)=x_min+dble(i)*h
    z=2.0d0*lambda*exp(-(x(i)-xe)/a)
    call chgm(-dble(n), alpha+1.0d0, z, hg)
    my_psi(i)=Nn*z**(0.50d0*alpha)*exp(-0.5d0*z)*Cn*hg
!    write(2,10) x, my_psi(i)
  end do

  my_sum=0.0d0
  do i=0, nSteps
    my_sum = my_sum + my_psi(i)*my_psi(i)
  end do

  my_psi=my_psi/sqrt(my_sum)
  write(*,*) my_sum

  do i=0,nSteps
    write(2,12) x(i), my_psi(i)
  end do

  close(2)
  deallocate(x, my_psi)
  12 format(f8.4,f20.10)

end program main

