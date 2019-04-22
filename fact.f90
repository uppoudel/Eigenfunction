real(8) function fact(n)
  implicit none
  integer(4) i, n
  real(8) prod
  prod = 1.0d0
  do i = 1, n
    prod = prod*dble(i)
  end do
  fact=prod
end function fact

