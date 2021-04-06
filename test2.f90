program test
  use number_types_mod, only: DP 
  use matrix_mod      , only: GMRES, HeissenbergSolver 

  implicit none

  real(DP) :: A(8,8), b(8), x(8), x0(8) = 0.0_DP
  real(DP) :: AA(9*8)
  integer :: JA(9*8), IA(9*8)
  integer :: k, i, j 

  
  A(1,:) = (/ 1.0,  2.0,  3.0,  0.0,  0.0,  0.0,  0.0,  0.0/)
  A(2,:) = (/ 2.0,  3.0, -5.0,  0.0,  0.0,  1.0,  0.0,  0.0/)
  A(3,:) = (/ 0.0, 13.0,  5.0,  5.0,  1.0,  0.0,  0.0, -2.0/)
  A(4,:) = (/ 0.0,  0.0,  4.0,  8.0,  2.0,  1.0,  0.0,  0.0/)
  A(5,:) = (/ 1.0,  0.0,  0.0,  5.0,  3.0,  7.0,  0.0,  0.0/)
  A(6,:) = (/ 4.0,  0.0,  5.0,  0.0,  3.0,  2.0,  1.0,  0.0/)
  A(7,:) = (/ 0.0,  0.0,  0.0,  4.0,  2.0,  2.0,  3.0,  0.0/)
  A(8,:) = (/ 0.0,  0.0,  0.0,  0.0,  0.0,  5.0,  4.0,  7.0/)

  b = (/12.0, 3.0, 0.0, -7.0, 2.0, 9.0, -5.0, 1.0/)
  k = 1
  do i = 1, 8
    IA(i) = k
    do j = 1, 8
      if (A(i,j) /= 0.0_DP) then
        AA(k) = A(i,j)
        JA(k) = j
        k = k+1
      end if 
    end do
  end do
  IA(9) = k 

  x(1:8) = GMRES(AA, JA, IA, b, x0, 10d-4, 8, 8)

  print *, x
end program test