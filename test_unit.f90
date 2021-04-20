module test_unit

  use number_types_mod , only: DP 
  use point_mod        , only: Point, idx, swapNodes
  use quadrilateral_mod, only: Square
  use param_mod        , only: Param
  use matrix_mod       , only: GMRES, swapColumns
  use ogpf
  use progress_mod
  use out_mod          , only: write_nodes
  use blas_sparse, only: duscr_begin, duscr_insert_entry , usds, uscr_end, dusmv, blas_upper_triangular, ussp, ussv 

  implicit none
  
  real(DP), parameter :: x1 = - 1.0_DP / sqrt(3.0), x2 = 1.0_DP / sqrt(3.0)
  real(DP), parameter :: w1 = 1.0_DP, w2 = 1.0_DP

  abstract interface 
  real(DP) function functBig(s, p, parametr, basis1, basis2) result(res)
  import :: Point, Square, Param, DP
        type(Square) :: s
        type(Point) :: p
        type(Param) :: parametr
        integer     :: basis1, basis2
      end function functBig
  end interface
  contains

  subroutine test_GMRES()
    real(DP) :: xxx(4), aaaaa(4,4), s(4)
    integer :: i, j
    integer :: AB, istat

    s =(/ 0.31746031746031772, 1.0910973084886129, 0.47204968944099379, -0.66321601104209804 /)
    xxx= 0.0_DP
    call duscr_begin(4, 4, AB, istat)
    do i = 1, 4
      do j = 1, 4
          aaaaa(i,j) = 2.0_DP*i-3.0_DP*j+1.0_DP*i*j
          if (mod(i+j,3)==0) aaaaa(i,j) = -1.0_Dp
          call duscr_insert_entry(AB, aaaaa(i,j), i, j, istat)
      end do
    end do

    call uscr_end(AB, istat)
    xxx = GMRES(AB, (/1.0_DP,3.0_DP,4.0_DP,-1.0_DP/), xxx, (10.0_dp)**(-10), 4, 4)

    print *, "TEST GMRES:"
    print *, s-xxx

    ! call usds(AB, istat)
    ! Solution
    !   0.31746031746031772     
    !   1.0910973084886129     
    !  0.47204968944099379     
    ! -0.66321601104209804  
  end subroutine test_GMRES

  ! Test the integral in a square that is not the reference
  subroutine test_Integral(parametr)
    type(Param), intent(in   ) :: parametr
    type(Point)  :: a1, a2, a3, a4
    type(Square) :: s 

    a1 = Point(2,3)
    a2 = Point(5,3)
    a3 = Point(5,5)
    a4 = Point(2,5)

    s = Square(a1, a2, a3, a4)
    
    print *, "TEST INTEGRAL:"
    print *, 84 - integralBig(s, testf, parametr, 1, 1)
  end subroutine test_Integral
  
  subroutine test_PointOrder(squares)
    type(Square), intent(in   ) :: squares(:)
    integer :: i, j, k

    print *, "TEST NODES:"
    do k = 1, size(squares)
      if (squares(k)%a(1)%x /= squares(k)%a(4)%x) then
        print *, "x coord in a1 a4", k, squares(k)%a(:)%x
        stop
      else if (squares(k)%a(2)%x /= squares(k)%a(3)%x) then
        print *, "x coord in a2 a3", k, squares(k)%a(:)%x
        stop
      else if (squares(k)%a(1)%y /= squares(k)%a(2)%y) then
          print *, "y coord in a1 a2", k, squares(k)%a(:)%y
          stop
      else if (squares(k)%a(3)%y /= squares(k)%a(4)%y) then
        print *, "x coord in a3 a4", k, squares(k)%a(:)%y
        stop
      end if
    end do
    PRINT *, "ALL OK"

  end subroutine test_PointOrder
  ! STUFF NEEDED FOR THE MAIN TESTS
  real(DP) function integralBig(self, f, parametr, k, l) result(res)
  type(square)       :: self
  procedure(functBig) :: f
  type(Param)         :: parametr
  integer             :: k, l
  res = f(self, Point(x1,x1), parametr, k ,l) * w1 * w1 * self%jacobDetInGaussPoint(1,1) + &
        f(self, Point(x1,x2), parametr, k ,l) * w1 * w2 * self%jacobDetInGaussPoint(1,2) + &
        f(self, Point(x2,x1), parametr, k ,l) * w2 * w1 * self%jacobDetInGaussPoint(2,1) + &
        f(self, Point(x2,x2), parametr, k ,l) * w2 * w2 * self%jacobDetInGaussPoint(2,2)  
end function integralBig

real(DP) function testf(s, p, parametr, k, l) result(res)
  type(Square) :: s  
  type(Point) :: p
  type(Point) :: ph
  type(Param) :: parametr 
  integer     :: k, l

  ph = s%refToMain(p)

  res =  ph%x * ph %y
 ! res = DX_basisF(l, p) * DX_basisF(k, p)
end function testf


real(DP) function basisF(A, p) result(res)
  integer :: A
  type(Point) :: p

  select case(A)
    case(1)
      res = linear_1D(1, p%x) * linear_1D(1, p%y)
    case(2)
      res = linear_1D(2, p%x) * linear_1D(1, p%y)
    case(3)
      res = linear_1D(2, p%x) * linear_1D(2, p%y)
    case(4)
      res = linear_1D(1, p%x) * linear_1D(2, p%y)
    case DEFAULT
      print *, "basisF basis can't have index", A, " in 2D"
      stop
  end select

end function basisF

real(DP) function DX_basisF(A, p) result(res)
  integer :: A
  type(Point) :: p

  select case(A)
    case(1)
      res = D_linear_1D(1, p%x) * linear_1D(1, p%y)
    case(2)
      res = D_linear_1D(2, p%x) * linear_1D(1, p%y)
    case(3)
      res = D_linear_1D(2, p%x) * linear_1D(2, p%y)
    case(4)
      res = D_linear_1D(1, p%x) * linear_1D(2, p%y)
    case DEFAULT
      print *, "basisF basis can't have index", A, " in 2D"
      stop
end select

end function DX_basisF

real(DP) function DY_basisF(A, p) result(res)
  integer :: A
  type(Point) :: p

  select case(A)
    case(1)
      res = linear_1D(1, p%x) * D_linear_1D(1, p%y)
    case(2)
      res = linear_1D(2, p%x) * D_linear_1D(1, p%y)
    case(3)
      res = linear_1D(2, p%x) * D_linear_1D(2, p%y)
    case(4)
      res = linear_1D(1, p%x) * D_linear_1D(2, p%y)
    case DEFAULT
      print *, "basisF basis can't have index", A, " in 2D"
      stop
end select

end function DY_basisF

real(DP) function linear_1D(A, p) result(res)
  integer  :: A
  real(DP) :: p

  if (A == 1) then
    res = (1.0_DP - p) / 2.0_DP
  else if (A == 2) then
    res = (1.0_DP + p) / 2.0_DP
  else 
    print *, "Linear Basis can't have index", A 
    stop
  end if
end function linear_1D

real(DP) function D_linear_1D(A, p) result(res)
  integer  :: A
  real(DP) :: p

  if (A == 1) then
    res = - 1.0_DP / 2.0_DP
  else if (A == 2) then
    res =   1.0_DP / 2.0_DP
  else 
    print *, "Linear Basis Derivative can't have index", A 
    stop
  end if
end function D_linear_1D
end module