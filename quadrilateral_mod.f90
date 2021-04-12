module quadrilateral_mod

  use number_types_mod, only: DP 
  use point_mod, only: Point
  use param_mod, only: Param
  
  implicit none
  
  private
  
  public :: basisF
  
  ! The reference square is the square (-1,-1), (1,-1), (1,1), (-1,1)
  ! We are going to assume Gaussian quadrature of order 2
  real(DP), parameter :: x1 = - 1.0_DP / sqrt(3.0), x2 = 1.0_DP / sqrt(3.0)
  real(DP), parameter :: w1 = 1.0_DP, w2 = 1.0_DP

  
  type, public :: Square
  type(Point) :: a(4)
  integer     :: aI(4)
  real(DP) :: jacobDetInGaussPoint(2,2)
  
  ! Coefficients of the equations that come from this element
  ! i,j index means that it is integrated over the i,j basis functions
  real(DP) :: K(4,4)
  real(DP) :: b(4)
  contains
  
  procedure :: integral
  procedure :: integralBig
  procedure :: refToMain 
  procedure :: jacobDet
  procedure :: nodesIndex
  procedure :: set_K
  ! procedure :: area => area_quadrilateral
  ! procedure :: basis
  end type
  
  interface Square
  procedure :: newSquare
  end interface
  
  
  abstract interface 
  real(DP) function funct(p) result(res)
  import :: Point, DP
  type(Point) :: p
end function funct

end interface

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
    
    real(DP) function integral(self, f) result(res)
    class(square) :: self
    procedure(funct) :: f
    res = f(Point(x1,x1)) * w1 * w1 * self%jacobDetInGaussPoint(1,1) + &
    f(Point(x1,x2)) * w1 * w2 * self%jacobDetInGaussPoint(1,2) + &
    f(Point(x2,x1)) * w2 * w1 * self%jacobDetInGaussPoint(2,1) + &
    f(Point(x2,x2)) * w2 * w2 * self%jacobDetInGaussPoint(2,2) 
  end function integral
  
  real(DP) function integralBig(self, f, parametr, k, l) result(res)
  class(square)       :: self
  procedure(functBig) :: f
  type(Param)         :: parametr
  integer             :: k, l
  res = f(self, Point(x1,x1), parametr, k ,l) * w1 * w1 * self%jacobDetInGaussPoint(1,1) + &
        f(self, Point(x1,x2), parametr, k ,l) * w1 * w2 * self%jacobDetInGaussPoint(1,2) + &
        f(self, Point(x2,x1), parametr, k ,l) * w2 * w1 * self%jacobDetInGaussPoint(2,1) + &
        f(self, Point(x2,x2), parametr, k ,l) * w2 * w2 * self%jacobDetInGaussPoint(2,2)  
end function integralBig

  type(Square) function newSquare(a1, a2, a3, a4)
    type(Point) :: a1, a2, a3, a4

    newSquare%a(1) = a1 
    newSquare%a(2) = a2 
    newSquare%a(3) = a3 
    newSquare%a(4) = a4 

    newSquare%jacobDetInGaussPoint(1,1) = newSquare%jacobDet(Point(x1,x1))
    newSquare%jacobDetInGaussPoint(1,2) = newSquare%jacobDet(Point(x1,x2))
    newSquare%jacobDetInGaussPoint(2,1) = newSquare%jacobDet(Point(x2,x1))
    newSquare%jacobDetInGaussPoint(2,2) = newSquare%jacobDet(Point(x2,x2))

  end function newSquare

  type(Point) function refToMain(self, p) result(res)
    class(Square) :: self
    type(Point)   :: p
    real(DP) :: x, y 
    integer  :: i 

    x = sum([(basisF(i,p) * self%a(i)%x, i = 1,4)])
    y = sum([(basisF(i,p) * self%a(i)%y, i = 1,4)])

    res = Point(x,y)

  end function refToMain

  real(DP) function jacobDet(self, p) result(res)
    class(Square) :: self
    type(Point)   :: p
    real(DP)      :: dxdpsi1, dxdpsi2, dydpsi1, dydpsi2
    integer       :: i

    dxdpsi1 = sum([(DX_basisF(i,p) * self%a(i)%x, i = 1,4)])
    dxdpsi2 = sum([(DY_basisF(i,p) * self%a(i)%x, i = 1,4)])

    dydpsi1 = sum([(DX_basisF(i,p) * self%a(i)%y, i = 1,4)])
    dydpsi2 = sum([(DY_basisF(i,p) * self%a(i)%y, i = 1,4)])
    res = abs(dxdpsi1 * dydpsi2 - dxdpsi2 * dydpsi1)
  end function jacobDet

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

  subroutine set_K(self, parametr)
    class(Square), intent(inout) :: self
    type(Param)  , intent(in   ) :: parametr

    integer  :: i, j

    do i = 1, 4
      do j = 1,4
        self%K(i,j) = self%integralBig(val, parametr, i, j) 
      end do
        self%b(i)  = self%integralBig(bval, parametr, i, i) 
    end do

  end subroutine set_K

  real(DP) function val(s, p, parametr, k, l) result(res)
    type(Square) :: s  
    type(Point) :: p
    type(Param) :: parametr 
    integer     :: k, l

    type(Point) :: ph
    real(DP) :: x, y
    real(DP) :: xMin, xMax, m
    real(DP) :: yMin, yMax, a
    
    m    = parametr%m 
    xMin = parametr%xMin 
    xMax = parametr%xMax 
    yMin = parametr%yMin 
    yMax = parametr%yMax
    a    = parametr%a
    
    ph = s%refToMain(p)
    x = ph%x
    y = ph%y
    res = &!m * u(m, y, yMin, yMax) * basisF(l,p) * DX_basisF(k,p) - &
              DX_basisF(l,p) * DX_basisF(k,p)                    + &
          !m * v(m, y, yMin, yMax) * basisF(l,p) * DY_basisF(k,p) - &
          a * DY_basisF(l,p) * DY_basisF(k,p)
  end function val

  real(DP) function u(m, y, yMin, yMax) result(res)
    real(DP) :: y, yMin, yMax, m 

    res = 6.0_DP * m * (yMax - y) * (y - yMin)
  end function u

  real(DP) function v(m, x, yMin, yMax) result(res)
    real(DP) :: x, yMin, yMax, m 

    res = 0.0_DP
  end function v

  real(DP) function bval(s, p, parametr, k, l) result(res)
    type(Square) :: s  
    type(Point) :: p
    type(Point) :: ph
    type(Param) :: parametr 
    integer     :: k, l

    ph = s%refToMain(p)
    res = parametr%q * 0.1_DP * exp(-sqrt((ph%x+10.0_DP)**2 + &
     100.0_DP *(ph%y-0.5_DP)**2)) * basisF(k, p)
  end function bval
  
  subroutine nodesIndex(self, node1, node2, node3, node4) 
    class(Square), intent(inout) :: self
    integer      , intent(in   ) :: node1, node2, node3, node4

    self%aI(1) = node1
    self%aI(2) = node2
    self%aI(3) = node3
    self%aI(4) = node4

  end subroutine nodesIndex
end module quadrilateral_mod