module point_mod

  use number_types_mod, only: DP 

  implicit none

  type :: Point
    real(DP) :: x, y 
  end type Point

  interface Point 
    procedure :: newPoint
  end interface

  contains

  type(Point) function newPoint(x,y)
    real(DP) :: x, y
    
    newPoint%x = x
    newPoint%y = y

  end function newPoint

  integer function idx(n,m) result(res)
    integer :: n
    integer :: m

    res = (n-1) * 121 + m
  end function idx

  subroutine swapNodes(n1, n2)
    type(Point), intent(inout) :: n1, n2
    type(Point) :: h 

    h = n1 
    n1 = n2 
    n2 = h
  end subroutine swapNodes
end module point_mod