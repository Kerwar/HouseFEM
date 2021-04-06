module out_mod
  use number_types_mod, only: DP
  use point_mod       , only: idx

  implicit none 

  contains

  subroutine write_nodes(x, y)
    real(DP), intent(in   ) :: x(:,:), y(:,:)
    integer :: i, j, n, m
    
    n = size(x,1)
    m = size(x,2)
    open(201, file = "x.sol")

    do i = 1, n 
      do j = 1, m
      write(201, "(F14.10)") x(i,j)
      end do
    end do 
    close(201) 

    open(301, file = "y.sol")

    do i = 1, n 
      do j = 1, m
      write(301, "(F14.10)") y(i,j )
      end do 
    end do 
    close(301) 

    open(401, file = "n.sol")

    do i = 1, n 
      do j = 1, m
      write(401, "(I14)") idx(i,j )
      end do 
    end do 
    close(401) 

  end subroutine write_nodes
end module out_mod