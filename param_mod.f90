module param_mod

  use number_types_mod, only: DP

  implicit none 

  type Param
    integer  :: nx, ny
    real(DP) :: xMin, xMax
    real(DP) :: yMin, Ymax, a
    real(DP) :: gamma, beta, q, m 

    real(DP) :: LeF, Lez

    contains 

    procedure :: read_param
  end type Param

  contains 

  subroutine read_param(self, file)
    class(Param)    , intent(  out) :: self
    character(len=*), intent(in   ) :: file
    
    open(101, file = file)
      
      read(101, *) self%nx  , self%ny

      read(101, *) self%xMin, self%xMax
      read(101, *) self%yMin, self%Ymax, self%a
      read(101, *) self%gamma, self%beta, self%q, self%m 
  
      read(101, *) self%LeF, self%Lez

    close(101)

    self%a = (1.0_DP/self%a)**2
  end subroutine read_param

end module param_mod