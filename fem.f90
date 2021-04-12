program FEM
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

! DIRICHLET BOUNDARY CONDITIONS RANT, READ IF IN NEED TO UNDERSTAND
! Let us recall that at the end we are solving a system of equation Ax=b.
! The nodes that are in the D.B already have they value set so there is no 
! need to calculate its solution. We just need to be careful because they 
! influence the other nodes in their elements. For that we are going to 
! rearrenge the A matrix in the following way:
! [A_int A_int_bc]
! [0       I     ]
! the block A_int contains the matrix of the coefficients of only the interior
! nodes. The block A_int_bc contains the columns of the matrix A for 
! the correspondent node that is in the boundary no change needed, this is 
! because they only have influence in other nodes therefor the rows can be 
! changed to the identity row and b to the Dirichlet values

  type(Point) , allocatable :: nodes(:)
  type(Square), allocatable :: elements(:)
  type(Param)               :: parametr
  type(gpf) :: gp 
  real(DP), allocatable :: xgrid(:), ygrid(:)
  real(DP), allocatable :: x(:,:), y(:,:), z(:,:)
  real(DP), allocatable :: A(:,:), b(:), AA(:)
  real(DP), allocatable :: x0(:), sol(:)
  integer , allocatable :: JA(:), IA(:), IDBC(:)
  integer :: ABlas, istat

  real(DP) :: help, tStart, tEnd, xxx(4), aaaaa(4,4)
  integer :: nEl, nN, nDBC
  integer :: i, j, k, l 

  call cpu_time(tStart) 

  ! CHUNK OF CODE FOR TEST GMRES IF NEEDED
  ! xxx= 0.0_DP
  ! call duscr_begin(4, 4, ABlas, istat)
  ! do i = 1, 4
  !   do j = 1, 4
  !       aaaaa(i,j) = 2.0_DP*i-3.0_DP*j+1.0_DP*i*j
  !       if (mod(i+j,3)==0) aaaaa(i,j) = -1.0_Dp
  !       call duscr_insert_entry(ABlas, aaaaa(i,j), i, j, istat)
  !   end do
  ! end do
  ! do i = 1,4
  !   print *, aaaaa(i,:)
  ! end do
  ! call  progress_bar(i, nN-nDBC)
  ! call uscr_end(ABlas, istat)
  ! xxx = GMRES(ABlas, (/1.0_DP,3.0_DP,4.0_DP,-1.0_DP/), xxx, (10.0_dp)**(-10), 4, 4)
  ! Solution
  !   0.31746031746031772     
  !   1.0910973084886129     
  !  0.47204968944099379     
  ! -0.66321601104209804  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! PARAMETER SET UP AND GRID SET UP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call parametr%read_param("datainput")

  xgrid = linspace(parametr%xMin, parametr%xMax, parametr%nx)
  ygrid = linspace(parametr%yMin, parametr%yMax, parametr%ny)

  call meshgrid(x, y, xgrid, ygrid)

  nN = parametr%nx * parametr%ny 
  allocate(nodes(nN), z(size(x,1), size(x,2)), IDBC(parametr%ny))
  call cpu_time(tEnd)
  print *, "Mesh dode", int(tEnd)/3600, ":", mod(int(tEnd),3600)/60, ":", mod(int(tEnd),60)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NODES SET UP AND BOUNDARY FLAGGING !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  l = 0
  nDBC = 0
  do i = 1, size(x, 1)
    do j = 1, size(x, 2)
      nodes(idx(i,j)) = Point(x(i,j),y(i,j))
      if (nodes(idx(i,j))%x == parametr%xMin) then
        l = l + 1
        IDBC(l) = idx(i,j)
        nDBC = nDBC + 1
      end if 
    end do
  end do
  print *, IDBC
  call cpu_time(tEnd)
  call write_nodes(x, y)
  print *, "Nodes node", int(tEnd)/3600, ":", mod(int(tEnd),3600)/60, ":", mod(int(tEnd),60)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ELEMENTS SET UP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nEl = 0
  allocate(elements((size(x, 1)-1)*(size(x, 2)-1)))
  do i = 1, size(x,1) -1
    do j = 1, size(x,2) -1
      nEl = nEl + 1
      elements(nEl) = Square(nodes(idx(i,j)), nodes(idx(i,j+1)), &
      nodes(idx(i+1,j+1)), nodes(idx(i+1,j)))  
      call elements(nEl)%nodesIndex(idx(i,j), idx(i+1,j), &
      idx(i+1,j+1), idx(i,j+1))
      call elements(nEL)%set_K(parametr)
    end do
  end do
  call cpu_time(tEnd)
  print *, "Elements done", int(tEnd)/3600, ":", mod(int(tEnd),3600)/60, ":", mod(int(tEnd),60)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! MATRIX SET UP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate(A(nN, nN), b(nN), x0(nN), sol(nN))

  A = 0.0_DP
  b = 0.0_DP
  ! do i = 1, size(x,1)
  !   do j = 1, size(x,2)
  !     if (nodes(idx(i,j))%x < -20) then
  !       x0 (idx(i,j))= 0.0_DP
  !     else if (nodes(idx(i,j))%x > 20) then
  !       x0 (idx(i,j))= 1.0_DP
  !     else
  !       x0 (idx(i,j))= (2/3*parametr%nx - i ) * 0.0_DP + &
  !       (i - 1/3 * parametr%nx ) * 3/parametr%nx
  !     end if
  !   end do
  ! end do
  x0 = 0.0_DP
  sol = 0.0_DP
  
  
  
  do k = 1, nEl 
    
    do i = 1, 4
      do j = 1, 4
        A(elements(k)%aI(i), elements(k)%aI(j)) = A(elements(k)%aI(i), elements(k)%aI(j)) + &
        elements(k)%K(i,j)
      end do
      b(elements(k)%aI(i)) = b(elements(k)%aI(i)) + elements(k)%b(i)
    end do
  end do
  
  call cpu_time(tEnd)
  print *, "Matrix done", int(tEnd)/3600, ":", mod(int(tEnd),3600)/60, ":", mod(int(tEnd),60)

  l = nN
  do k = 1, nDBC 
    call swapColumns(A, b, iDBC(k), l)
    call swapNodes(nodes(IDBC(k)), nodes(l))
    
    A(l,:) = 0.0_DP
    A(l,l) = 1.0_DP
    b(IDBC(k)) = b(l)
    ! Set the right DBC
    b(l) = 0.0_DP
    
    l = l - 1
  end do 

   open(801, file = "prematrix")
   do i = 1, nN
     write(801, *) A(i,1), A(i,nN), b(i)
   end do
   close(801)
  ! open(901, file = "postmatrix")
  ! do i = 1, nN
  !   write(901, *) A(i,1), A(i,nN)
  ! end do
  ! close(901)
  call cpu_time(tEnd)
  print *, "BC applied", int(tEnd)/3600, ":", mod(int(tEnd),3600)/60, ":", mod(int(tEnd),60)
  
  open(601, file = "order.sol")
  do i = 1, size(x,1)
    do j = 1, size(x,2) 
      write(601,*) nodes(idx(i,j) )
  end do
end do
close(601)

  k = 1
  call duscr_begin(nN-nDBC, nN-nDBC, ABlas, istat)
  do i = 1, nN -nDBC
    do j = 1, nN -nDBC
      if (abs(A(i,j)) >= 10D-8) then
        k = k+1
        call duscr_insert_entry(ABlas, A(i,j), i, j, istat)
        end if
    end do
      if (mod(nN-nDBC,i)==1000) call progress_bar(i, nN-nDBC)
  end do
  call  progress_bar(i, nN-nDBC)
  call uscr_end(ABlas, istat)

  !   k = 1
  !   do i = 1, size(AA) 
  !     if(IA(k+1) == i) k = k + 1
  !   end do 
    
  
  ! allocate(AA(k), JA(k), IA(k))

  
  ! k = 1
  ! do i = 1, nN - nDBC
  !   IA(i) = k
  !   do j = 1, nN - nDBC
  !     if (abs(A(i,j)) >= 10D-8) then
  !       AA(k) = A(i,j)
  !       JA(k) = j
  !       IA(k) = i 
  !       k = k+1
  !     end if 
  !   end do
  ! end do
  !IA(nN -nDBC +1) = k
  call cpu_time(tEnd)
  print *, "Everything ready", int(tEnd)/3600, ":", mod(int(tEnd),3600)/60, ":", mod(int(tEnd),60)
  sol = 0.0_DP
  sol(1:nN-nDBC) = GMRES(ABlas, b(1:NN-nDBC), x0(1:nN-nDBC), (10.0_dp)**(-5), nN-nDBC, 2*nN)
  

  l = nN
  do k = 1, nDBC 
    call swapNodes(nodes(IDBC(k)), nodes(l))
    help = sol(IDBC(k))
    sol(IDBC(k)) = sol(l)
    sol(l) = help
    l = l - 1
  end do 

  open(701, file = "z.sol")
  do i = 1, size(x, 1)
    do j = 1, size(x, 2)
      z(i,j) = sol(idx(i,j))
      write(701,*) x(i,j), y(i,j), z(i,j)
    end do
  end do
  close(701)
  call cpu_time(tEnd)
  print *, "Everything done", int(tEnd)/3600, ":", mod(int(tEnd),3600)/60, ":", mod(int(tEnd),60)
  call gp%title("Temperature in a channel")
  call gp%xlabel("x Position")
  call gp%ylabel("y Position")

  call gp%options('unset key')
  !call gp%options('unset surface')

  !call gp%contour(x,y,z, palette='jet')
  call gp%surf(x, y, z, lspec='t "default color spec"' ) ! color palette: gnuplot default
  
  deallocate(nodes, elements, A, b, x0, sol, x, y, z)

end program FEM

