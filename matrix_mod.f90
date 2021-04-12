module matrix_mod

  use number_types_mod , only: DP
  use point_mod        , only: Point
  use blas_sparse, only: duscr_begin, duscr_insert_entry , usds, uscr_end, &
  dusmv, blas_upper_triangular, ussp, ussv, duscr_insert_row
  !use blas, only: dtrsv
  !lapack, only: dgvem 

  implicit none

  contains 

  function GMRES(ABLAS, b, x0, tol, n, m_it) result(res)
    !real(DP) :: AA(:)
    !integer  :: JA(:), IA(n+1)
    real(DP) :: b(n), x0(n), tol 
    real(DP) :: res(n) 
    integer  :: n, m_it 
    real(DP) :: Q(n,m_it+1), H(m_it+1,m_it+1), sn(m_it), cs(m_it)
    real(DP) :: b_norm, r_norm, beta(m_it+1), error, y(m_it), ts, te
    integer :: i, j, k, k1, k2 

    integer :: ABlas, istat
    res = x0
    do j = 1, 1
      H = 0.0_DP
      cs = 0.0_DP
      sn = 0.0_DP 
      Q = 0.0_DP
      beta = 0.0_DP

      call dusmv(ABlas, res, Q(:,1), istat)

      Q(:,1) = b - Q(:,1) 


      b_norm = norm2(b)
      r_norm = norm2(Q(:,1))
      Q = Q / r_norm
      beta(1) = r_norm

      call cpu_time(ts)
      print *, "lap ready"
      do k = 1, m_it

        call Arnoldi(ABLAS, Q(:, 1:k+1), k, H(:,k))

        call apply_givens_rotation(H(1:k+1,k), cs, sn, k);

        beta(k+1) = -sn(k) * beta(k)
        beta(k)   =  cs(k) * beta(k)
        error = abs(beta(k+1)) / b_norm
        
        if (mod(k, 1000) == 0) then
          call cpu_time(te)
          print *, k, error, int(te)/3600, ":", mod(int(te),3600)/60, ":", mod(int(te),60)
        end if
        if (error < tol) exit
      end do 

      k = min(k, m_it)
      y(1:k) = HeissenbergSolver(H(1:k,1:k), beta(1:k))
      print *, k, error, int(te)/3600, ":", mod(int(te),3600)/60, ":", mod(int(te),60)
      open(501, file = "res.sol")
      do i = 1, n
        res(i) = res(i) + dot_product(Q(i, 1:k), y(1:k))
        write(501,*) res(i)
      end do    
      close(501)
    end do
    print *, error, k, n

    call usds(ABlas, istat)
  end function GMRES

  subroutine Arnoldi (ABLAS, Q, k, h)!(AA, JA, IA, Q, k, h)
    ! real(DP), intent(in   ) :: AA(:)
    ! integer , intent(in   ) :: JA(:), IA(:)
    integer , intent(in   ) :: ABLAS
    real(DP), intent(inout) :: h(:), Q(:,:)
    integer , intent(in   ) :: k
    integer :: i, k1, k2, istat

    call dusmv(ABLAS, Q(:,k), Q(:,k+1), istat)
    ! do i = 1, size(Q,1)
    !   k1 = IA(i)
    !   k2 = IA(i+1) - 1
    !   Q(i, k+1) = dot_product(AA(k1:k2), Q(JA(k1:k2), k))
    ! end do

    do i = 1, k   
      h(i) = dot_product(Q(:, k+1) , Q(:, i))
      Q(:, k+1) = Q(:, k+1) - h(i) * Q(:, i)
    end do
    h(k + 1) = norm2(Q(:, k+1))
    Q(:, k+1) = Q(:, k+1) / h(k + 1)
  end subroutine

  subroutine apply_givens_rotation(H, cs, sn ,k)
    real(DP), intent(inout) :: H(:), cs(:), sn(:)
    integer , intent(in   ) :: k
    real(DP) :: temp
    integer  :: i
    
    do i = 1, k-1
      temp   =  cs(i) * h(i) + sn(i) * h(i + 1)
      h(i+1) = -sn(i) * h(i) + cs(i) * h(i + 1)
      h(i)   = temp
    end do

    call givens_rotation(h(k), h(k + 1), cs(k), sn(k))

    h(k) = cs(k) * h(k) + sn(k) * h(k + 1)
    h(k + 1) = 0.0_DP
  end subroutine apply_givens_rotation

  subroutine givens_rotation(v1, v2, cs, sn)
    real(DP), intent(in    ) :: v1, v2
    real(DP), intent(   out) :: cs, sn
    real(DP) :: t 
    
    t = sqrt(v1**2 + v2**2)
    
    cs = v1 / t
    sn = v2 / t
    
  end subroutine givens_rotation

  function HeissenbergSolver(H, b) result(res)
    real(DP) :: H(:,:), b(:)
    real(DP) :: res(size(b))
    real(DP) :: c, A(size(b), size(b)), s(size(b))
    integer :: i, j, n, BBLAS, istat
    n = size(b)
    A = H
    s = b
    ! do i = 1, n-1
    !   c = A(i+1,i)/A(i,i)
    !   A(i+1, i+1:n) = A(i+1, i+1:n) - c * A(i, i+1:n)  
    !   s(i+1) = s(i+1) - c * s(i)
    ! end do
    call dtrsv("U", "N", "N", n, A, n, s, 1)
    !call duscr_begin(n, n, BBLAS,istat)
    !call ussp(BBLAS,blas_upper_triangular,istat)
    !do i = 1,n 
    !    call duscr_insert_row(BBLAS,i,A(i,i:n),[(j, j = i,n)],istat)
    !  if (mod(i,100)== 0) print *, i, n
    !end do
    !call uscr_end(BBLAS,istat )
    !call usSV(BBLAS,s, istat)
    ! call tbsv()
    ! call dgtrsm("L", "U", "N", "N", size(b), 1, 1.0, A, size(b), s, size(b))
    res = s 
    ! res(n) = s(n) / A(n,n)

    !call usds(BBLAS, istat)
    ! do i = n-1, 1, -1
    !  res(i) = (s(i) - dot_product(res(i+1:n), A(i,i+1:n))) / A(i,i)
    ! end do
  end function HeissenbergSolver

  subroutine swapColumns(A, b, j1, j2)
    real(DP), intent(inout) :: A(:,:), b(:)
    integer , intent(inout) :: j1, j2
    real(DP), dimension(size(A,1)) :: c1, c2 
    integer :: i
    c1 = A(:, j1)
    c2 = A(:, j2)

    A(:, j1) = c2
    A(:, j2) = c1

    c1 = A(j1, :)
    c2 = A(j2, :)
    
    A(j1, :) = c2
    A(j2, :) = c1

    do i = 1, size(b)
      b(i) = b(i) ! - BC value
    end do
  end subroutine swapColumns

end module