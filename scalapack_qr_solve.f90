program test_scalapack

    use, intrinsic :: iso_fortran_env, only : dp => REAL64
    use mpi

    implicit none

    logical :: succeeded
    integer :: ierr, rank, nranks, ictxt
    
    integer :: mg, ng, ml, nl, mb, nb
    real(dp), allocatable :: A(:,:), B(:,:), X(:)

    call mpi_init(ierr)
    call blacs_pinfo(rank, nranks)
    call blacs_get(-1, 0, ictxt)
    call blacs_gridinit(ictxt, 'C', nranks, 1)
    call process_arguments()
    call random_seed()

    allocate(A(ml, nl))
    allocate(B(ml,1))
    allocate(X(nl), source=0.0_dp)
    call random_number(A)
    call random_number(B)
    write(*,*) "shape(A) = ", shape(A)
    write(*,*) "shape(B) = ", shape(B)
    write(*,*) "shape(X) = ", shape(X)

    write(*,*) "Solving ..."
    call SP_Matrix_QR_Solve()

    write(*,*) "X after :::"
    write(*,*) "", X

    call blacs_gridexit(ictxt)
    call mpi_finalize(ierr)

    contains

    subroutine SP_Matrix_QR_Solve()
        integer, dimension(9) :: descA, descB
        real(dp), dimension(:), allocatable :: tau, work

        if (rank+1 < nranks .and. mod(ml, mb) /= 0) &
            ! "SP_Matrix_QR_Solve: nrows is not a multiple of blocksize: "//ml//" "//mb
            call exit(3)
    
        call descinit (descA, mg, ng, mb, nb, 0, 0, ictxt, size(A, 1), ierr)
        call descinit (descB, mg,  1, mb,  1, 0, 0, ictxt, size(B, 1), ierr)

        call ScaLAPACK_pdgeqrf_wrapper(descA, A, tau, work)
        call ScaLAPACK_pdormqr_wrapper(descA, A, descB, B, tau, work)
        call ScaLAPACK_pdtrtrs_wrapper(descA, A, descB, B, .true.)
    
        call MatrixD_to_array1d(B, X)
    
    end subroutine SP_Matrix_QR_Solve

    subroutine ScaLAPACK_pdgeqrf_wrapper(desc, data, tau, work)
        integer, dimension(9), intent(in) :: desc
        real(dp), intent(inout), dimension(:,:) :: data
        real(dp), intent(out), dimension(:), allocatable :: tau
        real(dp), intent(out), dimension(:), allocatable :: work
      
        integer :: m, n, k, lwork, info
      
        m = mg
        n = ng
        k = min(m, n)
    
        call reallocate_real1d(tau, k)
        call reallocate_real1d(work, 1)
        call pdgeqrf(m, n, data, 1, 1, desc, tau, work, -1, info)
        lwork = work(1)
        call reallocate_real1d(work, lwork)
        call pdgeqrf(m, n, data, 1, 1, desc, tau, work, lwork, info)
      end subroutine ScaLAPACK_pdgeqrf_wrapper

    subroutine ScaLAPACK_pdormqr_wrapper(descA, A_data, descC, C_data, tau, work)
        integer, dimension(9), intent(in) :: descA, descC
        real(dp), intent(inout), dimension(:,:) :: A_data, C_data
        real(dp), intent(inout), dimension(:), allocatable :: tau
        real(dp), intent(inout), dimension(:), allocatable :: work
        
        integer :: m, n, k, lwork, info
        
        m = mg
        n = 1
        k = size(tau)
        
        call reallocate_real1d(work, 1)
        call pdormqr('L', 'T', m, n, k, A_data, 1, 1, descA, &
            tau, C_data, 1, 1, descC, work, -1, info)
        lwork = work(1)
        call reallocate_real1d(work, lwork)
        call pdormqr('L', 'T', m, n, k, A_data, 1, 1, descA, &
            tau, C_data, 1, 1, descC, work, lwork, info)
    end subroutine ScaLAPACK_pdormqr_wrapper
      
    subroutine ScaLAPACK_pdtrtrs_wrapper(descA, A_data, descB, B_data, cheat_nb_A)
        integer, dimension(9), intent(inout) :: descA, descB
        real(dp), intent(inout), dimension(:,:) :: A_data ! distributed triangular matrix
        real(dp), intent(inout), dimension(:,:)  :: B_data !
        logical, intent(in) :: cheat_nb_A
        
        integer, parameter :: mb_ = 5, nb_ = 6
        integer :: n, nrhs, info, nb
        
        n = min(mg, ng)
        nrhs = 1
        
        if (cheat_nb_A) then
            nb = descA(nb_)
            descA(nb_) = descA(mb_)
        end if
        
        ! A(lda,n+), B(ldb,nrhs+)
        call pdtrtrs('U', 'N', 'N', n, nrhs, A_data, 1, 1, descA, &
            B_data, 1, 1, descB, info)
        
        if (cheat_nb_A) then
            descA(nb_) = nb
        end if
        
    end subroutine ScaLAPACK_pdtrtrs_wrapper

    subroutine MatrixD_to_array1d(matrix, array)
        real(dp), intent(inout), dimension(:,:) :: matrix
        real(dp), intent(out), dimension(:) :: array
      
        integer :: nrows
      
        nrows = min(ng, size(array, 1))
        array(:nrows) = matrix(:nrows,1)
        array(nrows+1:) = 0.0_dp
    end subroutine MatrixD_to_array1d

    subroutine process_arguments()
        integer :: i, val
        character(len=32) :: arg
     
        mg = 1000
        ng = 100
        mb = 0
        nb = 100

        do i = 1, command_argument_count()
            call get_command_argument(i, arg)
            if (arg(3:3) /= "=") then
                ! "Argument does not conform to format __=..."
                call exit(1)
            end if
            read(arg(4:), '(i16)') val
            select case(arg(1:2))
                case('mg')
                    mg = val
                case('ng')
                    ng = val
                case('mb')
                    mb = val
                case('nb')
                    nb = val
                case default
                    ! "Argument not found."
                    call exit(2)
            end select
        end do

        write(*,*) "=== Parse arguments ==="
        write(*,*) "mg,ng = ", mg, " ", ng
        write(*,*) "mb,nb = ", mb, " ", nb

        ml = (mg + nranks - 1) / nranks
        nl = ng
        if (mb == 0) mb = ml
        if (nb == 0) nb = nl
        ml = increase_to_multiple(ml, mb)
        nl = increase_to_multiple(nl, nb)
        
        write(*,*) "=== Converted parameters ==="
        write(*,*) "mg,ng = ", mg, " ", ng
        write(*,*) "mb,nb = ", mb, " ", nb
        write(*,*) "ml,nl = ", ml, " ", nl
    end subroutine process_arguments

    function increase_to_multiple(a, m) result(res)
        integer, intent(in) :: a, m
        integer :: res
        res = (a / m) * m
        if (res < a) res = res + m
    end function increase_to_multiple

    subroutine reallocate_real1d(array, d1)
        real(dp), allocatable, dimension(:), intent(inout) :: array
        integer, intent(in)    :: d1

        if (allocated(array)) then
           if (size(array) == d1) return
           deallocate(array)
        end if
        allocate(array(d1))
      end subroutine reallocate_real1d

end program test_scalapack
