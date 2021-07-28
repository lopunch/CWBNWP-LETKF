module eigen
    implicit none

    integer                              :: lwork, liwork
    integer, allocatable, dimension(:)   :: iwork
    real*8,  allocatable, dimension(:)   :: work, eval
    real*8,  allocatable, dimension(:,:) :: evect

contains

    subroutine set_optimal_workspace_for_eigen(n)
        implicit none
        integer, intent(in)     :: n
        integer                 :: info

        allocate(work(1), iwork(1), eval(n), evect(n,n))

        call dsyevd('v', 'l', n, evect, n, eval, work, -1, iwork, -1, info)
        lwork  = int(work(1))
        liwork = iwork(1)

        deallocate(work, iwork)

        !optimal space for work 
        allocate(work(lwork), iwork(liwork))
    end subroutine set_optimal_workspace_for_eigen

    subroutine inverse_matrix(a, n)
        implicit none

        real*8, dimension(n,n), intent(in out) :: a
        integer,                intent(in)     :: n

        !local
        real*8, dimension(n,n)                 :: tmp
        integer                                :: i, info

        call dcopy(n*n, a, 1, evect, 1)
        call dsyevd('v', 'l', n, evect, n, eval, work, lwork, iwork, liwork, info)

        do concurrent(i=1:n)
            eval(i)  = 1d0 / eval(i)
            tmp(:,i) = evect(:,i) * eval(i)
        end do

        call dgemm('n', 't', n, n, n, 1d0, tmp, n, evect, n, 0d0, a, n)

    end subroutine inverse_matrix

    subroutine sqrt_matrix(a, n)
        implicit none

        real*8, dimension(n,n), intent(out) :: a
        integer,                intent(in)  :: n

        !local
        real*8, dimension(n,n)              :: tmp
        integer                             :: i

        do concurrent(i=1:n)
            tmp(:,i) = evect(:,i) * sqrt(eval(i))
        end do

        call dgemm('n', 't', n, n, n, 1d0, tmp, n, evect, n, 0d0, a, n)

    end subroutine sqrt_matrix

    subroutine destroy_eigen_array
        implicit none
        deallocate(eval, evect, work, iwork)
    end subroutine destroy_eigen_array

end module eigen
