!John Spritzer
!Project 1, read data in order to solve Ax = b, where A is an n × n square matrix and b is an n-
!vector. The code reads and allocates memory for any nxn (within memory size).

module GaussE
    implicit none
contains

    !***********************************************************
    !***** Guassing Elimination *********************
    subroutine GaussElim(am,an,bm,bn,A,At,b,bc,iSing)

        implicit none

        logical,intent(out) :: iSing
        integer :: am,an,bm,bn,n
        real, intent(in) :: A(am,an),b(bm,bn)
        real, intent(out) :: At(am,an)
        real, intent(out) :: bc(bm,bn)
        real s(am)
        real c, pivot, store
        integer i, j, k, l, u, v
10      format(7f5.2)

        At(1:am,1:an) = A(1:am,1:an)
        bc = b
        ! Gaussian Elimiaton with pivoting
        do k=1, am
            do i=k,am
                s(i) = 0.0
                do j=k,am
                    s(i) = max(s(i),abs(At(i,j)))
                end do
            end do
            pivot = abs(At(k,k)/s(k)) ! find largest pivot
            l = k
            do j=k+1,am
                if((At(j,k)/s(j)) > pivot) then
                    pivot = (At(j,k)/s(j))
                    l = j
                end if
            end do
            if(pivot .eq. 0.0) then ! Check if sigular
                iSing = .true.
                return
            end if

            if (l .ne. k) then
                do j=k,am
                    store = At(k,j)
                    At(k,j) = At(l,j)
                    At(l,j) = store
                end do
                do j=k,bn
                    store = bc(k,j) !pivot rows k and l
                    bc(k,j) = bc(l,j)
                    bc(l,j) = store
                enddo
            end if

            ! uncomment to view intermediate matrixies
            !       print *, " "
            !       print *, "intermediate "
            !       do i=1,am
            !           write(*,10) (A(i,j), j=1,am)
            !     enddo

            ! simple gaussian elimation
            do i=k+1,am
                c= At(i,k)/At(k,k)
                do u = 1, bn
                    bc(i,u)=bc(i,u)-c*bc(k,u)
                enddo
                At(i,k) = 0.0
                do j=k+1,am+1
                    At(i,j) = At(i,j)-c*At(k,j)
                end do
            end do
        enddo
    end subroutine GaussElim
end module GaussE
