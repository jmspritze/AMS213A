!John Spritzer
!Project 1, create LU decompostion routine

module LU
    implicit none
contains

    !***********************************************************
    !***** LU Decompoition *********************
    subroutine LUD(am,an,bm,bn,A,At,L,U,b,bc,iSing)

        implicit none

        logical,intent(out) :: iSing
        integer :: am,an,bm,bn,n
        real, intent(in) :: A(am,an),b(bm,bn)
        real, intent(out) :: At(am,an),L(am,an),U(am,an)
        real, intent(out) :: bc(bm,bn)
        real s(am)
        real c, pivot, store
        integer i, j, k, ll, uu, v
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
            ll = k
            do j=k+1,am
                if((At(j,k)/s(j)) > pivot) then
                    pivot = (At(j,k)/s(j))
                    ll = j
                end if
            end do
            if(pivot .eq. 0.0) then ! Check if sigular
                iSing = .true.
                return
            end if

            if (ll .ne. k) then
                do j=k,am
                    store = At(k,j)
                    At(k,j) = At(ll,j)
                    At(ll,j) = store
                end do
                do j=k,bn
                    store = bc(k,j) !pivot rows k and l
                    bc(k,j) = bc(ll,j)
                    bc(ll,j) = store
                enddo
            end if

            ! uncomment to view intermediate matrixies
            !       print *, " "
            !       print *, "intermediate "
            !       do i=1,am
            !           write(*,10) (A(i,j), j=1,am)
            !     enddo

            ! Forward substitution
            do i=k+1,am
                c= At(i,k)/At(k,k)
                L(i,k) = c
                do uu = 1, bn
                    bc(i,uu)=bc(i,uu)-c*bc(k,uu)
                enddo
                At(i,k) = 0.0
                do j=k+1,am+1
                    At(i,j) = At(i,j)-c*At(k,j)
                end do
            end do
        enddo

        ! L matrix
        do i = 1,am
            L(i,i) =1.0
        enddo
        ! U matrix
        do j = 1,am
            do i =1,j
                U(i,j) = At(i,j)
            enddo
        enddo
    end subroutine LUD
end module LU
