module LinAl
    implicit none
contains

!***********************************************************************
!************************** Backsubstitution ***************************
!   Complete the backsubstion solution for Gaussian elimination. At is
!   the input matrix from the Gaussian Elimiation w/ pivots. x is the
!   solution vector and b is the product of Ax.
!***********************************************************************
    subroutine backsub(At,x,b,am,an,bm,bn)
        implicit none
        real,intent(in) :: At(am,an),b(bm,bn)
        real,intent(out) :: x(bm,bn)
        real :: bc(bm), xc(bm)
        real :: c
        integer :: i,j,k,am,an,bm,bn,n
10      format (10f5.2)
        xc=0.
        x=0.
        c=0.

        do n = 1, bn
            bc(:) = b(:,n)
            xc(bm) = bc(bm)/at(am,an)
            do i=bm-1,1,-1

                do j=i+1,am
                    c= c + at(i,j)*xc(j)
                end do
                xc(i) = (bc(i)- c)/at(i,i)
                c=0.0
            end do
            x(:,n) = xc(:)
            c=0.
        enddo
        print*, ' '
        print*, 'Solution x = '
        do k=1,bm
            write(*,*) (x(k,j), j=1,bn)
        enddo
    end subroutine backsub

!***********************************************************************
!************************** Matrix Trace *******************************
!   Computes the trace of a matrix. Sums the diagonal elements of an
!   mxm matrix. A is the matrix input.
!***********************************************************************
    subroutine trace(A,am,an)
        implicit none
        integer,intent(in) ::am,an
        real,intent(in) :: A(am,an)
        real :: traceVal
        integer :: i
11      format(f10.5)

        do i =1,am
            traceVal = traceVal + A(i,i)
        enddo

        print*, ' '
        print*, 'Trace of matrix = '
        write(*,11) traceVal
    end subroutine trace

!***********************************************************************
!************************** Euclidian Norm *****************************
!   Calulates the Euclidian Norm of a vector. Take the sqare root of
!   the sum of squares of an nx1 vector. V is the imput vector.
!***********************************************************************
    subroutine norm(V,am,j)
        implicit none
        integer,intent(in) ::am,j
        real,intent(in) :: V(am)
        real :: squareSum, normval
        integer :: i
11      format(f10.5)

        squareSum = 0.
        do i =1,am
            squareSum = squareSum + V(i)**2
        enddo
        normval = sqrt(squareSum)

        print*," "
        print*, 'Euclidian Norm of column ',j
        write(*,11) normval
    end subroutine norm

!***********************************************************************
!************************** Residual Matrix ****************************
!   Calulates the residual Error = (b-Ax) of a vector solution times the
!   original A matrix subracted from the original b.
!   A input matrix, x solution vector, b vector, emat is the error matrix
!***********************************************************************
    subroutine resid(A,x,b,emat,order)
        implicit none
        real,intent(in) :: A(order,order),b(order),x(order)
        real,intent(out) :: emat(order)
        real :: xt(order)
        integer,intent(in) :: order
        integer :: i,j,k
10      format(f10.5)
        emat=0.
        xt = matmul(A,x)
        do i = 1, order
            do j = 1, order
                emat(i) =  b(i) - xt(i)
            enddo
        enddo

        print*, ' '
        print*, 'Error Matrix = '
        do i = 1,order
            write(*,*) emat(i)
        enddo
    end subroutine resid

!***********************************************************************
!************************** LU Decompoition ****************************
!    LU factorization by Gaussian elimination (with partial pivoting)
!    am, an, bm, bn are sizes of matrices A a and b. A is the orignal
!    matrix and At is the temp matrix to work on. L is the lower and U
!    the upper triangle of A. Ising is the boolean value for a singular
!    matrix when found.
!***********************************************************************
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

!***********************************************************************
!********************* Guassing Elimination  ***************************
! Calcultes guassian elimination with pivots, A is the origial matrix,
! At is the temp matrix, b is the column vector and Ising is the
! logical result if the matrix is signular.
!***********************************************************************
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

!***********************************************************************
!*********************** Write Data to Screen  *************************
!   Writes matrix or vector to screen. A is a mxn matrix or mx1 vector.
!***********************************************************************
    subroutine writeToScreen(A,msize,nsize)
        implicit none
        integer :: msize,nsize,i,j,k
        real,intent(in) :: A(msize,nsize)
        do i=1,msize
            write(*,*) (A(i,j), j=1, nsize)
        enddo
    end subroutine writeToScreen

!***********************************************************************
!********************* Cholesky Factorization **************************
!  Compute Cholesky Factoriztion At*A*x = At*b. A is a vandermond matrix
!  of order n.
!***********************************************************************
    subroutine Chol(n,A,iSing)
        implicit none
        logical,intent(out) :: iSing
        integer :: n
        real, intent(inout) :: A(n,n)
        integer i, j, k

        do i = 1,n
            A(i,i) = sqrt(A(i,i))
            A(i+1:n,i) = A(i+1:n,i)/A(i,i)
            do j = i+1,n
                A(j:n,j) = A(j:n,j)- A(j,i)*A(j:n,i)
            enddo
        enddo
    end subroutine Chol

!***********************************************************************
!********************** Cholesky Backsubstitution **********************
!    Calulates the Cholesky backsubstition by first calulating the
!    Forward Substitution, solving Ly = b then Backward substitution,
!    Solving L*x = y. Input L is At*A, X is the solution vector, b is
!    the original b vector, n is the order/size of the columns and rows.
!***********************************************************************
    subroutine cholbacksub(L,x,b,y,n)
        implicit none
        integer :: i,j,k,n
        real :: sum1
        real,intent(in) :: L(n,n),b(n)
        real,intent(out) :: x(n),y(n)

        ! Forward Substitution, solving Ly = b
        do i = 1,n
            sum1 = b(i)
            do j = 1, i-1
                sum1 = sum1 - y(j)*L(i,j)
            end do
            y(i) = sum1/L(i,i)
        end do

        ! Backward substitution, Solving L*x = y
        do i = n, 1, -1
            if( L(i,i).lt. EPSILON(1.)) then
                stop
            end if
            do k = i + 1, n
                y(i) = y(i) - L(k,i)*x(k)
            enddo
            x(i) = y(i)/L(i,i)
        enddo
    end subroutine  cholbacksub

!***********************************************************************
!************************** Vandermond Matrix **************************
!    Produces a vandermond matrix the imput value of order. D is the
!    data from file (mx2). A is the produced vandermond matrix, and
!    order os the X^n order to produce.
!***********************************************************************
    subroutine vanderm(D, A, order, msize, nsize)
        real, dimension(:,:), allocatable :: A
        integer :: order, msize, nsize
        real, intent(in) :: D(msize,nsize)
        integer i, j, k
        ALLOCATE(A(msize, order))

        do i = 1,msize
            A(i,1) = 1
        enddo
        k =1
        do j = 2,order
            do i = 1, msize
                A(i,j) = D(i,1)**k
            enddo
            k = k+1
        enddo
    end subroutine Vanderm

!***********************************************************************
!********************* Build and write output file *********************
!    Builds output of the calulated polynomial and the original x value
!    from imported data. Additionally, the rms error is calulated by
!    taking the square root of the sum of differences between calc'd x
!    and the original data column. D is the original data, x is the sol
!   calc'd, T is temp, and W is output to write.
!***********************************************************************
    subroutine buildOut(D,x,msize, order)
        implicit none
        integer, intent(in) :: msize, order
        real, intent(in) :: D(msize,2),x(order)
        real :: T(msize),xt(msize), W(msize,2), rms, summ
        integer i, j, k
        print*," "
        T=0.
        rms = 0.
        summ = 0.

        print*,"writing data to output.txt..."

        !initalize all x1 values
        do i = 1, msize
            T(i) = x(1)
        enddo

        do i = 1, msize
            k = 1
            do j = 2, order
                T(i) = T(i) + x(j)*D(i,1)**k
                k = k+1
            enddo
        enddo

do i = 1, msize
   summ = summ + T(i) - D(i,2)
 enddo
 rms = sqrt(summ)
 print*, " "
 print*,"The RMS error between the Data and the interpolation: "
 print*,"RMS = ",rms
        ! create two column data file
        W(:,1) = D(:,1)
        W(:,2) = T(:)

        ! write file
        OPEN(UNIT=12, FILE="output.txt", ACTION="write", STATUS="replace")
        do i=1,msize
            write(12,*) (W(i,j), j=1,2)
        enddo
    end subroutine buildOut
end module LinAl
