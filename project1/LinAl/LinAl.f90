program LinAl

    use writeScreen
    use GaussE
    use LU

    implicit none
    integer :: am,an,bm,bn,io,i,j,k
    integer :: n
    real :: z,w,c,p,r,q
    logical :: iSing = .false.
    !*********** Functions ***********************************
    real, external :: trace
    real, external :: norm
    real, external :: backsub
    real, external :: resid
    real, external :: FSolve
    !************* Arrays ***********************************
    real, dimension(:,:), allocatable :: A,At,b,bc,x,emat,y
    real, dimension(:), allocatable :: V,xr
    real, dimension(:,:), allocatable :: L,U, lw
    !**********************************************************
    !********Get Number of Columns and Rows for MxN matrix *******
    open(10,file='Amat.dat',form='formatted',status = 'old')
    ! Read the matrix dimensions
    read(10,*) am,an
    open(12,file='Bmat.dat',form='formatted',status = 'old')
    ! Read the matrix dimensions
    read(12,*) bm,bn
    !**** ALLOCATE ALL MATRICIES ********************************
    ALLOCATE(A(am,an))
    ALLOCATE(At(am,an))
    ALLOCATE(x(bm,bn))
    ALLOCATE(b(bm,bn))
    ALLOCATE(bc(bm,bn))
    ALLOCATE(V(am))
    ALLOCATE(emat(bm,bn))
    A=0.
    b=0.
    ALLOCATE(L(am,an))
    ALLOCATE(Lw(am,an))
    ALLOCATE(U(am,an))
    ALLOCATE(y(bm,bn))
    !************************************************************
10  format(7f5.2)
    L=0.
    U=0.
    !y=0.
    x=0.
    At=0.
    emat=0.

    do i=1,am
        read(10,*) (A(i,j),j=1,am)
    enddo
    close(10)

    open(12,file = 'Bmat.dat')
    do k=1,bm
        read(12,*) (b(k,j),j=1,bn)
    enddo
    close(2)

    print*," "
    print*, "************************ Matrix Trace **************************"

    z = trace(A,am,an) !matrix trace

    print*," "
    print*, "************************ Matrix Norm ***************************"
    !matrix norm
    do  j = 1,am
        do i = 1,an
            V(i) = A(i,j)
        enddo
        w = norm(V,am,j)
    enddo
    print*," "
    print*," "
    print*, "************************ GAUSSING ELIM *************************"

    print*," "
    print*, "************************ Matrix A,b ****************************"
    call writeToScreen(am,an,bm,bn,A,b) ! Write Data to screen
    call GaussElim(am,an,bm,bn,A,At,b,bc,iSing)
    if(iSing .eqv. .true.) then
        write(*,*) "Matrix is singular"
        stop
    end if
    print*," "
    print*, "************************ Reduced Matrix A,b *********************"
    call writeToScreen(am,an,bm,bn,At,bc) ! Write Data to screen
    print*," "
    print*, "************************ Backsub ********************************"
    p = backsub(At,x,bc,am,an,bm,bn)
    print*," "
    print*, "************************ Residual Ax-B **************************"

    r = resid(A,x,b,emat,am,an,bm,bn)

    print*," "
    print*, "************************ Error matrix norm **********************"
        !matrix norm
    V=0
    do  j = 1,bm
        do i = 1,bn
            V(i) = emat(i,j)
        enddo
        w = norm(V,bm,j)
    enddo

    print*," "
    print*," "
    print*, "************************ LU DECOMP ******************************"
    print*," "
    print*, "************************ Reduced Matrix A,b *********************"
    call LUD(am,an,bm,bn,A,At,L,U,b,bc,iSing)

    print*, ' '
    print*, 'U Matrix = '
    do i=1,am
        write(*,10) (U(i,j), j=1,am)
    enddo

    print*, ' '
    print*, 'L Matrix = '
    do i=1,am
        write(*,10) (L(i,j), j=1,am)
    enddo
    q = backsub(U,x,bc,am,an,bm,bn)
    print*," "
    print*, "************************ Residual Ax-B **************************"

    r = resid(A,x,b,emat,am,an,bm,bn)

    print*," "
    print*, "************************ Error matrix norm **********************"
        !matrix norm
    V=0
    do  j = 1,bm
        do i = 1,bn
            V(i) = emat(i,j)
        enddo
        w = norm(V,bm,j)
    enddo


end program LinAl

!********** Backsubstitution function ************************
real function backsub(At,x,b,am,an,bm,bn)
    implicit none
    real,intent(in) :: At(am,an),b(bm,bn)
    real,intent(out) :: x(bm,bn)
    real :: bc(bm), xc(bm)
    real :: c
    integer :: i,j,k,am,an,bm,bn,n
10  format (10f5.2)
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
        write(*,10) (x(k,j), j=1,bn)
    enddo
end function backsub

!!***** Compute Trace of matrix A*******************************
real function trace(A,am,an)
    implicit none
    integer,intent(in) ::am,an
    real,intent(in) :: A(am,an)
    real :: traceVal
    integer :: i
11  format(f10.5)

    do i =1,am
        traceVal = traceVal + A(i,i)
    enddo

    print*, ' '
    print*, 'Trace of matrix = '
    write(*,11) traceVal

end function trace

!!!***** Compute Euclidian Norm of Column vector*******************************
real function norm(V,am,j)
    implicit none
    integer,intent(in) ::am,j
    real,intent(in) :: V(am)
    real :: normVal, squareSum
    integer :: i
11  format(f10.5)

    do i =1,am
        squareSum = squareSum + V(i)**2
    enddo
    normVal = sqrt(squareSum)
    squareSum = 0.
    print*," "
    print*, 'Euclidian Norm of column ',j
    write(*,11) normVal
end function norm

!!***** Compute residual of matrix Ax - b*******************************
real function resid(A,x,b,emat,am,an,bm,bn)
    implicit none
    real,intent(in) :: A(am,an),b(bm,bn),x(bm,bn)
    real,intent(out) :: emat(bm,bn)
    real :: Ax(bm,bn)
    integer,intent(in) ::am,an,bm,bn
    integer :: i,j,k
10  format(f10.5)
    emat=0.

    ! A*x
    do i = 1, bm
        do j = 1, bn
            Ax(i,j) = dot_product(A(i,:),x(:,j))
        enddo
    enddo

    ! Ax-b
    do i = 1, bm
        do j = 1, bn
            emat(i,j) = Ax(i,j) - b(i,j)
        enddo
    enddo

    print*, ' '
    print*, 'Error Matrix = '
    do i = 1,bm
        write(*,*) (emat(i,j), j=1,bn)
    enddo

end function resid
