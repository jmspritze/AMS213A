!***********************************************************************
!************************** Program 1 **********************************
! Cholesky solution of the LS problem, use 'make' to compile program,
! then use ./Chol_LS to run the program. The inputs to the program are
! the file 'atkins.dat' and user input to determine the order of
! polynomial to fit. The output of the program prints the x solution
! vector, the residual error, rms of the calculated vs original data,
 ! and prints the fitted data to 'output.txt'.
!***********************************************************************
program driver1

    use LinAl

    implicit none
    integer :: msize,nsize,bmsize,bnsize
    integer :: n,io,i,j,k, order
    real :: u,v,w
    logical :: iSing = .false.

    !************* Arrays ***********************************
    real, dimension(:,:), allocatable :: A,At,data
    real, dimension(:), allocatable :: b,x,y
    real, dimension(:,:), allocatable :: L,emat
    !**********************************************************
    msize = 21
    nsize = 2

    !**** ALLOCATE ALL MATRICIES ********************************
    ALLOCATE(data(msize, nsize))
    ALLOCATE(L(order,order))
    ALLOCATE(At(order,msize))
    ALLOCATE(x(msize))
    ALLOCATE(b(msize))
    ALLOCATE(y(msize))
    ALLOCATE(emat(msize,nsize))
    !************************************************************
    x=0.
    y=0.

    open(10,file='atkinson.dat',form='formatted',status = 'old')
    do i=1,msize
        read(10,*) (data(i,j),j=1,nsize)
    enddo
    close(10)

    do i = 1,msize
        b(i) = data(i,2)
    enddo

    print*," "
    print*, "Enter the order of the polynomial to fit the data: "
    read(*,*) order

    !order + 1 add to order for the 1 in the vandermond matrix
    order = order+1
    call vanderm(Data, A, order, msize, nsize)

    L = matmul(transpose(A),A)
    b = matmul(transpose(A),b)

   At = L
    call Chol(order,L,iSing)

    if(iSing .eqv. .true.) then
        write(*,*) "Matrix is singular"
        stop
    end if

    call cholbacksub(L,x,b,y,order)
    print*," "
    print*, 'Solution X (1, x, x^2, ...x^n) = '
    call writeToScreen(x,order,1) ! Write Data to screen

    call resid(At,x,b,emat,order) ! error check on Ax = b
    call buildOut(data,x,msize,order)

    DEALLOCATE(Data,A,L,b,y,x)
end program driver1

