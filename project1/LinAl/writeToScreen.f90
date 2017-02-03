
!AMS213A-Spritzer
!Homework 1, write to screen subroutine.

module writeScreen
    implicit none

contains

    !***********************************************************
    !***** Write Data to Screen  *******************************
    subroutine writeToScreen(am,an,bm,bn,A,b)

        implicit none
        integer :: am,an,bm,bn,i,j,k
        real,intent(in) :: A(am,an)
        real,intent(in) :: b(bm,bn)
10      format(7f5.2)
        print*, ' '
        print*, 'Matrix = '
        do i=1,am
            write(*,10) (A(i,j), j=1,am)
        enddo

        print*, ' '
        print*, 'Vector = '
        do k=1,bm
            write(*,10) (b(k,j), j=1,bn)
        enddo

    end subroutine writeToScreen

end module writeScreen
