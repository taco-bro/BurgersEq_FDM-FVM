Program main
! For Demo
!  | 5 2 0 0 1 | x1 = 14       1
!  | 1 5 2 0 0 | x2 = 17       2
!  | 0 1 5 2 0 | x3 = 25       3
!  | 0 0 1 5 2 | x4 = 33       5
!  | 2 0 0 1 5 | x5 = 31       5
!
implicit none
integer,parameter      :: n=5
real*8, dimension(1:n) :: a,b,c,d
real*8, dimension(0:n+1) :: u
integer :: i

        do i=1,n
        a(i) = 1
        b(i) = 5
        c(i) = 2
        end do

        d = [14,17,25,33,31]
       call cyclic(a,b,c,c(n),a(1),d,u(1:n),n)
        write(*,*) (u(i),i=1,n)

END

!---------------------------------------------------------------------    
subroutine cyclic(a,b,c,alpha,beta,r,u_new,n)
    integer i,n
    real*8, dimension(n) ::  a,b,c,r,u_new,bb,u,z
    real*8               ::  alpha,beta,fact,gamma
      gamma=-b(1)
      bb(1)=b(1)-gamma
      bb(n)=b(n)-alpha*beta/gamma
      do i=2,n-1
        bb(i)=b(i)
      enddo
      call tridag(a,bb,c,r,u_new,n)
      u(1)=gamma
      u(n)=alpha
      do i=2,n-1
        u(i)=0.
      enddo
      call tridag(a,bb,c,u,z,n)
      fact=(u_new(1)+beta*u_new(n)/gamma)/(1.+z(1)+beta*z(n)/gamma)
      do i=1,n
        u_new(i)=u_new(i)-fact*z(i)
      enddo

end subroutine
!*********************************************************************
subroutine tridag(a,b,c,r,u,n)
    integer j,n
    real*8, dimension(n) ::  a,b,c,r,u,gam
    real*8               ::  alpha,beta,fact,gamma,bet
      bet=b(1)
      u(1)=r(1)/bet
      do j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        u(j)=(r(j)-a(j)*u(j-1))/bet
      enddo
      do j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
      enddo

end subroutine
!*********************************************************************
