program oneD_V_multigrid

!Variable declaration==========================================================
integer                                         :: i,j,k,l,level,lvl,N,iter,mv
double precision                                :: dx,tol
double precision, dimension(:), allocatable     :: x
double precision, dimension(:), allocatable     :: U,res, er
double precision, dimension(:,:), allocatable   :: f,e

!User input====================================================================
print*, "Enter the number of fine grid points"
read*, N
10 print*, "Enter the number of levels"
read*, level

if ((N-1)/2**level<=1)then
    print*,"Number of levels entered are not compatible."
    goto 10
end if


20 print*, "Enter the number of iterations at each node of the V cycles"
read*, iter
print*, "Enter the acceptable tolerance level"
read*, tol
print*, "Enter number of V-cycles"
read*, mv

!Assigning values==============================================================
dx=1.0/(N-1)
pi=3.14159265

!Grid generation===============================================================
allocate(x(N))
x(:)=0.0
do i=1,N-1
    x(i+1)=x(i)+dx
end do

!Allocation and initialisation of solution=====================================
allocate(U(N+1))
allocate(res(N+1))
allocate(er(N+1))
allocate(f((((N-1)/2))+1,5))
allocate(e((((N-1)/2))+1,5))

U(:)=0.0
res(:)=0.0

do k=1,mv
!Application of scheme=====================================================
    do l=1,iter                      !/////Multiple inner iterations/////
        do i=2,64
            U(i)=0.5*(U(i+1)+U(i-1))-((dx)**2)*0.25*(sin(pi*x(i))+sin(16*pi*x(i)))
        end do
        do i=2,64
            res(i)=0.5*(sin(pi*x(i))+sin(16*pi*x(i)))-(U(i+1)+U(i-1)-2*U(i))/(dx**2)
        end do
    end do


    !Restrcition===============================================================
    f(:,:)=0.0
    e(:,:)=0.0
    lvl=1
    do j=1,level
        do i=2,(N-1)/2**lvl
            f(i,lvl)=0.25*(res(2*i-2)+2*res(2*i-1)+res(2*i))
            e(i,lvl)=0.5*(e(i+1,lvl)+e(i-1,lvl))-f(i,lvl)*0.5*(dx*2**lvl)**2
        end do
        res(:)=0.0
        do l=1,iter                     !/////Multiple inner iterations/////
            do i=2,(N-1)/2**lvl
                res(i)=f(i,lvl)-(e(i+1,lvl)+e(i-1,lvl)-2*e(i,lvl))/(dx*2**lvl)**2
            end do
        end do
        lvl=lvl+1
    end do

    !Prolongation==============================================================
    lvl=lvl-1
    do j=1,level
        er(:)=0.0
        do i=1,(N-1)/2**lvl+1
            er(2*i-1)=e(i,lvl)
        end do
        do i=1,(N-1)/2**lvl
            er(2*i)=0.5*(e(i,lvl)+e(i+1,lvl))
        end do
        if (lvl==1)then
            GOTO 30
        end if
        lvl=lvl-1
        do i=1,(N-1)/2**lvl+1
            e(i,lvl)=e(i,lvl)+er(i)
        end do
        do l=1,iter                     !/////Multiple inner iterations/////
            do i=2,(N-1)/2**lvl
                e(i,lvl)=0.5*(e(i+1,lvl)+e(i-1,lvl))-f(i,lvl)*0.5*(dx*2**lvl)**2
            end do
        end do
    end do
30  do i=1,65
        U(i)=U(i)+er(i)                 !/////Solution updation at finest level/////
    end do


    !Re-iterate on finest mesh==============================================

    do l=1,iter                      !/////Multiple inner iterations/////
        do i=2,64
            U(i)=0.5*(U(i+1)+U(i-1))-((dx)**2)*0.25*(sin(pi*x(i))+sin(16*pi*x(i)))
        end do
        do i=2,64
            res(i)=0.5*(sin(pi*x(i))+sin(16*pi*x(i)))-(U(i+1)+U(i-1)-2*U(i))/(dx**2)
        end do
    end do

    if (maxval(res)<tol)then
        EXIT
    end if


end do
print*,k,maxval(res)
deallocate(U)
deallocate(er)
deallocate(f)
deallocate(e)
deallocate(x)

if(maxval(res)>tol)then
    print*, "Solution didn't converge. Please use different set of parameters for the V-cycle"
    deallocate(res)
    GOTO 20
end if
deallocate(res)




end program oneD_V_multigrid
