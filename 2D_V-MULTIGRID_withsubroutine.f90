program twoD_V_multigrid

!Variable declaration==========================================================

integer                                             :: i,j,k,l,g,level,lvl,N,M,iter,mv
double precision                                    :: dx,dy,tol
double precision, dimension(:), allocatable         :: x,y
double precision, dimension(:,:), allocatable       :: T,Tcalc,res,er
double precision, dimension(:,:,:), allocatable     :: f,e

!User input====================================================================
print*, "Enter the grid dimension"
read*,M
N=M

10 print*, "Enter number of levels"
read*, level

if ((N-1)/2**level<=1)then
    print*,"Number of levels entered are not compatible."
    goto 10
end if

print*, "Enter the number of iterations at each node of the V cycles"
read*, iter
print*, "Enter the acceptable tolerance level"
read*, tol
print*, "Enter number of V-cycles"
read*, mv

!Assigning values==============================================================
dx=1.0/(N-1)
dy=1.0/(M-1)

!Grid generation===============================================================
allocate(x(N))
allocate(y(M))
x(:)=0.0
y(:)=0.0
do i=1,M-1
    x(i+1)=x(i)+dx
end do
do i=1,N-1
    y(i+1)=y(i)+dy
end do

!Allocation and initilisation of solution======================================
allocate(T(M+1,N+1))
allocate(Tcalc(M+1,N+1))
allocate(res(M+1,N+1))
allocate(er(M+1,N+1))
allocate(f((((N-1)/2))+1,((M-1)/2)+1,level))
allocate(e((((N-1)/2))+1,((M-1)/2)+1,level))

T(:,:)=0.0
res(:,:)=0.0




call cpu_time(start)
do l=1,mv


   !Appplication of scheme====================================================
   do k=1,iter
        do j=2,N-1
            do i=2,M-1
                T(i,j)= ((dy**2)*(T(i+1,j) + T(i-1,j)) + (dx**2)*(T(i,j+1) + T(i,j-1))+ 2*(dx*dy)**2*&
                        ((1-6*x(i)**2)*(y(j)**2)*(1-y(j)**2) + (1-6*y(j)**2)*(x(i)**2)*(1-x(i)**2)))/(2*(dx**2+dy**2))
            end do
        end do
    end do
    do j=2,N-1
        do i=2,M-1
            res(i,j)=-2*((1-6*x(i)**2)*(y(j)**2)*(1-y(j)**2) + (1-6*y(j)**2)*(x(i)**2)*(1-x(i)**2))&
                        -(T(i+1,j)+T(i-1,j)-2*(T(i,j)))/(dx**2)-(T(i,j+1)+T(i,j-1)-2*(T(i,j)))/(dy**2)
        end do
    end do



    !Restriction==========================================================================================================
    !=====================================================================================================================
    call restriction(lvl,level,M,N,iter,dx,res,e,f)


    !Prolongation=========================================================================================================
    !=====================================================================================================================
    call prolongation(lvl,level,M,N,iter,dx,e,f,er,T)


    !Re-iterate on finest mesh==============================================
    do k=1,iter
        do j=2,N-1
            do i=2,M-1
                T(i,j)= ((dy**2)*(T(i+1,j) + T(i-1,j)) + (dx**2)*(T(i,j+1) + T(i,j-1))+ 2*(dx*dy)**2*&
                        ((1-6*x(i)**2)*(y(j)**2)*(1-y(j)**2) + (1-6*y(j)**2)*(x(i)**2)*(1-x(i)**2)))/(2*(dx**2+dy**2))
            end do
        end do
    end do
    do j=2,N-1
        do i=2,M-1
            res(i,j)=-2*((1-6*x(i)**2)*(y(j)**2)*(1-y(j)**2) + (1-6*y(j)**2)*(x(i)**2)*(1-x(i)**2))&
                        -(T(i+1,j)+T(i-1,j)-2*(T(i,j)))/(dx**2)-(T(i,j+1)+T(i,j-1)-2*(T(i,j)))/(dy**2)
        end do
    end do



    if (maxval((abs(res)))<tol)then
        EXIT
    end if
end do
call cpu_time(finish)
print*,"Time=",finish-start

!Analytical Solution===============================================================



print*,l,maxval(abs(res))
do j=2,N-1
    do i=2,M-1
        Tcalc(i,j)= (x(i)**2)*(y(j)**2)*(1-(x(i))**2)*((y(j))**2-1)
        res(i,j)=(Tcalc(i,j)-T(i,j))*100/Tcalc(i,j)
    end do
end do

open(10,file='2.dat')
write(10,*),"VARIABLES=""X"" ""Y"" ""T"" ""Tanalytical"" ""Residual"""
write(10,*),"ZONE"
write(10,*),"i=",M
write(10,*),"j=",N
do i=1,N
    do j=1,M
        write(10,*),x(i), y(j), T(i,j),Tcalc(i,j),res(i,j)
    end do
end do

close(10)


deallocate(T)
deallocate(Tcalc)
deallocate(er)
deallocate(f)
deallocate(e)
deallocate(x)
deallocate(y)
deallocate(res)

end program

subroutine restriction(lvl,level,M,N,iter,dx,res,e,f)
    integer  ::i,j,k,g
    integer,intent(in)::M,N,level,iter
    integer,intent(inout)::lvl
    double precision,intent(in)::dx
    double precision, dimension(M+1,N+1),intent(inout) ::res
    double precision, dimension((N-1)/2+1,((M-1)/2)+1,level) ,intent(out)  :: f,e
    f(:,:,:)=0.0
    e(:,:,:)=0.0
    lvl=1
    do k=1,level
        do j=2,(N-1)/2**lvl
            do i=2,(M-1)/2**lvl
                f(i,j,lvl)=(res(2*i-2,2*j-1)+res(2*i,2*j-2)+res(2*i-2,2*j)+res(2*i,2*j)+2*&
                            (res(2*i-1,2*j-2)+res(2*i-1,2*j)+res(2*i-2,2*j-1)+res(2*i,2*j-1))+4*res(2*i-1,2*j-1))/16
            end do
        end do
        do g=1,iter                                             !/////Multiple inner iterations\\\\\
            do j=2,(N-1)/2**lvl
                do i=2,(M-1)/2**lvl
                    e(i,j,lvl)=0.25*(e(i+1,j,lvl)+e(i-1,j,lvl)+e(i,j+1,lvl)+e(i,j-1,lvl))-0.25*f(i,j,lvl)*(dx*2**lvl)**2
                end do
            end do
        end do
        res(:,:)=0.0
        do j=2,(N-1)/2**lvl
            do i=2,(M-1)/2**lvl
                res(i,j)=f(i,j,lvl)-(e(i+1,j,lvl)+e(i-1,j,lvl)+e(i,j+1,lvl)+e(i,j-1,lvl)-4*(e(i,j,lvl)))/(dx*2**lvl)**2
            end do
        end do
        lvl=lvl+1
    end do
end subroutine

subroutine prolongation(lvl,level,M,N,iter,dx,e,f,er,T)
    integer ::i,j,k,g
    integer, intent(in)::level,M,N,iter
    integer, intent(inout)::lvl
    double precision,intent(in)::dx
    double precision, dimension(M+1,N+1),intent(inout)::er
   double precision, dimension((N-1)/2+1,((M-1)/2)+1,level) ,intent(inout)  :: f,e
   double precision,dimension(M+1,N+1),intent(out)::T

    lvl=lvl-1
    do k=1,level
        er(:,:)=0.0
        do j=1,(N-1)/2**lvl+1
            do i=1,(M-1)/2**lvl+1
                er(2*i-1,2*j-1)=e(i,j,lvl)
            end do
        end do
        do j=1,(N-1)/2**lvl+1
            do i=1,(N-1)/2**lvl
                er(2*i,2*j-1)=0.5*(e(i,j,lvl)+e(i+1,j,lvl))
            end do
        end do
        do j=1,(N-1)/2**lvl
            do i=1,(N-1)/2**lvl+1
                er(2*i-1,2*j)=0.5*(e(i,j,lvl)+e(i,j+1,lvl))
            end do
        end do
        do j=1,(N-1)/2**lvl
            do i=1,(N-1)/2**lvl
                er(2*i,2*j)=0.25*(e(i,j,lvl)+e(i+1,j,lvl)+e(i,j+1,lvl)+e(i+1,j+1,lvl))
            end do
        end do
        if (lvl==1)then
            GOTO 30
        end if
        lvl=lvl-1
        do j=1,(N-1)/2**lvl+1
            do i=1,(M-1)/2**lvl+1
                e(i,j,lvl)=e(i,j,lvl)+er(i,j)
            end do
        end do
        do g=1,iter                     !/////Multiple inner iterations\\\\\
            do j=2,(N-1)/2**lvl
                do i=2,(M-1)/2**lvl
                    e(i,j,lvl)=0.25*(e(i+1,j,lvl)+e(i-1,j,lvl)+e(i,j+1,lvl)+e(i,j-1,lvl))-f(i,j,lvl)*0.25*(dx*2**lvl)**2
                end do
            end do
        end do
    end do

30  do j=1,N                            !/////Solution update at finest level\\\\\
        do i=1,M
            T(i,j)=T(i,j)+er(i,j)
        end do
    end do

end subroutine
