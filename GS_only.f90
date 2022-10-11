
 program twoD_V_multigrid

!Variable declaration==========================================================

integer                                             :: i,j,k,l,g,level,lvl,N,M,iter
double precision                                    :: dx,dy,tol,o,start,finish
double precision, dimension(:), allocatable         :: x,y
double precision, dimension(:,:), allocatable       :: T,res


!User input====================================================================
print*, "Enter the grid dimension"
read*,M
N=M
print*,N,M
tol=1e-10


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
allocate(res(M+1,N+1))

T(:,:)=0.0
res(:,:)=0.0


print*, "What will happen"

call cpu_time(start)

do k=1,1000000
        do j=2,N-1
            do i=2,M-1
                T(i,j)= ((dy**2)*(T(i+1,j) + T(i-1,j)) + (dx**2)*(T(i,j+1) + T(i,j-1))+ 2*(dx*dy)**2*&
                        &((1-6*x(i)**2)*y(j)**2*(1-y(j)**2) + (1-6*y(j)**2)*x(i)**2*(1-x(i)**2)))/(2*(dx**2+dy**2))

            end do
        end do

    do j=2,N-1
        do i=2,M-1
            res(i,j)=-2*((1-6*x(i)**2)*y(j)**2*(1-y(j)**2) + (1-6*y(j)**2)*x(i)**2*(1-x(i)**2))&
                        &-(T(i+1,j)+T(i-1,j)-2*(T(i,j)))/(dx**2)-(T(i,j+1)+T(i,j-1)-2*(T(i,j)))/(dy**2)
        end do
    end do
    if(maxval(abs(res))<tol)then
       EXIT
    end if
end do
call cpu_time(finish)
print*,"Time=",finish-start,k

open(10,file='Poisso1_GS.dat')
write(10,*),"VARIABLES=""X"" ""Y"" ""T"""
write(10,*),"ZONE"
write(10,*),"i=",M
write(10,*),"j=",N
do i=1,N
    do j=1,M
        write(10,*),x(i), y(j), T(i,j)
    end do

end do


!print*,k,maxval(abs(res))
deallocate(T)

deallocate(x)
deallocate(y)
deallocate(res)

end program
