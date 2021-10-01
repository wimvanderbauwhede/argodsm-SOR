program sor_main
use sor_params
use sor_routines
real, dimension(0:im+1,0:jm+1,0:km+1) :: p0
real, dimension(0:im+1,0:jm+1,0:km+1) :: p1
real, dimension(0:im+1,0:jm+1,0:km+1) :: rhs 
integer :: iter
integer, parameter :: niters = 1200
integer :: i,j,k


do i = 0,im+1
do j = 0,jm+1
do k = 0,km+1
    rhs(i,j,k) = 1.0
    p0(i,j,k) = 1.0
    end do
end do
end do

do iter = 1,niters
call sor (p0,p1,rhs)
p0=p1
end do

print *, p0(im/2,jm/2,km/2)

end program sor_main

