!------------------------------------------------------------------------------------------!
! Authors : Çağrı METIN                                                                    !
!                             Ege University                                               !
!           Faculty of Engineering, Mechanical Engineering Department                      !
!                                                                                          !
! This program is free software: you can redistribute it and/or modify                     !
! it under the terms of the GNU General Public License as published by                     !
! the Free Software Foundation, either version 3 of the License, or                        !
! (at your option) any later version.                                                      !
!                                                                                          !
! This program is distributed in the hope that it will be useful,                          !
! but WITHOUT ANY WARRANTY; without even the implied warranty of                           !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                            !
! GNU General Public License for more details.                                             !
!                                                                                          !
! You should have received a copy of the GNU General Public License                        !
! along with this program.  If not, see <http://www.gnu.org/licenses/>.                    !
!------------------------------------------------------------------------------------------!

program main
    implicit none
    integer, parameter :: nx = 41, ny = 41, nt = 700, nit = 50
    real(8), parameter :: xmin = 0.d0, ymin = 0.d0, xmax = 2.d0, ymax = 2.d0
    real(8)            :: dx, dy
    real(8), dimension(nx,ny) :: u, un, v, vn
    real(8), dimension(nx,ny) :: p, pn, b
    real(8), dimension(nx) :: x
    real(8), dimension(ny) :: y
    real(8), parameter :: rho = 1.d0, vis = 0.1d0, dt = 0.001d0
    integer :: i, j, k

    dx = (xmax - xmin) / (nx-1)
    dy = (ymax - ymin) / (ny-1)

    call  gridGen(x, y, dx, dy, nx, ny)

    u = 0.d0
    v = 0.d0
    p = 0.d0
    b = 0.d0

    do k = 1, nt   ! Loop over time
        un = u
        vn = v

        call buildUpB(b, rho, dt, u, v, dx, dy, nx, ny)
        call presPoisson(p, dx, dy, b, nx, ny, nit)

        call NavierStokes(u, un, v, vn, p, nx, ny, dx, dy, dt, rho, vis)
        call boundaryConditions(u, v, nx, ny)

    enddo

    call writeFile(u, v, x, y, nx, ny)


    contains


    subroutine buildUpB(b, rho, dt, u, v, dx, dy, nx, ny)
        integer, intent(in) :: nx, ny
        real(8), dimension(nx,ny), intent(inout) :: b
        real(8), dimension(nx,ny), intent(in) :: u, v
        real(8), intent(in) :: rho, dt, dx, dy
        integer :: i, j
        
        do j = 2, ny-1
            do i = 2, nx-1
               b(i,j) = rho*( (u(i+1,j)-u(i-1,j))/2.d0/dx + &
               & (v(i,j+1)-v(i,j-1))/2.d0/dy ) / dt - &
               & ( (u(i+1,j)-u(i-1,j))/2.d0/dx )**2 - &
               & 2.d0*(u(i,j+1)-u(i,j-1))/2.d0/dy * (v(i+1,j)-v(i,j-1))/2.d0/dx - &
               & ( (v(i,j+1)-v(i,j-1))/2.d0/dy)**2
            enddo
        enddo


    end subroutine buildUpB


    subroutine presPoisson(p, dx, dy, b, nx, ny, nit)
        integer, intent(in) :: nx, ny, nit
        real(8), dimension(nx, ny), intent(inout) :: p
        real(8), dimension(nx, ny) :: pn
        real(8), dimension(nx,ny), intent(in) :: b
        real(8), intent(in) :: dx, dy
        integer :: i, j, k

        pn = p
        do k = 1, nit
            pn = p
            do j = 2, ny-1
                do i = 2, nx-1 
                    p(i,j) = ((pn(i+1,j) + pn(i-1,j))*dy**2 + &
                           & (pn(i,j+1) + pn(i,j-1))*dx**2 - &
                           & b(i,j) * dx**2 * dy**2) / (dx**2 + dy**2) / 2.d0
                enddo
             enddo
        enddo

        do j = 1, ny
            p(nx,j) = p(nx-1,j)   ! dp/dy = 0 at x = 2
            p(1,j) = p(2,j)       ! dp/dy = 0 at x = 0
        enddo

        do i = 1, nx
            p(i,1) = p(i,2)       ! dp/dx = 0 at y = 0
            p(i,ny) = 0.d0        ! p = 0 at y = 2
        enddo

    end subroutine presPoisson


    subroutine gridGen(x, y, dx, dy, nx, ny)
        integer, intent(in) :: nx, ny
        real(8), dimension(nx), intent(inout) :: x
        real(8), dimension(ny), intent(inout) :: y
        real(8), intent(in) :: dx, dy
 
        x(1) = 0.d0
        do i = 2, nx
            x(i) = x(i-1) + dx
        end do

        y(1) = 0.d0
        do i = 2, ny
            y(i) = y(i-1) + dy
        end do
    end subroutine gridGen


    subroutine NavierStokes(u, un, v, vn, p, nx, ny, dx, dy, dt, rho, vis)
        integer, intent(in) :: nx, ny
        real(8), dimension(nx,ny), intent(inout) :: u, v
        real(8), dimension(nx,ny), intent(in) :: un, vn, p
        real(8), intent(in) :: dx, dy, dt
        real(8), intent(in) :: rho, vis

        integer :: i, j

        do j = 2, ny-1
            do i = 2, nx-1
                u(i,j) = un(i,j) - un(i,j)*dt/dx*( un(i,j)-un(i-1,j) ) - &
                       & vn(i,j)*dt/dy*( un(i,j)-un(i,j-1) ) - &
                       & 1.d0/rho*( p(i+1,j)-p(i-1,j) )*dt/2.d0/dx + &
                       & vis*dt/dx**2*( un(i+1,j)-2.d0*un(i,j)+un(i-1,j) ) + &
                       & vis*dt/dy**2*( un(i,j+1)-2.d0*un(i,j)+un(i,j-1) )

                v(i,j) = vn(i,j) - un(i,j)*dt/dx*( vn(i,j)-vn(i-1,j) ) - &
                       & vn(i,j)*dt/dy*( vn(i,j)-vn(i,j-1) ) - &
                       & 1.d0/rho*( p(i,j+1)-p(i,j-1) )*dt/2.d0/dy + &
                       & vis*dt/dx**2*( vn(i+1,j)-2.d0*vn(i,j)+vn(i-1,j) ) + &
                       & vis*dt/dy**2*( vn(i,j+1)-2.d0*vn(i,j)+vn(i,j-1) )
            enddo
        enddo 
    end subroutine NavierStokes


    subroutine boundaryConditions(u, v, nx, ny)
        integer, intent(in) :: nx, ny
        real(8), dimension(nx,ny), intent(inout) :: u, v
        integer :: i, j

        do i = 1, nx
            u(i,1) = 0.d0
            u(i,ny) = 1.d0

            v(i,1) = 0.d0
            v(i,ny) = 0.d0
        enddo

        do j = 1, ny
           u(1,j) = 0.d0
           u(nx,j) = 0.d0

           v(1,j) = 0.d0
           v(nx,j) = 0.d0
        enddo

    end subroutine boundaryConditions


    subroutine writeFile(u, v, x, y, nx, ny)
        integer, intent(in) :: nx, ny
        real(8), dimension(nx,ny), intent(in) :: u, v
        real(8), dimension(nx), intent(in) :: x
        real(8), dimension(ny), intent(in) :: y

        open(7, file='pressure.txt')
        do j = 1, ny
            do i = 1, nx
                write(7, *), x(i), ',', y(j), ',',  p(i,j)
            end do
        end do
        close(7)

        open(8, file='u_vel.txt')
        do j = 1, ny
            do i = 1, nx
                write(8, *), x(i), ',', y(j), ',',  u(i,j)
            end do
        end do
        close(8)

        open(9, file='v_vel.txt')
        do j = 1, ny
            do i = 1, nx
                write(9, *), x(i), ',', y(j), ',',  v(i,j)
            end do
        end do
        close(9)
    end subroutine writeFile

end program
