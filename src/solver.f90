module PoissonSolver

    use mpi
    implicit none

    contains

    subroutine eval_residual(nx, ny, dx, dy, x, y, btop, bbot, blft, bryt, temperature, residual, source_term)

        implicit none
        
        interface

            subroutine source_term(nx, ny, dx, dy, x, y, val)

                integer, intent(in)  :: nx, ny
                real*8,  intent(in)  :: dx, dy, x(nx, ny), y(nx, ny)
                real*8,  intent(out) :: val(nx, ny)

            end subroutine source_term

        end interface

        integer, intent(in)  :: nx, ny

        real*8,  intent(in)  :: dx, dy,              &
                                temperature(nx, ny), &
                                x(nx, ny),           &
                                y(nx, ny),           &
                                btop(nx),            &
                                bbot(nx),            &
                                blft(ny),            &
                                bryt(ny)

        real*8,  intent(out) :: residual(nx, ny)

        !-------------------------------------------------------------------------

        real*8 :: val(nx, ny)

        call source_term(nx, ny, dx, dy, x, y, val)

        residual(2:nx-1,2:ny-1) = - (2.0D0/dx**2 + 2.0D0/dy**2) * temperature(2:nx-1, 2:ny-1) +&
                                    (1.0D0/dx**2) * (temperature(1:nx-2,2:ny-1) + temperature(3:nx,2:ny-1)) +&
                                    (1.0D0/dy**2) * (temperature(2:nx-1,1:ny-2) + temperature(2:nx-1,3:ny)) +&
                                    val(2:nx-1, 2:ny-1)

        residual( 1,2:ny-1) = - (2.0D0/dx**2 + 2.0D0/dy**2) * temperature( 1, 2:ny-1) +&
                                (1.0D0/dx**2) * (          blft(2:ny-1) + temperature( 2,2:ny-1)) +&
                                (1.0D0/dy**2) * (temperature( 1,1:ny-2) + temperature( 1,3:ny))   +&
                                val( 1, 2:ny-1)

        residual(nx,2:ny-1) = - (2.0D0/dx**2 + 2.0D0/dy**2) * temperature(nx, 2:ny-1) +&
                                (1.0D0/dx**2) * (temperature(nx-1,2:ny-1) +         bryt(2:ny-1)) +&
                                (1.0D0/dy**2) * (temperature(nx,1:ny-2)   + temperature(nx,3:ny))   +&
                                val(nx, 2:ny-1)


        residual(2:nx-1, 1) = - (2.0D0/dx**2 + 2.0D0/dy**2) * temperature(2:nx-1, 1) +&
                                (1.0D0/dx**2) * (temperature(1:nx-2, 1) + temperature(3:nx, 1)) +&
                                (1.0D0/dy**2) * (          bbot(2:nx-1) + temperature(2:nx-1,2))   +&
                                val(2:nx-1,  1)


        residual(2:nx-1,ny) = - (2.0D0/dx**2 + 2.0D0/dy**2) * temperature(2:nx-1, ny) +&
                                (1.0D0/dx**2) * (temperature(1:nx-2,  ny) + temperature(3:nx, ny)) +&
                                (1.0D0/dy**2) * (temperature(2:nx-1,ny-1) + btop(2:nx-1))   +&
                                val(2:nx-1, ny)

        residual( 1, 1) = - (2.0D0/dx**2 + 2.0D0/dy**2) * temperature( 1, 1) +&
                            (1.0D0/dx**2) * (blft(1) + temperature(2,1))     +&
                            (1.0D0/dy**2) * (bbot(1) + temperature(1,2))     +&
                            val( 1, 1)

        residual(nx, 1) = - (2.0D0/dx**2 + 2.0D0/dy**2) * temperature(nx, 1) +&
                            (1.0D0/dx**2) * (temperature(nx-1,1) + bryt(1))  +&
                            (1.0D0/dy**2) * (temperature(  nx,2) + bbot(nx)) +&
                            val(nx, 1)


        residual( 1,ny) = - (2.0D0/dx**2 + 2.0D0/dy**2) * temperature( 1,ny)  +&
                            (1.0D0/dx**2) * (blft(ny) + temperature( 2,ny))   +&
                            (1.0D0/dy**2) * (btop(1)  + temperature( 1,ny-1)) +&
                            val( 1,ny)


        residual(nx,ny) = - (2.0D0/dx**2 + 2.0D0/dy**2) * temperature(nx,ny)  +&
                            (1.0D0/dx**2) * (bryt(ny) + temperature(nx-1,ny)) +&
                            (1.0D0/dy**2) * (btop(nx) + temperature(nx,ny-1)) +&
                            val(nx,ny)

    end subroutine eval_residual

    subroutine communicate(nx, ny, nzx, nzy, zonex, zoney, variable, btop, bbot, blft, bryt)

        implicit none

        integer,  intent(in)  :: nx, ny, nzx, nzy, zonex, zoney
        real*8,   intent(in)  :: variable(nx, ny)
        real*8,   intent(out) :: btop(nx), &
                                 bbot(nx), &
                                 blft(ny), &
                                 bryt(ny)

        integer :: this_proc, ierr, stat(MPI_STATUS_SIZE)

        call mpi_comm_rank(mpi_comm_world, this_proc, ierr)

        if (mod(zonex+zoney,2)==0) then
            if (zonex/=1) then
                call mpi_send(temperature( 1,:), ny, MPI_DOUBLE_PRECISION, this_proc - 1, 0, mpi_comm_world, ierr)
                call mpi_recv(blft,              ny, MPI_DOUBLE_PRECISION, this_proc - 1, 1, mpi_comm_world, stat, ierr)
            end if
            if (zonex/=nzx) then
                call mpi_send(temperature(nx,:), ny, MPI_DOUBLE_PRECISION, this_proc + 1, 10, mpi_comm_world, ierr)
                call mpi_recv(bryt,              ny, MPI_DOUBLE_PRECISION, this_proc + 1, 11, mpi_comm_world, stat, ierr)
            end if
            if (zoney/=1) then
                call mpi_send(temperature(:, 1), nx, MPI_DOUBLE_PRECISION, this_proc - nzx, 100, mpi_comm_world, ierr)
                call mpi_recv(bbot,              nx, MPI_DOUBLE_PRECISION, this_proc - nzx, 101, mpi_comm_world, stat, ierr)
            end if
            if (zoney/=nzy) then
                call mpi_send(temperature(:,ny), nx, MPI_DOUBLE_PRECISION, this_proc + nzx, 1000, mpi_comm_world, ierr)
                call mpi_recv(btop,              nx, MPI_DOUBLE_PRECISION, this_proc + nzx, 1001, mpi_comm_world, stat, ierr)
            end if
        else
            if (zonex/=nzx) then
                call mpi_recv(bryt,              ny, MPI_DOUBLE_PRECISION, this_proc + 1, 0, mpi_comm_world, stat, ierr)
                call mpi_send(temperature(nx,:), ny, MPI_DOUBLE_PRECISION, this_proc + 1, 1, mpi_comm_world, ierr)
            end if
            if (zonex/=1) then
                call mpi_recv(blft,              ny, MPI_DOUBLE_PRECISION, this_proc - 1, 10, mpi_comm_world, stat, ierr)
                call mpi_send(temperature( 1,:), ny, MPI_DOUBLE_PRECISION, this_proc - 1, 11, mpi_comm_world, ierr)
            end if
            if (zoney/=nzy) then
                call mpi_recv(btop,              nx, MPI_DOUBLE_PRECISION, this_proc + nzx, 100, mpi_comm_world, stat, ierr)
                call mpi_send(temperature(:,ny), nx, MPI_DOUBLE_PRECISION, this_proc + nzx, 101, mpi_comm_world, ierr)
            end if
            if (zoney/=1) then
                call mpi_recv(bbot,              nx, MPI_DOUBLE_PRECISION, this_proc - nzx, 1000, mpi_comm_world, stat, ierr)
                call mpi_send(temperature(:, 1), nx, MPI_DOUBLE_PRECISION, this_proc - nzx, 1001, mpi_comm_world, ierr)
            end if
        end if

    end subroutine communicate

    subroutine iterate(nx, ny, temperature)

        

    end subroutine iterate

end module PoissonSolver
