program numerical_project

    use global_variables

    call init_simulation()

    open(unit=logfile,file='log.dat',status='unknown')
    
    
    residual = 1
    s(:,:) = 0.0

    select case(solver)
    case(gs)
        write(logfile,*) 'nx = ', nx, ' ny = ', ny, ' Solver = Gauss Seidel'
        write(logfile,*) 'Iteration        Residual'
        do while(residual>tol .and. iter<maxiter)
            call gauss_seidel(nx,ny,dx,w,s,ukp1)
            call calc_residual(nx,ny,dx,dy,ukp1,residual)
            write(logfile,*) iter,residual
            uk = ukp1 
            iter = iter + 1
        end do
    case(multiG)
        call init_multigrid(nx,ny)
        write(logfile,*) 'nx = ', nx, ' ny = ', ny,' Solver = Multigrid, Smoother = ', it_smoother
        write(logfile,*) 'Iteration        Residual'
        do while(residual>tol .and. iter<maxiter)
            call multigrid(nx,ny,uk,nlev,it_smoother,w,ukp1)
            call calc_residual(nx,ny,dx,dy,ukp1,residual)
            write(logfile,*) iter,residual
            uk = ukp1 
            iter = iter + 1
        end do
        call end_multigrid()
    endselect

    u(:,:) = ukp1(:,:)

    close(logfile)

    call write_postproc(nx,ny,x,y,u)
    call deallocate_arrays()

end program numerical_project

subroutine init_simulation()
    use global_variables
    call read_input()
    call allocate_arrays()

    do i=1,nx
        x(i) = (i-1)*dx
    enddo

    do i = 1,ny
        y(i) = (i-1)*dy
    enddo

    ! Initialising Domain
    call random_generator(nx,ny,-1.0,1.0,u_init)
    open(15,file='random_number.dat',status='unknown')
    do i=1,ny 
        write(15,*) u_init(:,i)
    enddo
    ! close(15)
    ! open(17,file='random_number.dat',status='old')
    ! do i=1,ny 
    !     read(17,*) u_init(:,i)
    ! enddo
    ! close(17)

    u_init(:,:) = u_init(:,:)*0.0
    
    ! Setting Boundary Conditions
    u_init(:,1) = sin(kw*x(:))
    u_init(:,ny) = 0
    u_init(1,:) = sin(kw*y(:))
    u_init(nx,:) = 0

    write(*,*) sin(kw*x(2)),kw*x(2)
    u(:,:) = u_init(:,:)
    ukp1(:,:) = u(:,:)
    uk(:,:) = u(:,:)

    open(15,file='u_init.dat',status='unknown')
    do i=ny,1,-1 
        write(15,*) u_init(:,i)
    enddo
    close(15)

     
endsubroutine init_simulation

subroutine read_input()
    use global_variables
    use fnames

    open(input, file='input.dat', status='old')
    read(input,*) 
    read(input,*)
    read(input,*)
    read(input,*) nx, ny 
    read(input,*)
    read(input,*) kw 
    read(input,*)
    read(input,*)
    read(input,*) solver 
    read(input,*)
    read(input,*) it_smoother  
    read(input,*)
    read(input,*) w
    read(input,*)
    read(input,*) tol, maxiter
    read(input,*)
    read(input,*) nlev
    nx = nx + 1
    ny = ny + 1

    dx = 2*pi/(nx-1)
    dy = 2*pi/(ny-1)

    write(*,*) 'nx = ', nx, ' ny = ', ny, ' Solver = ', solver, &
    ' Smoother = ', it_smoother, ' w = ', w, ' tol = ', tol, ' maxiter = ', maxiter, ' nlev = ', nlev, &
    ' kw = ', kw, ' dx = ', dx, ' dy = ', dy


endsubroutine read_input

subroutine write_postproc(nx,ny,x,y,data)
    integer :: i,j,nx,ny
    real,dimension(nx,ny),intent(in) :: data 
    real,dimension(nx), intent(in) :: x 
    real,dimension(nx), intent(in) :: y

    open(12, file='data.dat', status='unknown')
    write(12,*) 'TITLE = "Post Processing Tecplot"'
    write(12,*) 'VARIABLES = "X", "Y", "U"'
    write(12,*) 'ZONE T="BIG ZONE", I=',nx,', J=',ny,', DATAPACKING=POINT'
    
    DO j=1,ny
        do i = 1,nx
            write(12,*) x(i), y(j), data(i,j)
        enddo
    enddo
    close(12)
endsubroutine

subroutine random_generator(nx,ny,a,b,x)
    implicit none
    integer, intent(in) :: nx,ny
    real,intent(in) :: a,b
    real,dimension(nx,ny),intent(out) :: x
    real,dimension(nx,ny) :: u
    call random_number(u)
    x = (b-a)*u + a
end subroutine random_generator


subroutine allocate_arrays()
    use global_variables

    allocate(u(nx,ny),ukp1(nx,ny),uk(nx,ny),u_init(nx,ny),s(nx,ny))
    allocate(x(nx),y(ny))

endsubroutine allocate_arrays

subroutine deallocate_arrays()
    use global_variables

    deallocate(u,ukp1,uk,u_init,s)
    deallocate(x,y)
endsubroutine deallocate_arrays

subroutine calc_residual(nx,ny,dx,dy,data,residual_out)
    integer :: i,j
    integer, intent(in) :: nx, ny 
    real, intent(out) :: residual_out
    real, intent(in) :: dx, dy
    real :: residual
    real, dimension(nx,ny), intent(in) :: data

    residual = 0

    ! L1 Norm
    do j=2,ny-1
        do i = 2,nx-1
            residual = residual + abs((1/(dx**2))*(data(i+1,j) -2*data(i,j) + data(i-1,j)) &
                                    + (1/(dy**2))*(data(i,j+1) -2*data(i,j) + data(i,j-1)))
        enddo
    enddo
    residual = residual/((nx-2)*(ny-2))

    ! L2 norm
    ! do j=2,ny-1
    !     do i = 2,nx-1
    !         residual = residual + ((1/(dx**2))*(data(i+1,j) -2*data(i,j) + data(i-1,j)) &
    !                              + (1/(dy**2))*(data(i,j+1) -2*data(i,j) + data(i,j-1)))**2
    !     enddo
    ! enddo

    ! residual = sqrt(residual)/((nx-2)*(ny-2))

    residual_out = residual
endsubroutine