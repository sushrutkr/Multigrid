module multigrid_mod
    implicit none
    integer :: i,j,n
    integer :: mx1, my1, mx2, my2, smoother ! number of grid points in x and y direction
    real :: dxm, dym, w ! weight for smoother !!! using the same variable again might be a problem !!!
    real, parameter :: pi = 3.141592653589793238462643383279502884197
    real, dimension(:,:), allocatable :: f, gs_src, gsout 
    real, dimension(:,:,:), allocatable :: res, err
endmodule

subroutine init_multigrid(nx,ny,nlev)
    use multigrid_mod
    integer, intent(in) :: nx, ny, nlev
    allocate(res(nx,ny,nlev))
    allocate(err(nx,ny,nlev))
    allocate(f(nx,ny))
    allocate(gs_src(nx,ny),gsout(nx,ny))
end subroutine init_multigrid

subroutine multigrid(nx,ny,uk,nlev,solver,relax,ukp1)
    use multigrid_mod
    integer, intent(in) :: nx, ny, nlev, solver
    real, intent(in) :: relax
    real, dimension(nx,ny), intent(in) :: uk 
    real, dimension(nx,ny), intent(out) :: ukp1
    real, dimension(nx,ny) :: source

    w = relax
    smoother = solver
    mx1 = nx 
    my1 = ny
    dxm = 2*pi/(mx1-1)
    dym = 2*pi/(my1-1) 

    f(:,:) = uk(:,:)
    source(:,:) = 0.0
    res(:,:,:) = 0.0
    err(:,:,:) = 0.0

    call gauss_seidel(mx1,my1,dxm,w,source,f)
    ! open(15,file='gs1.dat',status='unknown')
    ! do i=my1,1,-1
    !     write(15,*) f(:,i)
    ! enddo
    ! close(15)

    call calc_r_multigrid()
    ! open(15,file='res1.dat',status='unknown')
    ! do i=my1,1,-1
    !     write(15,*) res(:,i,1)
    ! enddo
    ! close(15)

    ! Restriction
    do n=2,nlev
        call fine_to_coarse()
    enddo

    ! Prolongation
    do n=nlev-1,1,-1 
        call coarse_to_fine()
    enddo
  
    f(:,:) = f(:,:) - err(:,:,1)
    ukp1(:,:) = f(:,:)
endsubroutine multigrid

subroutine fine_to_coarse()
    use multigrid_mod

    mx2 = (mx1-1)/2 + 1
    my2 = (my1-1)/2 + 1
    dxm = 2*pi/(mx2-1)
    dym = 2*pi/(my2-1)
    call restriction()

    ! if (n==3) then 
    !     open(15,file='res2.dat',status='unknown')
    !     do i=my2,1,-1
    !         write(15,*) res(:,i,n)
    !     enddo
    !     close(15)
    ! endif
    !Solving Error Equation
    gsout(:,:) = 0.0
    gs_src(:,:) = 0.0
    gs_src(1:mx2,1:my2) = res(1:mx2,1:my2,n)
    do k = 1,n
        call gauss_seidel(mx2,my2,dxm,w,gs_src(1:mx2,1:my2),gsout(1:mx2,1:my2))
    end do
    err(1:mx2,1:my2,n) = gsout(:,:)
    call calc_r_error(mx2,my2)

    ! write(*,*) dxm

    ! if(n==2) then
    !     open(15,file='err2.dat',status='unknown')
    !     do i=my2,1,-1
    !         write(15,*) err(:,i,n)
    !     enddo
    !     close(15)
    !     open(15,file='res2.dat',status='unknown')
    !     do i=my2,1,-1
    !         write(15,*) res(:,i,n)
    !     enddo
    !     close(15)
    ! endif

    mx1 = mx2
    my1 = my2
    ! write(*,*) 'Restriction ',n
    ! write(*,*) mx1,my1
endsubroutine fine_to_coarse

subroutine coarse_to_fine()
    use multigrid_mod
    mx1 = (mx2-1)*2 + 1
    my1 = (my2-1)*2 + 1
    dxm = 2*pi/(mx1-1)
    dym = 2*pi/(my1-1)
    
    call prolongation()

    ! if (n==2) then
    !     write(*,*) 'Prolongation ',n
    !     write(*,*) res(2,2,n)
    !     write(*,*) err(1,2,n), err(3,2,n), err(2,1,n), err(2,3,n)
    ! endif
    ! if(n == 1) then
    !     open(15,file='err1.dat',status='unknown')
    !     do i=my1,1,-1
    !         write(15,*) err(:,i,n)
    !     enddo
    !     close(15)
    ! endif

    ! Check for current configuration - Target to 43(or42) iteration (Benchmark)
    gsout(:,:) = err(1:mx1,1:my1,n)
    gs_src(:,:) = 0.0
    gs_src(1:mx1,1:my1) = res(1:mx1,1:my1,n)
    do k = 1,1
        call gauss_seidel(mx1,my1,dxm,w,gs_src(1:mx1,1:my1),gsout(1:mx1,1:my1))
    end do
    err(1:mx1,1:my1,n) = gsout(:,:)

    ! if(n == 1) then
    !     open(15,file='errgsp.dat',status='unknown')
    !     do i=my1,1,-1
    !         write(15,*) err(:,i,n)
    !     enddo
    !     close(15)
    ! endif

    mx2 = mx1
    my2 = my1
    ! write(*,*) 'Prolongation',n
    ! write(*,*) mx1,my1

endsubroutine coarse_to_fine

subroutine restriction()
    use multigrid_mod
    integer :: i1, j1, incx, incy

    incx = 0
    incy = 0
    do j = 2, my2 - 1
        do i = 2, mx2 - 1
            i1 = 2*i - 1
            j1 = 2*j - 1
            res(i,j,n) = (res(i1,j1,n-1) &
                          + 0.5*(res(i1-1,j1,n-1) + res(i1+1,j1,n-1) + res(i1,j1+1,n-1) + res(i1,j1-1,n-1)) &
                          + 0.25*(res(i1-1,j1-1,n-1) + res(i1-1,j1+1,n-1) + res(i1+1,j1-1,n-1) + res(i1+1,j1+1,n-1)))/(4)
        enddo
    enddo
endsubroutine restriction

subroutine prolongation()
    use multigrid_mod
    integer :: i1, j1

    j1 = 1
    do j=2,mx1-1,2
        j1 = j1 + 1
        i1 = 1
        do i=2,mx1-1,2
            i1 = i1 + 1
            !element in fine is related to top right coarse element
            ! n+1 because value of n is interpolated from level (n+1) (eg level 1 obtained from level 2)
            err(i,j,n) = err(i,j,n) + (err(i1,j1,n+1) + err(i1-1,j1,n+1) + err(i1,j1-1,n+1) + err(i1-1,j1-1,n+1))/4
            ! res(i,j,n) = (res(i1,j1,n+1) + res(i1-1,j1,n+1) + res(i1,j1-1,n+1) + res(i1-1,j1-1,n+1))/4
        enddo
    enddo

    j1 = 0
    do j=2,my1-1,2
        j1 = j1 + 1
        i1 = 1
        do i=3,mx1-2,2
            i1 = i1 + 1
            !element in fine is related to top coarse element
            err(i,j,n) = err(i,j,n) + (err(i1,j1,n+1) + err(i1,j1+1,n+1))/2
            ! res(i,j,n) = (res(i1,j1,n+1) + res(i1,j1+1,n+1))/2
        enddo
    enddo

    j1 = 1
    do j=3,my1-2,2  
        j1 = j1 + 1
        i1 = 0
        do i=2,mx1-1,2
            i1 = i1 + 1
            !element in fine is related to left coarse element
            err(i,j,n) = err(i,j,n) + (err(i1,j1,n+1) + err(i1+1,j1,n+1))/2
            ! res(i,j,n) = (res(i1,j1,n+1) + res(i1+1,j1,n+1))/2
        enddo
    enddo

    do j=3,my1-2,2
        do i=3,mx1-2,2
            i1 = (i+1)/2
            j1 = (j+1)/2
            err(i,j,n) = err(i,j,n) + err(i1,j1,n+1)
            ! res(i,j,n) = res(i1,j1,n+1)
        enddo
    enddo

endsubroutine prolongation

subroutine gauss_seidel(nx,ny,del,w,source,data)
    integer :: i,j
    integer, intent(in) :: nx,ny
    real, intent(in) :: del,w
    real, dimension(nx,ny), intent(inout) :: data
    real, dimension(nx,ny), intent(in) :: source
    real, dimension(nx,ny) :: data_old

    data_old(:,:) = data(:,:)

    ! open(15,file='gss1.dat',status='unknown')
    ! do i=ny,1,-1
    !     write(15,*) data(:,i)
    ! enddo
    ! close(15)

    ! open(15,file='gsssource.dat',status='unknown')
    ! do i=ny,1,-1
    !     write(15,*) source(:,i)
    ! enddo
    ! close(15)
    ! i = 2
    ! j = 2
    ! if(ny==5) then
    !     write(*,*) 'ny=5',data(2,2), 0.25*(data(i+1,j)+data(i-1,j)+data(i,j+1)+data(i,j-1)),(del**2)*source(i,j)/4
    ! endif

    do j = 2, ny-1
        do i = 2, nx-1
            data(i,j) = 0.25*(data(i+1,j)+data(i-1,j)+data(i,j+1)+data(i,j-1)) - (del**2)*source(i,j)/4
        end do  
    end do

    data(:,:) = (1-w)*data_old(:,:) + w*data(:,:)
    ! open(15,file='gss2.dat',status='unknown')
    ! do i=ny,1,-1
    !     write(15,*) data(:,i)
    ! enddo
    ! close(15)

endsubroutine gauss_seidel

subroutine calc_r_multigrid()
    use multigrid_mod
    res = 0

    do j=2,my1-1
        do i = 2,mx1-1
            res(i,j,1) = ((1/(dxm**2))*(f(i+1,j) -2*f(i,j) + f(i-1,j)) &
                         + (1/(dym**2))*(f(i,j+1) -2*f(i,j) + f(i,j-1)))
        enddo
    enddo
endsubroutine

subroutine calc_r_error(nxr,nyr)
    use multigrid_mod
    integer, intent(in) :: nxr, nyr

    do j=2,nxr-1
        do i = 2,nyr-1
            res(i,j,n) =  -(((1/(dxm**2))*(err(i+1,j,n) -2*err(i,j,n) + err(i-1,j,n)) &
                         + (1/(dym**2))*(err(i,j+1,n) -2*err(i,j,n) + err(i,j-1,n)))  - res(i,j,n))
        enddo
    enddo
endsubroutine

subroutine end_multigrid()
    use multigrid_mod
    deallocate(res)
    deallocate(err)
    deallocate(f)
    deallocate(gs_src)
    deallocate(gsout)
end subroutine end_multigrid