module multigrid_mod
    implicit none
    integer :: i,j
    integer :: mx1, my1, mx2, my2, smoother ! number of grid points in x and y direction
    real :: dxm, dym, w ! weight for smoother !!! using the same variable again might be a problem !!!
    real, parameter :: pi = 3.141592653589793238462643383279502884197
    real, dimension(:,:), allocatable :: f
endmodule

subroutine init_multigrid(nx,ny)
    use multigrid_mod
    integer, intent(in) :: nx, ny
    mx1 = nx 
    my1 = ny
    mx2 = (mx1-1)/2 + 1
    my2 = (my1-1)/2 + 1
    ! allocate(res1(mx1,my1))
    ! allocate(res2(mx2,my2))
    ! allocate(err2(mx2,my2))
    ! allocate(err1(mx1,my1))

    allocate(f(nx,ny))
end subroutine init_multigrid

subroutine end_multigrid()
    use multigrid_mod
    deallocate(f)
    ! deallocate(res1,res2,err1,err2)
end subroutine end_multigrid

subroutine multigrid(nx,ny,uk,nlev,solver,relax,ukp1)
    use multigrid_mod
    integer, intent(in) :: nx, ny, nlev,solver
    real, intent(in) :: relax
    real, dimension(nx,ny), intent(in) :: uk 
    real, dimension(nx,ny), intent(out) :: ukp1
    real, dimension(nx,ny) :: source
    real, dimension(mx1,my1) :: res1, err1
    real, dimension(mx2,my2) :: res2, err2

    w = relax
    smoother = solver
    mx1 = nx 
    my1 = ny
    dxm = 2*pi/(mx1-1)
    dym = 2*pi/(my1-1) 

    f(:,:) = uk(:,:)
    source(:,:) = 0.0

    call gauss_seidel(mx1,my1,dxm,w,source,f)
    ! open(15,file='gs1.dat',status='unknown')
    ! do i=my1,1,-1
    !     write(15,*) f(:,i)
    ! enddo
    ! close(15)

    call calc_r_multigrid(res1)
    ! open(15,file='res1.dat',status='unknown')
    ! do i=my1,1,-1
    !     write(15,*) res1(:,i)
    ! enddo
    ! close(15)
    
    mx2 = (mx1-1)/2 + 1
    my2 = (my1-1)/2 + 1
    dxm = 2*pi/(mx2-1)
    dym = 2*pi/(my2-1)
    res2(:,:) = 0.0
    err2(:,:) = 0.0

    ! Restriction
    call restriction(res1,res2)
    ! open(15,file='res2.dat',status='unknown')
    ! do i=my2,1,-1
    !     write(15,*) res2(:,i)
    ! enddo
    ! close(15)

    !Solving Error Equation
    do k = 1,1
        call gauss_seidel(mx2,my2,dxm,w,res2,err2) ! Maybe multiple gauss seidel iterations?? - Could be done but thats a choice to make    
    end do
    ! open(15,file='err2.dat',status='unknown')
    ! do i=my2,1,-1
    !     write(15,*) err2(:,i)
    ! enddo
    ! close(15)

    ! Prolongation
    err1(:,:) = 0.0
    call prolongation(err2,err1)
    ! open(15,file='err1.dat',status='unknown')
    ! do i=my1,1,-1
    !     write(15,*) err1(:,i)
    ! enddo
    ! close(15)

    f(:,:) = f(:,:) - err1(:,:)
    ! open(15,file='fc.dat',status='unknown')
    ! do i=my1,1,-1
    !     write(15,*) f(:,i)
    ! enddo
    ! close(15)
    ukp1(:,:) = f(:,:)
endsubroutine multigrid

subroutine restriction(res11,res22)
    use multigrid_mod
    integer :: i1, j1, incx, incy
    real, dimension(mx1,my1), intent(in) :: res11
    real, dimension(mx2,my2), intent(inout) :: res22
    ! This is one method, other could be the increment one in notes - requires testing
    ! for restriction trapezoidal rule can be used 
    incx = 0
    incy = 0
    do j = 2, my2 - 1
        do i = 2, mx2 - 1
            i1 = 2*i - 1
            j1 = 2*j - 1
            res22(i,j) = (res11(i1,j1) &
                       + 0.5*(res11(i1-1,j1) + res11(i1+1,j1) + res11(i1,j1+1) + res11(i1,j1-1)) &
                       + 0.25*(res11(i1-1,j1-1) + res11(i1-1,j1+1) + res11(i1+1,j1-1) + res11(i1+1,j1+1)))/(4)
        enddo
    enddo
endsubroutine restriction


subroutine prolongation(err22,err11)
    use multigrid_mod
    integer :: i1, j1, incx, incy
    real, dimension(mx2,my2), intent(in) :: err22
    real, dimension(mx1,my1), intent(inout) :: err11

    
    j1 = 1
    do j=2,mx1-1,2
        j1 = j1 + 1
        i1 = 1
        do i=2,mx1-1,2
            i1 = i1 + 1
            !element in fine is related to top right coarse element
            err11(i,j) = (err22(i1,j1) + err22(i1-1,j1) + err22(i1,j1-1) + err22(i1-1,j1-1))/4
        enddo
    enddo

    j1 = 0
    do j=2,my1-1,2
        j1 = j1 + 1
        i1 = 1
        do i=3,mx1-2,2
            i1 = i1 + 1
            err11(i,j) = (err22(i1,j1) + err22(i1,j1+1))/2
        enddo
    enddo

    j1 = 1
    do j=3,my1-2,2  
        j1 = j1 + 1
        i1 = 0
        do i=2,mx1-1,2
            i1 = i1 + 1
            err11(i,j) = (err22(i1,j1) + err22(i1+1,j1))/2
        enddo
    enddo

    do j=3,my1-2,2
        do i=3,mx1-2,2
            i1 = (i+1)/2
            j1 = (j+1)/2
            err11(i,j) = err22(i1,j1)
        enddo
    enddo


endsubroutine prolongation

! subroutine itsolv
!     use multigrid_mod

!     ! write(*,*) smoother

!     select case(smoother)
!     case(1)
!         call gauss_seidel(mx1,my1,f,0.0)
!     case(2)
!         !call srj(nx,ny,uk,ukp1)
!     end select

! endsubroutine itsolv

subroutine gauss_seidel(nx,ny,del,w,source,data)
    integer :: i,j
    integer, intent(in) :: nx,ny
    real, intent(in) :: del,w
    real, dimension(nx,ny), intent(inout) :: data
    real, dimension(nx,ny), intent(in) :: source
    real, dimension(nx,ny) :: data_old

    data_old(:,:) = data(:,:)

    do j = 2, ny-1
        do i = 2, nx-1
            data(i,j) = 0.25*(data(i+1,j)+data(i-1,j)+data(i,j+1)+data(i,j-1)) - (del**2)*source(i,j)/4
        end do  
    end do

    data(:,:) = (1-w)*data_old(:,:) + w*data(:,:)

endsubroutine gauss_seidel

! subroutine srj(nx,ny,uk,ukp1)
!     use multigrid_mod

!     integer, intent(in) :: nx,ny
!     real, dimension(:,:), intent(in) :: uk

! endsubroutine srj

subroutine calc_r_multigrid(res11)
    use multigrid_mod
    real, dimension(mx1,my1), intent(inout) :: res11
    res11 = 0

    do j=2,my1-1
        do i = 2,mx1-1
            res11(i,j) = ((1/(dxm**2))*(f(i+1,j) -2*f(i,j) + f(i-1,j)) &
                       + (1/(dym**2))*(f(i,j+1) -2*f(i,j) + f(i,j-1)))
        enddo
    enddo
endsubroutine

subroutine allocate_multigrid()
    use multigrid_mod
endsubroutine

subroutine deallocate_multigrid()
    use multigrid_mod
endsubroutine