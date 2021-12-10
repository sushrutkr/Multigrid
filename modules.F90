module global_variables
    integer :: i,j,k,nx,ny,nlev, maxiter, it_smoother, solver
    integer, parameter :: gs=1, line_sor=2, srj=3, multiG=2
    integer, parameter :: vcycle=1, wcycle=2, fcycle=3 
    real :: residual, tol, dx, dy, kw , lx, ly, w
    real :: start_time, end_time
    real, parameter :: pi = 3.1415927
    real, dimension(:,:), allocatable :: u, ukp1, uk,u_init,s
    real, dimension(:), allocatable :: x, y 

endmodule

module fnames
    integer, parameter :: input=12, postproc=13, logfile=14
endmodule 