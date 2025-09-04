program Main
    use Utils
    character(100) :: filename
    real(dp), allocatable, dimension(:,:) :: x,y,l2_res
    real(dp), dimension(2,2) :: b_forces
    real(dp), allocatable, dimension(:,:,:) :: q_mesh
    integer :: imax,jmax,n
    real(dp) :: M,p,T,v2,v4,CFL

    CFL = 0.1_dp
    M = 0.3_dp
    p = 101325.0_dp
    T = 295.0_dp
    
    v2 = 0.3
    v4 = 0.02
    n = 15000
    allocate(l2_res(n,4))

    filename = "ortho_mesh1.dat"
    call read_tecplot_mesh(filename,x,y)
    imax = size(x, dim=1)
    jmax = size(y, dim=2)
    allocate(q_mesh(-1:imax+1,-1:jmax+1,4))
    call init_mesh(q_mesh,M,p,T,x,y)
    call runge_kutta(q_mesh,x,y,CFL,M,p,T,n,v2,v4,l2_res)
    call write_py_sparse_q(x,y,q_mesh(1:imax-1,1:jmax-1,:),imax,jmax,"data")
    call write_res("res.dat",l2_res)
    call bump_forces(q_mesh,b_forces,x,y)
    write(*,*) b_forces(1,:)
    write(*,*) b_forces(2,:)
    deallocate(q_mesh)
    
end program Main