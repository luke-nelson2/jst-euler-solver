module Utils
    
    use iso_fortran_env, only: dp => real64
    use ieee_arithmetic
    implicit none
    
    real(dp) :: gamma = 1.4
    real(dp) :: R = 287.052874
    real(dp) :: pi = 3.1415926
    
contains
    
    subroutine mesh_area(x,y,area)
        
        real(dp), allocatable, intent(in), dimension(:,:) :: x,y
        real(dp), allocatable, intent(out), dimension(:,:) :: area
        integer :: imax, jmax
        
        imax = size(x, dim=1)
        jmax = size(y, dim=2)
        
        area = 0.5_dp * ((x(2:,2:)-x(1:imax-1,1:jmax-1))*(y(1:imax-1,2:)&
        -y(2:,1:jmax-1))-(y(2:,2:)-y(1:imax-1,1:jmax-1))*(x(1:imax-1,2:)-x(2:,1:jmax-1)))
    
    end subroutine mesh_area

    
    subroutine mesh_delta(x,y,dxh,dxv,dyh,dyv)
        
        real(dp), allocatable, dimension(:,:), intent(in) :: x,y
        real(dp), allocatable, dimension(:,:), intent(out) :: dxh,dxv,dyh,dyv
        integer :: imax, jmax
        
        imax = size(x, dim=1)
        jmax = size(y, dim=2)

        allocate(dxh(imax-1,0:jmax-1))
        allocate(dxv(0:imax-1,jmax-1))
        allocate(dyh(imax-1,0:jmax-1))
        allocate(dyv(0:imax-1,jmax-1))
        
        dxh = x(2:,:) - x(:imax-1,:)
        dxv = x(:,2:) - x(:,:jmax-1)
        dyh(:,:) = y(2:,:) - y(:imax-1,:)
        dyv = y(:,2:) - y(:,:jmax-1)
    
    end subroutine mesh_delta

    
    subroutine cell_flux(q,f,g)
        
        real(dp), dimension(4), intent(in) :: q
        real(dp), dimension(4), intent(out) :: f, g
        real(dp) :: p
        
        p = (gamma-1)*(q(4)-0.5_dp*(q(2)**2+q(3)**2)/q(1))
        
        f(1) = q(2)
        f(2) = q(2)**2/q(1)+p
        f(3) = q(2)*q(3)/q(1)
        f(4) = q(2)/q(1)*(q(4)+p)

        g(1) = q(3)
        g(2) = q(3)*q(2)/q(1)
        g(3) = q(3)**2/q(1) + p
        g(4) = q(3)/q(1)*(q(4)+p)

    end subroutine cell_flux

    
    subroutine mesh_face_lengths(x,y,lenh,lenv)
        
        real(dp), allocatable, dimension(:,:), intent(in) :: x,y
        real(dp), allocatable, dimension(:,:), intent(out) :: lenh,lenv
        real(dp), allocatable, dimension(:,:) :: dxh,dxv,dyh,dyv
        integer :: icmax, jcmax

        icmax = size(x, dim=1) - 1
        jcmax = size(x, dim=2) - 1

        call mesh_delta(x,y,dxh,dxv,dyh,dyv)

        allocate(lenh(icmax,0:jcmax))
        allocate(lenv(0:icmax,jcmax))

        lenh = sqrt(dxh**2+dyh**2)
        
        lenv = sqrt(dxv**2+dyv**2)

    end subroutine mesh_face_lengths

    
    subroutine mesh_normals(x,y,nh,nv)
        
        real(dp), allocatable, dimension(:,:), intent(in) :: x,y
        real(dp), allocatable, dimension(:,:,:), intent(out) :: nh,nv
        real(dp), allocatable, dimension(:,:) :: dxh,dxv,dyh,dyv,lenh,lenv
        integer :: icmax, jcmax
        
        icmax = size(x, dim=1) - 1
        jcmax = size(x, dim=2) - 1
        allocate(nh(icmax,0:jcmax,2))
        allocate(nv(0:icmax,jcmax,2))

        call mesh_delta(x,y,dxh,dxv,dyh,dyv)
        call mesh_face_lengths(x,y,lenh,lenv)

        ! Consensus here is that norms on horizontals point up and
        ! on verticals they point right
        nh(:,:,1) = -dyh(:,:)/lenh(:,:)
        nh(:,:,2) = dxh(:,:)/lenh(:,:)
        nv(:,:,1) = dyv(:,:)/lenv(:,:)
        nv(:,:,2) = -dxv(:,:)/lenv(:,:)

    end subroutine mesh_normals

    
    subroutine mesh_pressure(q_mesh, p_mesh)
        
        real(dp), allocatable, dimension(:,:,:), intent(in) :: q_mesh
        real(dp), allocatable, dimension(:,:), intent(out) :: p_mesh
        integer :: icmax, jcmax
        icmax = size(q_mesh, dim=1) - 4
        jcmax = size(q_mesh, dim=2) - 4

        allocate(p_mesh(-1:icmax+2,-1:jcmax+2))
        
        p_mesh(:,:) = (gamma-1.0_dp)*(q_mesh(:,:,4)-0.5_dp*(q_mesh(:,:,2)**2+q_mesh(:,:,3)**2)/q_mesh(:,:,1))
    
    end subroutine mesh_pressure

    
    subroutine mesh_residuals(q_mesh,dxh,dxv,dyh,dyv,r_mesh)
        
        real(dp), allocatable, dimension(:,:,:), intent(in) :: q_mesh
        real(dp), allocatable, dimension(:,:), intent(in) :: dxh,dxv,dyh,dyv
        real(dp), allocatable, dimension(:,:,:), intent(out) :: r_mesh

        real(dp), allocatable, dimension(:,:,:) :: f_mesh, g_mesh
        real(dp), allocatable, dimension(:) :: temp_f, temp_g,temp_q
        integer :: icmax, jcmax, i, j,k

        ! Perimeter of ghost cells of width 2
        icmax = size(q_mesh, dim=1) - 4
        jcmax = size(q_mesh, dim=2) - 4

        allocate(f_mesh(0:icmax+1,0:jcmax+1,4))
        allocate(g_mesh(0:icmax+1,0:jcmax+1,4))
        allocate(r_mesh(icmax,jcmax,4))
        allocate(temp_q(4))
        allocate(temp_f(4))
        allocate(temp_g(4))

        do j=0, jcmax+1
            
            do i=0, icmax+1
                
                temp_q = q_mesh(i,j,:)
                call cell_flux(temp_q,temp_f,temp_g)
                f_mesh(i,j,:) = temp_f
                g_mesh(i,j,:) = temp_g
            
            end do
        
        end do
        
        do concurrent (k=1:4)
            
            r_mesh(:,:,k) = 0.5_dp*( & 
            (f_mesh(1:icmax,1:jcmax,k)+f_mesh(1:icmax,2:jcmax+1,k))*-dyh(1:icmax,1:jcmax) &
            - (g_mesh(1:icmax,1:jcmax,k)+g_mesh(1:icmax,2:jcmax+1,k))*-dxh(1:icmax,1:jcmax) & ! NORTH
            + (f_mesh(1:icmax,1:jcmax,k)+f_mesh(0:icmax-1,1:jcmax,k))*-dyv(0:icmax-1,1:jcmax) &
            - (g_mesh(1:icmax,1:jcmax,k)+g_mesh(0:icmax-1,1:jcmax,k))*-dxv(0:icmax-1,1:jcmax) & ! WEST
            + (f_mesh(1:icmax,1:jcmax,k)+f_mesh(2:icmax+1,1:jcmax,k))*dyv(1:icmax,1:jcmax) & 
            - (g_mesh(1:icmax,1:jcmax,k)+g_mesh(2:icmax+1,1:jcmax,k))*dxv(1:icmax,1:jcmax) & ! EAST
            + (f_mesh(1:icmax,1:jcmax,k)+f_mesh(1:icmax,0:jcmax-1,k))*dyh(1:icmax,0:jcmax-1) &
            - (g_mesh(1:icmax,1:jcmax,k)+g_mesh(1:icmax,0:jcmax-1,k))*dxh(1:icmax,0:jcmax-1)) ! SOUTH
        
        end do
    
    end subroutine mesh_residuals

    
    ! Assuming alpha = 0
    ! We'll say the far field is a ghost cell
    subroutine inlet_boundary(q_mesh,M,p,T)
        
        real(dp), allocatable, dimension(:,:,:), intent(inout) :: q_mesh
        real(dp), intent(in) :: M,p,T
        real(dp) :: c_ff, v_ff
        real(dp), allocatable, dimension(:) :: p0_ff, p_i2, r1, r2, u1,c1,m1,rho1,p1,e1
        integer :: icmax,jcmax

       
        icmax = size(q_mesh, dim=1) - 4
        jcmax = size(q_mesh, dim=2) - 4

        allocate(p0_ff(-1:jcmax+2))
        allocate(r1(-1:jcmax+2))
        
        p0_ff(:) = p*(1.0_dp+0.5_dp*(gamma-1.0_dp)*M**2)**(gamma/(gamma-1.0_dp))
        c_ff = sqrt(gamma*R*T)
        v_ff = M*c_ff
        r1(:) = v_ff + 2.0_dp*c_ff/(gamma-1.0_dp)
        
        p_i2 = (gamma-1.0_dp)*(q_mesh(1,:,4)-0.5_dp*(q_mesh(1,:,2)**2+q_mesh(1,:,3)**2)/q_mesh(1,:,1))
        
        r2 = sqrt(q_mesh(1,:,2)**2 + q_mesh(1,:,3)**2)/q_mesh(1,:,1)-2.0_dp*sqrt(gamma*p_i2/q_mesh(1,:,1))/(gamma-1.0_dp)
        
        u1 = 0.5_dp*(r1+r2)
        c1 = 0.25_dp*(gamma-1)*(r1-r2)
        m1 = u1/c1
        
        p1 = p0_ff/(1.0_dp+0.5_dp*(gamma-1)*m1**2)**(gamma/(gamma-1))
        rho1 = gamma*p1/c1**2
        e1 = p1/(rho1*(gamma-1)) + 0.5_dp*u1**2
        
        q_mesh(0,:,1) = rho1
        q_mesh(0,:,2) = rho1*u1*0.866_dp
        
        q_mesh(0,:,3) = rho1*u1*0.5_dp
        q_mesh(0,:,4) = rho1*e1

        q_mesh(-1,:,:) = q_mesh(0,:,:)

    end subroutine inlet_boundary


    subroutine wall_boundary(q_mesh,nh)

        real(dp), allocatable, dimension(:,:,:), intent(inout) :: q_mesh
        real(dp), allocatable, dimension(:,:,:), intent(in) :: nh
        real(dp), allocatable, dimension(:,:) :: flow_norm_top1,flow_norm_bot1,flow_norm_top2,flow_norm_bot2
        integer :: icmax, jcmax, i

        icmax = size(q_mesh, dim=1)-4
        jcmax = size(q_mesh, dim=2)-4

        allocate(flow_norm_top1(icmax,2))
        allocate(flow_norm_bot1(icmax,2))
        allocate(flow_norm_top2(icmax,2))
        allocate(flow_norm_bot2(icmax,2))

        do concurrent (i=1:icmax)
            
            flow_norm_top1(i,:) = dot_product(q_mesh(i,jcmax,2:3),nh(i,jcmax,:))*nh(i,jcmax,:)
            flow_norm_bot1(i,:) = dot_product(q_mesh(i,1,2:3),-nh(i,0,:))*(-nh(i,0,:))

            flow_norm_top2(i,:) = dot_product(q_mesh(i,jcmax-1,2:3),nh(i,jcmax,:))*nh(i,jcmax,:)
            flow_norm_bot2(i,:) = dot_product(q_mesh(i,2,2:3),-nh(i,0,:))*(-nh(i,0,:))
        
        end do

        ! Top boundary
        q_mesh(1:icmax,jcmax+1,1) = q_mesh(1:icmax,jcmax,1)
        q_mesh(1:icmax,jcmax+1,2:3) = q_mesh(1:icmax,jcmax,2:3) - 2.0_dp* flow_norm_top1(:,:)
        q_mesh(1:icmax,jcmax+1,4) = q_mesh(1:icmax,jcmax,4)
        ! Outer ghost cell
        q_mesh(1:icmax,jcmax+2,1) = q_mesh(1:icmax,jcmax-1,1)
        q_mesh(1:icmax,jcmax+2,2:3) = q_mesh(1:icmax,jcmax-1,2:3) -2.0_dp* flow_norm_top2(:,:)
        q_mesh(1:icmax,jcmax+2,4) = q_mesh(1:icmax,jcmax-1,4)
        
        ! Bot boundary
        q_mesh(1:icmax,0,1) = q_mesh(1:icmax,1,1)
        q_mesh(1:icmax,0,2:3) = q_mesh(1:icmax,1,2:3) -2.0_dp* flow_norm_bot1(:,:)
        q_mesh(1:icmax,0,4) = q_mesh(1:icmax,1,4)
        ! Outer ghost cell
        q_mesh(1:icmax,-1,1) = q_mesh(1:icmax,2,1)
        q_mesh(1:icmax,-1,2:3) = q_mesh(1:icmax,2,2:3) -2.0_dp* flow_norm_bot2(:,:)
        q_mesh(1:icmax,-1,4) = q_mesh(1:icmax,2,4)

    end subroutine wall_boundary


    subroutine outlet_boundary(q_mesh, p_ff)

        real(dp), allocatable, dimension(:,:,:), intent(inout) :: q_mesh
        real(dp), allocatable, dimension(:) :: p_imax,p_imaxn, c_imax,c_imaxn, v_imax, r1_imax, vt_imax, e_imax,u_imax,rho_imax
        real(dp), intent(in) :: p_ff
        real(dp) :: m_avg,c_avg,v_avg
        integer :: icmax, jcmax

        icmax = size(q_mesh,dim=1) - 4
        jcmax = size(q_mesh,dim=2) - 4

        ! check if subsonic
        p_imax = (gamma-1.0_dp)*(q_mesh(icmax+1,:,4)-0.5_dp*(q_mesh(icmax+1,:,2)**2+q_mesh(icmax+1,:,3)**2)/q_mesh(icmax+1,:,1))
        p_imaxn = (gamma-1.0_dp)*(q_mesh(icmax,:,4)-0.5_dp*(q_mesh(icmax,:,2)**2+q_mesh(icmax,:,3)**2)/q_mesh(icmax,:,1))
        c_imax = sqrt(gamma*p_imax/q_mesh(icmax+1,:,1)) !
        c_imaxn = sqrt(gamma*p_imaxn/q_mesh(icmax,:,1)) !
        v_imax = sqrt(q_mesh(icmax+1,:,2)**2 + q_mesh(icmax+1,:,3)**2)/q_mesh(icmax+1,:,1) !
        c_avg = sum(c_imax)/size(c_imax)
        v_avg = sum(v_imax)/size(v_imax)
        m_avg = v_avg/c_avg
        
        if (m_avg >= 1.0_dp) then
            p_imax = p_imaxn
        else 
            p_imax(:) = p_ff
            ! write(*,*) 'ello'
        end if
        
        rho_imax = q_mesh(icmax,:,1) * (p_imax/p_imaxn)**(1.0_dp/gamma) !
        r1_imax = sqrt(q_mesh(icmax,:,2)**2+q_mesh(icmax,:,3)**2)/q_mesh(icmax,:,1) + 2*c_imaxn/(gamma-1.0_dp)

        c_imax = sqrt(gamma*p_imax/rho_imax)
        vt_imax = r1_imax - 2*c_imax/(gamma-1.0_dp)

        v_imax = q_mesh(icmax,:,3)/q_mesh(icmax,:,1)
        u_imax = sqrt(vt_imax**2-v_imax**2)
        e_imax = p_imax/(gamma-1.0_dp)/rho_imax + 0.5_dp*(vt_imax**2)
        
        q_mesh(icmax+1,:,1) = rho_imax
        q_mesh(icmax+1,:,2) = rho_imax*u_imax
        q_mesh(icmax+1,:,3) = rho_imax*v_imax
        
        q_mesh(icmax+1,:,4) = rho_imax*e_imax

        q_mesh(icmax+2,:,:) = q_mesh(icmax+1,:,:)
        
    end subroutine outlet_boundary


    subroutine mesh_switches(q_mesh,s2_zeta,s2_eta,v2)

        real(dp), allocatable, dimension(:,:,:), intent(in) :: q_mesh
        real(dp), allocatable, dimension(:,:), intent(out) :: s2_zeta, s2_eta
        real(dp), intent(in) :: v2
        real(dp), allocatable, dimension(:,:) :: p_mesh
        real(dp) :: epsilon
        integer :: icmax,jcmax

        icmax = size(q_mesh,dim=1)-4
        jcmax = size(q_mesh,dim=2)-4
        epsilon = 0.00001

        call mesh_pressure(q_mesh,p_mesh)

        allocate(s2_zeta(0:icmax+1,-1:jcmax+2))
        allocate(s2_eta(-1:icmax+2,0:jcmax+1))

        ! s2_zeta of size i: 102 > 100 ... 100,24
        ! s2_eta of size j: 24 > 22 ... 102, 22
        s2_zeta = v2*(abs(p_mesh(1:icmax+2,:)-2.0_dp*p_mesh(0:icmax+1,:)+p_mesh(-1:icmax,:)) & 
        & /(p_mesh(1:icmax+2,:)+p_mesh(0:icmax+1,:)+p_mesh(-1:icmax,:)))
        s2_eta = v2*(abs(p_mesh(:,-1:jcmax) - 2.0_dp * p_mesh(:,0:jcmax+1) +p_mesh(:,1:jcmax+2)) & 
        / (p_mesh(:,-1:jcmax) + 2.0_dp * p_mesh(:,0:jcmax+1) +p_mesh(:,1:jcmax+2)))

    end subroutine mesh_switches


    subroutine mesh_dissipation(q_mesh,p_mesh,nh,nv,lh,lv,v2,v4,a_mesh,diss_mesh,dt_max)

        real(dp), allocatable, dimension(:,:,:), intent(in) :: q_mesh,nh,nv
        real(dp), allocatable, dimension(:,:,:), intent(out) :: diss_mesh
        real(dp), allocatable, dimension(:,:), intent(out) :: dt_max
        real(dp), allocatable, dimension(:,:), intent(in) :: p_mesh, lh, lv, a_mesh
        real(dp), allocatable, dimension(:,:) :: u_mesh, v_mesh, s2_zeta, s2_eta, c_mesh, s2d_zeta, s2d_eta, & 
        & s4d_zeta, s4d_eta
        real(dp), allocatable, dimension(:,:,:) :: lm_mesh, dzq_mesh, deq_mesh, d3zq_mesh, d3eq_mesh
        real(dp), intent(in) :: v2, v4

        integer :: icmax,jcmax,i,j
        icmax = size(q_mesh, dim=1) - 4
        jcmax = size(q_mesh, dim=2) - 4

        allocate(u_mesh(-1:icmax+2,-1:jcmax+2))
        allocate(v_mesh(-1:icmax+2,-1:jcmax+2))

        u_mesh = q_mesh(:,:,2)/q_mesh(:,:,1)
        v_mesh = q_mesh(:,:,3)/q_mesh(:,:,1)
        
        allocate(c_mesh(-1:icmax+2,-1:jcmax+2))
        allocate(lm_mesh(icmax,jcmax,4))
        
        c_mesh(:,:) = sqrt(gamma*p_mesh(:,:)/q_mesh(:,:,1))

        ! SOUTH, NORTH, EAST, WEST
        do concurrent(j=1:jcmax, i=1:icmax)
            
            lm_mesh(i,j,4) = abs(u_mesh(i,j)*nh(i,j-1,1)+v_mesh(i,j)*nh(i,j-1,2)) + c_mesh(i,j)
            lm_mesh(i,j,3) = abs(u_mesh(i,j)*nh(i,j,1)+v_mesh(i,j)*nh(i,j,2)) + c_mesh(i,j) ! NORTH
            lm_mesh(i,j,2) = abs(u_mesh(i,j)*nv(i-1,j,1)+v_mesh(i,j)*nv(i-1,j,2)) + c_mesh(i,j) ! WEST
            lm_mesh(i,j,1) = abs(u_mesh(i,j)*nv(i,j,1)+v_mesh(i,j)*nv(i,j,2)) + c_mesh(i,j) ! EAST
        
        end do
        
        call mesh_switches(q_mesh,s2_zeta,s2_eta,v2)

        allocate(s2d_zeta(0:icmax,-1:jcmax+2))
        allocate(s2d_eta(-1:icmax+2,0:jcmax))
        
        s2d_zeta = 0.5_dp*(s2_zeta(0:icmax,:)+s2_zeta(1:icmax+1,:))
        s2d_eta = 0.5_dp*(s2_eta(:,0:jcmax)+s2_eta(:,1:jcmax+1))

        do concurrent(i=0:icmax,j=-1:jcmax+2)
            
            s2d_zeta(i,j) = max(s2d_zeta(i,j),0.02)
        
        end do

        do concurrent(i=-1:icmax+2,j=0:jcmax)
            s2d_eta(i,j) = max(s2d_zeta(i,j),0.02)
        end do

        allocate(dzq_mesh(-1:icmax+1,-1:jcmax+2,4))
        allocate(deq_mesh(-1:icmax+2,-1:jcmax+1,4))

        dzq_mesh = q_mesh(0:icmax+2,:,:) - q_mesh(-1:icmax+1,:,:)
        deq_mesh = q_mesh(:,0:jcmax+2,:) - q_mesh(:,-1:jcmax+1,:)

        allocate(d3zq_mesh(0:icmax,-1:jcmax+2,4))
        allocate(d3eq_mesh(-1:icmax+2,0:jcmax,4))
        
        d3zq_mesh = q_mesh(2:icmax+2,:,:) - 3.0_dp * q_mesh(1:icmax+1,:,:) + 3.0_dp * q_mesh(0:icmax,:,:) - q_mesh(-1:icmax-1,:,:)
        d3eq_mesh = q_mesh(:,2:jcmax+2,:) - 3.0_dp * q_mesh(:,1:jcmax+1,:) + 3.0_dp * q_mesh(:,0:jcmax,:) - q_mesh(:,-1:jcmax-1,:)

        allocate(s4d_zeta(0:icmax,-1:jcmax+2))
        allocate(s4d_eta(-1:icmax+2,0:jcmax))

        do concurrent(i=0:icmax, j=-1:jcmax+2)
            
            s4d_zeta(i,j) = max(0.15_dp*v4,v4-s2d_zeta(i,j))
        
        end do

        do concurrent(i=-1:icmax+2, j=0:jcmax)
            
            s4d_eta(i,j) = max(0.15_dp*v4,v4-s2d_eta(i,j))
        
        end do

        allocate(diss_mesh(1:icmax,1:jcmax,4))
        
        do concurrent(i=1:icmax, j=1:jcmax)
            
            diss_mesh(i,j,:) = & 
            &   s2d_zeta(i,j)*lv(i,j)*lm_mesh(i,j,1)*dzq_mesh(i,j,:) &
            & - s2d_zeta(i-1,j)*lv(i-1,j)*lm_mesh(i,j,2)*dzq_mesh(i-1,j,:) &
            & + s2d_eta(i,j)*lh(i,j)*lm_mesh(i,j,3)*deq_mesh(i,j,:) &
            & - s2d_eta(i,j-1)*lh(i,j-1)*lm_mesh(i,j,4)*deq_mesh(i,j-1,:) & 
            & - s4d_zeta(i,j)*lv(i,j)*lm_mesh(i,j,1)*d3zq_mesh(i,j,:) &
            & + s4d_zeta(i-1,j)*lv(i-1,j)*lm_mesh(i,j,2)*d3zq_mesh(i-1,j,:) &
            & - s4d_eta(i,j)*lh(i,j)*lm_mesh(i,j,3)*d3eq_mesh(i,j,:) &
            & + s4d_eta(i,j-1)*lh(i,j-1)*lm_mesh(i,j,4)*d3eq_mesh(i,j-1,:)
        
        end do

        allocate(dt_max(icmax,jcmax))

        do concurrent(i=1:icmax, j=1:jcmax)
            
            dt_max(i,j) = a_mesh(i,j)/(lm_mesh(i,j,1)*lv(i,j) & ! EAST
            & + lm_mesh(i,j,2)*lv(i-1,j) & ! WEST
            & + lm_mesh(i,j,3)*lh(i,j) & ! NORTH 
            & + lm_mesh(i,j,4)*lh(i-1,j) )! SOUTH
        
        end do

    end subroutine mesh_dissipation


    subroutine read_mach_area(filename, n, mach, area_ratio1)
        
        character(len=*), intent(in) :: filename
        integer, intent(in) :: n
        real(dp), intent(out) :: mach(n), area_ratio1(n)
    
        integer :: i, ios
        open(unit=10, file=filename, status='old', action='read', iostat=ios)
    
        do i = 1, n
            
            read(10, *, iostat=ios) mach(i), area_ratio1(i)
        
        end do
    
        close(10)
    
    end subroutine


    subroutine area_ratio(mach,ratio)

        real(dp), intent(in) :: mach
        real(dp), intent(out) :: ratio
        ratio = (1.0_dp/mach) * ((1.0_dp + (gamma - 1.0_dp)/2.0_dp * mach**2)**((gamma + 1.0_dp)/(2.0_dp*(gamma - 1.0_dp)))) &
        & / ((gamma + 1.0_dp)/2.0_dp)**((gamma + 1.0_dp)/(2.0_dp*(gamma - 1.0_dp)))

    end subroutine area_ratio


    ! finds mach number at a point for initialization
    subroutine point_mach(i,icmax,mach_array,ratio_array,mach_inlet,mach_point)

        integer, intent(in) :: i, icmax
        real(dp), intent(in) :: mach_inlet
        real(dp), intent(out) :: mach_point
        real(dp) :: point_a !ratio_con, ratio_div
        real(dp) :: x_pos, ri, ricmax, throat_l, throat_a, point_l
        real(dp) :: area_star, ratio_inlet, point_ratio, sub_root, super_root, area_star_inlet
        real(dp), allocatable, dimension(:), intent(in) :: mach_array, ratio_array
        integer :: arr_size, j

        ri = real(i, kind=dp)
        ricmax = real(icmax, kind=dp)
        x_pos = 5.0_dp*ri/ricmax
        throat_l = 0.76
        throat_a = throat_l

        if (x_pos>2 .and. x_pos<3) then
            
            point_l = 1.0_dp - 0.24_dp*sin(pi*(x_pos-2.0_dp))
            point_a = point_l
        
        else 
            
            point_a = 1.0
        
        end if

        call area_ratio(mach_inlet,ratio_inlet)

        area_star_inlet = 1.0_dp / ratio_inlet

        if (throat_a <= area_star_inlet) then
            
            area_star = throat_a
        
        else
            
            area_star = area_star_inlet
        
        end if

        point_ratio = point_a / area_star
        arr_size = size(mach_array)

        do j=1, arr_size-1
            
            if (point_ratio > ratio_array(j) .and. point_ratio < ratio_array(j+1)) then
                
                super_root = mach_array(j) + (point_ratio - ratio_array(j)) &
                & * (mach_array(j+1) - mach_array(j)) / (ratio_array(j+1) - ratio_array(j))
            
            else if(point_ratio < ratio_array(j) .and. point_ratio > ratio_array(j+1)) then
                
                sub_root = mach_array(j) + (point_ratio - ratio_array(j)) &
                & * (mach_array(j+1) - mach_array(j)) / (ratio_array(j+1) - ratio_array(j))
            
            end if
        end do

        if (ri/ricmax <= 0.5_dp) then
            
            mach_point = sub_root
            
        else
           
            if (throat_a <= area_star_inlet) then
                
                mach_point = super_root
            
            else 
                
                mach_point = sub_root
            
            end if

        end if

    end subroutine point_mach
    

    subroutine init_mesh(q_mesh,M,inlet_p,T,x,y)

        real(dp), allocatable, dimension(:,:,:), intent(inout) :: q_mesh
        real(dp), allocatable, dimension(:,:), intent(in) :: x, y
        real(dp), allocatable, dimension(:,:,:) :: nh, nv
        real(dp), allocatable, dimension(:) :: mach_i, mach_arr, ratio_arr
        real(dp), intent(in) :: M,inlet_p,T
        real(dp) :: c, u, inlet_rho, e1, e2,smoothdp,div, rho0, p0, rho_ratio, p_ratio
        real(dp) :: point_rho, point_p, pt_mach, point_vt, point_e, point_c
        integer :: icmax, jcmax, i, j, smooth

        icmax = size(q_mesh, dim=1)-4
        jcmax = size(q_mesh, dim=2)-4

        call mesh_normals(x,y,nh,nv)

        rho_ratio = (1.0_dp + (gamma - 1.0_dp)/2.0_dp * M**2)**(-1.0_dp/(gamma - 1.0_dp))

        p_ratio = (1.0_dp + (gamma - 1.0_dp)/2.0_dp * M**2)**(-gamma/(gamma - 1.0_dp))

        inlet_rho = inlet_p/R/T
        p0 = inlet_p/p_ratio
        rho0 = inlet_rho/rho_ratio

        c = sqrt(gamma*R*T)
        u = M*c
        inlet_rho = inlet_p/R/T
        e1 = inlet_p/(gamma-1)/inlet_rho
        e2 = e1 + 0.5_dp*u**2

        allocate(mach_arr(250))
        allocate(ratio_arr(250))

        call read_mach_area("ratio.dat", 250,mach_arr,ratio_arr)

        do i=-1,icmax
            
            if (i<1) then
                
                pt_mach = M
            
            else 
                
                call point_mach(i,icmax,mach_arr,ratio_arr,M,pt_mach) 
            
            end if
            
            rho_ratio = (1.0_dp + (gamma - 1.0_dp)/2.0_dp * pt_mach**2)**(-1.0_dp/(gamma - 1.0_dp))
            p_ratio = (1.0_dp + (gamma - 1.0_dp)/2.0_dp * pt_mach**2)**(-gamma/(gamma - 1.0_dp))
            point_rho = rho0*rho_ratio
            point_p = p0*p_ratio
            point_c = sqrt(gamma*point_p/point_rho)
            point_vt = pt_mach*point_c
            point_e = point_p/(point_rho*(gamma-1))+0.5_dp*point_vt**2
            
            do j=-1,jcmax+2
                
                q_mesh(i,j,:) = (/point_rho, point_rho*point_vt, 0.0_dp, point_rho*point_e/)
            
            end do
        end do

        q_mesh(icmax+1,:,:) = q_mesh(icmax,:,:)
        q_mesh(icmax+2,:,:) = q_mesh(icmax,:,:)

        q_mesh(-1,:,:) = q_mesh(0,:,:)
        
        smooth = 1
        smoothdp = 1.0
        div = smooth*4 + 2.0_dp

        do j=1, jcmax/2
            
            do i=smooth+1, icmax-smooth-1
                
                q_mesh(i,j,2) = q_mesh(i,j,2)*(sum(nh(i-smooth:i+smooth,j,2)) + sum(nh(i-smooth:i+smooth,j+1,2)))/div
                q_mesh(i,j,3) = -q_mesh(i,j,2)*(sum(nh(i-smooth:i+smooth,j,1)) + sum(nh(i-smooth:i+smooth,j+1,1)))/div
            
            end do
        
        end do

        do j=jcmax/2+1, jcmax
            
            do i=smooth+1, icmax-smooth-1
                
                q_mesh(i,j,2) = q_mesh(i,j,2)*(sum(nh(i-smooth:i+smooth,j,2)) + sum(nh(i-smooth:i+smooth,j-1,2)))/div
                q_mesh(i,j,3) = -q_mesh(i,j,2)*(sum(nh(i-smooth:i+smooth,j,1)) + sum(nh(i-smooth:i+smooth,j-1,1)))/div
            
            end do

        end do

        ! write(*,*) q_mesh(22,1,2:3)

    end subroutine init_mesh


    subroutine runge_kutta(q_mesh,x,y,CFL,M,p,T,n,v2,v4,l2_res)

        real(dp), allocatable, dimension(:,:,:), intent(inout) :: q_mesh
        real(dp), allocatable, dimension(:,:), intent(in) :: x,y
        real(dp) :: CFL
        real(dp), intent(in) :: M,p,T, v2,v4
        integer, intent(in) :: n

        real(dp), allocatable, dimension(:,:) :: a_mesh, dxh,dxv,dyh,dyv, pressure, lh,lv,dt_max
        real(dp), allocatable, dimension(:,:,:) :: nh,nv,residuals,dissipation,rhs

        real(dp), dimension(4) :: alpha, init_l2_res

        real(dp), dimension(n,4), intent(out) :: l2_res

        logical :: has_nan

        integer :: itr,sitr,icmax,jcmax,i,j,k,n_cells

        icmax = size(q_mesh, dim=1) - 4
        jcmax = size(q_mesh, dim=2) - 4
        n_cells = icmax*jcmax

        alpha = (/0.25_dp, 1.0_dp/3.0_dp, 0.5_dp, 1.0_dp/)

        call mesh_area(x,y,a_mesh)
        call mesh_delta(x,y,dxh,dxv,dyh,dyv)
        call mesh_normals(x,y,nh,nv)
        call mesh_face_lengths(x,y,lh,lv)
        
        allocate(rhs(icmax,jcmax,4))

        do itr=1, n
            if (mod(itr,100)==0) then
                write(*,*) "Itr:", itr
            end if

            call mesh_pressure(q_mesh,pressure)
            
            call mesh_dissipation(q_mesh,pressure,nh,nv,lh,lv,v2,v4,a_mesh,dissipation,dt_max)
            if (mod(itr,20000)==0) then
                CFL = CFL*1.1
            end if
            
            ! write(*,*) dissipation(2,:,1)
            do sitr=1, 4

                call inlet_boundary(q_mesh,M,p,T)
                
                call outlet_boundary(q_mesh,p)
                
                call wall_boundary(q_mesh,nh)
                
                call mesh_residuals(q_mesh,dxh,dxv,dyh,dyv,residuals)

                do concurrent(k=1:4, j=1:jcmax, i=1:icmax)
                    rhs(i,j,k) = alpha(sitr)*CFL*dt_max(i,j)*(residuals(i,j,k)-dissipation(i,j,k))/a_mesh(i,j)
                end do
                
                q_mesh(1:icmax,1:jcmax,:) = q_mesh(1:icmax,1:jcmax,:) - rhs
            
            end do
            
            if (any(q_mesh /= q_mesh)) then
                
                print *, "NaN found", itr
                stop 1
            
            end if

            if (itr==1) then
                
                do k=1,4
                    
                    init_l2_res(k) = sqrt(sum(residuals(:,:,k)**2)/n_cells)
                    l2_res(1,k) = 1.0_dp 
                
                end do
            
            else 
                
                do k=1,4
                    
                    l2_res(itr,k) = sqrt(sum(residuals(:,:,k)**2)/n_cells)/init_l2_res(k)
                
                end do
            
            end if
            
        end do
        
        write(*,*) "Finished!"

    end subroutine runge_kutta


    subroutine bump_forces(q_mesh,b_forces,x,y)

        real(dp), allocatable, dimension(:,:,:), intent(in) :: q_mesh
        real(dp), allocatable, dimension(:,:), intent(in) :: x,y
        real(dp), allocatable, dimension(:,:,:) :: nh,nv
        real(dp), allocatable, dimension(:,:) :: face_lengths, dv
        real(dp), dimension(2,2), intent(out) :: b_forces
        real(dp), allocatable,dimension(:,:) :: p_mesh
        real(dp), allocatable, dimension(:,:) :: top_force, bot_force
        integer :: icmax,jcmax,i,bump_start,bump_length,k

        icmax = size(q_mesh,dim=1) -4
        jcmax = size(q_mesh,dim=2) -4
        bump_length = icmax/5
        bump_start = bump_length*2
        
        call mesh_normals(x,y,nh,nv)
        
        call mesh_face_lengths(x,y,face_lengths,dv)
        
        call mesh_pressure(q_mesh,p_mesh)
        allocate(top_force(bump_start:bump_start+bump_length,2))
        allocate(bot_force(bump_start:bump_start+bump_length,2))

        do i=bump_start, bump_start+bump_length
            top_force(i,:) = (1.875_dp*p_mesh(i,1)-1.25_dp*p_mesh(i,jcmax-1)+0.375_dp*p_mesh(i,jcmax-2))& 
            & *nh(i,jcmax,:)*face_lengths(i,jcmax)
            bot_force(i,:) = (1.875_dp*p_mesh(i,1)-1.25_dp*p_mesh(i,2)+0.375_dp*p_mesh(i,3))*nh(i,0,:) &
            & *face_lengths(i,jcmax)
        end do
        
        do k=1,2
            b_forces(1,k) = sum(bot_force(:,k))
            b_forces(2,k) = sum(top_force(:,k))
        end do

    end subroutine bump_forces
    

    subroutine read_tecplot_mesh(filename, x, y)
        !-----------------------------------------------------------------------
        ! Reads a simple Tecplot file (POINT format) into x and y arrays.
        ! Assumes correct file format and existence. Minimal error checks.
        !-----------------------------------------------------------------------
    
        character(len=*), intent(in) :: filename
        real(8), allocatable, dimension(:,:), intent(out) :: x, y
    
        ! Local variables
        integer :: unit_num
        integer :: ni, nj, i, j
        character(len=256) :: line        ! Buffer to read lines
        integer :: idx_i, idx_j
    
        ! Use a free unit number (Fortran 2008+) or choose one (e.g., 10)
        open(newunit=unit_num, file=filename, status='old', action='read')
    
        ! --- Read Header (assuming correct format) ---
        read(unit_num, '(A)') line ! Read and discard TITLE
        read(unit_num, '(A)') line ! Read and discard VARIABLES
        read(unit_num, '(A)') line ! Read ZONE line

        ! --- Extract I and J from ZONE line (simple parsing) ---
        ! Find 'I=' and 'J='
        idx_i = index(line, 'I=')
        idx_j = index(line, 'J=')
    
        ! Read the integer values starting right after '='
        read(line(idx_i+2:), *) ni
        read(line(idx_j+2:), *) nj
    
        ! --- Allocate Arrays ---
        allocate(x(ni, nj), y(ni, nj))
    
        ! --- Read Data ---
        ! Tecplot POINT format lists data with I varying fastest
        do j = 1, nj
            do i = 1, ni
                read(unit_num, *) x(i, j), y(i, j)
            end do
        end do
    
        ! --- Clean up ---
        close(unit_num)
    
    end subroutine read_tecplot_mesh


    subroutine write_res( filename, l2_res )
        character(len=*), intent(in)  :: filename
        real(dp), allocatable, dimension(:,:), intent(in)  :: l2_res
     
        integer :: n, i, u
        n = size(l2_res,1)
        u = 10
     
        open(unit=u, file=filename, status='replace', action='write')
        do i = 1, n
           write(u,'(4ES22.14)') l2_res(i,1), l2_res(i,2), l2_res(i,3), l2_res(i,4)
        end do
        close(u)
     end subroutine write_res
     


    subroutine write_py_sparse_q(x, y, q, imax, jmax, basename)
        implicit none
    
        integer, intent(in) :: imax, jmax
        real(kind=8), dimension(imax, jmax), intent(in) :: x, y
        ! q remains cell-centered
        real(kind=8), dimension(imax-1, jmax-1, 4), intent(in) :: q
        character(len=*), intent(in) :: basename
    
        integer :: i, j, k, file_unit
        integer :: icells, jcells
        character(len=20) :: var_name
        character(len=255) :: filename
        character(len=20), parameter :: nan_placeholder = '       NaN         ' ! Padded
    
        ! Short names suitable for headers/filenames
        character(len=8), dimension(4) :: var_names = &
            [character(len=8) :: "Density", "U", "V", "Energy"]
    
        icells = imax - 1
        jcells = jmax - 1
    
        do k = 1, 4
            var_name = trim(var_names(k))
            filename = trim(basename) // "_" // trim(adjustl(var_name)) // ".txt" ! Suffix
    
            open(newunit=file_unit, file=filename, status='replace', action='write')
    
            ! Write header line including NODAL dimensions
            write(file_unit, '(A, A, A, A, I0, A, I0)') '# X_node Y_node ', &
                trim(var_name), '(at node i,j if cell exists)', &
                ' INodes=', imax, ' JNodes=', jmax
    
            ! Loop through ALL nodes (I inner -> column-major write order)
            do j = 1, jmax
                do i = 1, imax
                    ! Check if a cell exists for which this node is the bottom-left corner
                    if (i <= icells .and. j <= jcells) then
                        ! Write X, Y, and the corresponding Q value
                        write(file_unit, '(3E20.10)') x(i,j), y(i,j), q(i,j,k)
                    else
                        ! Write X, Y, and the placeholder string for Q
                        write(file_unit, '(2E20.10, A)') x(i,j), y(i,j), nan_placeholder
                    end if
                end do
                ! Optional: Add a blank line between 'rows' of nodes (j-levels) for visual clarity?
                ! write(file_unit, *)
            end do
    
            close(file_unit)
        end do
    
    end subroutine write_py_sparse_q
      

end module Utils

