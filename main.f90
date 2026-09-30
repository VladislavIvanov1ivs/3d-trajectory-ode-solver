program main_move_programm
    use vector_m
    use int_m
    use forces_m
    implicit none
    
    type(state) :: y, y_prev
    real(8):: t, dt, tmax, t_prev
    integer:: nsteps, i 
    integer:: scenario, scenario_prev 
    integer:: method 
    integer:: stop_point
    integer :: ground_choice

    ! Меню параметров
    real(8):: gval, k1, k2, mu 

    ! Начальные условия 
    real(8):: V0, az_deg, el_deg, x0, y0, z0
    character(len=*), parameter :: outfile = "traektor.txt"
    real(8):: a_x, a_y, a_z
    logical :: hit_ground, use_ground_stop
    real(8), parameter :: y_ground = 0d0
    
    scenario_prev=-1
    ! Выбор сценария
    print*, "Hello you are in move_program!!!"
    print*, "Please, enter start position x0, y0, z0"
    read(*,*) x0, y0, z0 
    y%r = v3(x0, y0, z0)
    print*, "Now enter start value of speed V0"
    read(*,*) V0
    print*, "Do you want enter angles value in rad(1) or deg(2)?"
    read*, stop_point
    print*, "Now you shoud enter azimut angle (xy_plot)"
    print*, "Azimut angle positive in this direction: from X+ to Y+"
    read*, az_deg
    print*, "Now you shoud enter el angle (Z ax.)"
    print*, "El angle positive in this direction:Z+"
    read*, el_deg

    if (stop_point==1) then
        y%v = from_angle(V0, az_deg, el_deg)
    else
        y%v =from_angle_deg(V0, az_deg, el_deg)
    end if

    stop_point=1
    use_ground_stop = .false.
    call disable_all()
    do while (stop_point==1)
        print*, " 0 - No forces (a==0)"
        print*, " 1 - Near the Easrt (set gravity_const)"
        print*, " 2 - Linear resistance ( a = -k * v )"
        print*, " 3 - Quad resistance ( a = -k * |v| * v )"
        print*, " 4 - Orbit move ( a = -mu * r / |r|^3 )"
        print*, " 5 - There is a craving (for exmpl rocket)"
        print*, "Enter your choice: "
        read*, scenario
        if (scenario == scenario_prev) then
            print*, "You picked the same choise"
        else
            select case(scenario)
                case(0)
                    stop_point = 0
                case(1)
                    print*, "Enter g (positive) == "
                    read(*,*) gval 
                    call set_gravity(v3(0d0, -abs(gval), 0d0))
                    call enable_gravity(.true.)
                case(2)
                    print*, "Enter k1, for linear resistance"
                    read(*,*) k1
                    call set_drag_linear(k1)
                    call enable_drag_linear(.true.)
                case(3)
                    print*, "Enter k2, for linear resistance"
                    read(*,*) k2
                    call set_drag_quad(k2)
                    call enable_drag_quad(.true.)
                case(4)
                    print*, "Enter mu, (for exmpl 3.986e14 for Earth): mu =="
                    read(*,*) mu
                    call set_orbit_mu(mu)
                    call enable_orbit_mu(.true.)
                case(5)
                    print*, "Enter a_x, a_y, a_z for craving"
                    read(*,*) a_x, a_y, a_z
                    call set_thrust(v3(a_x , a_y, a_z))
                    call enable_thrust(.true.)
                case default
                    print*, "oops wrong number, try again!"
                end select
        end if  
        scenario_prev= scenario
        Print*, "!!! if you have selected everything you wanted, press 0, or press 1 !!!"
        read*, stop_point
    end do

    Print*, "First cond are saved"
    print*, "Do you want to enable ground collision at y = 0? (1 = yes, 0 = no)"
    read(*,*) ground_choice
    if (ground_choice == 1) then
        use_ground_stop = .true.
    else
        use_ground_stop = .false.
    end if

    print*, "What type of method do you wnat?"
    print*, "1 - Euler "
    print*, "2 - RK4 "
    read(*,*) method

    print*, "Enter dt and tmax (for exmpl 0.01 20)"
    read(*,*) dt, tmax 

    if( dt<= 0d0 .or. tmax <= 0d0) stop "dt and tmax must be positive"
    nsteps= int(tmax/dt)

    if (use_ground_stop .and. y%r%y < y_ground) then
        print*, "Start point is below the ground level y =", y_ground
        y%r%y = y_ground
        if (y%v%y < 0d0) y%v%y = 0d0
    end if

    open(unit=10 , file= outfile, status="replace", action="write")
    write(10, '(a)') "t,x,y,z,vx,vy,vz"
    t=0d0
    hit_ground = .false.

    do i= 0, nsteps
        write(10, '(f12.6,1x,6(es18.10,1x))') t, y%r%x,  y%r%y,  y%r%z,  y%v%x, y%v%y, y%v%z
        if (i == nsteps) exit

        y_prev = y
        t_prev = t

        select case(method)
            case(1)
                call step_euler(t, y, dt, accel_total)
            case(2)
                call step_rk4(t, y, dt, accel_total)
            case default
                stop "Uknow method!"
        end select

        if (use_ground_stop) then
            call apply_ground_stop(t_prev, y_prev, t, y, y_ground, hit_ground)
            if (hit_ground) then
                write(10, '(f12.6,1x,6(es18.10,1x))') t, y%r%x,  y%r%y,  y%r%z,  y%v%x, y%v%y, y%v%z
                exit
            end if
        end if
    end do
    close(10)
    print*, "Done!!!. Output file written to", outfile
    read(*,*) scenario
contains
    subroutine apply_ground_stop(t0, y0, t1, y1, ground_y, hit_ground)
        real(8), intent(in) :: t0, ground_y
        type(state), intent(in) :: y0
        type(state), intent(inout) :: y1
        real(8), intent(inout) :: t1
        logical, intent(out) :: hit_ground
        real(8) :: alpha, dy

        hit_ground = .false.

        if (y0%r%y > ground_y .and. y1%r%y <= ground_y) then
            dy = y1%r%y - y0%r%y
            if (abs(dy) > 1d-14) then
                alpha = (ground_y - y0%r%y) / dy
            else
                alpha = 1d0
            end if
            alpha = max(0d0, min(1d0, alpha))

            t1   = t0   + alpha*(t1   - t0)
            y1%r = y0%r + alpha*(y1%r - y0%r)
            y1%v = y0%v + alpha*(y1%v - y0%v)
            y1%r%y = ground_y
            if (y1%v%y < 0d0) y1%v%y = 0d0
            hit_ground = .true.
        end if
    end subroutine apply_ground_stop
end program main_move_programm
