module forces_m
    use vector_m
    implicit none
    private
    ! Настройка параметров
    public:: set_gravity, set_drag_linear, set_drag_quad, set_orbit_mu, set_thrust
    ! Объявление моделей ускорения
    public:: enable_gravity, enable_drag_linear, enable_drag_quad, enable_orbit_mu, enable_thrust
    public:: disable_all
    public:: a_gravity, a_drag_linear, a_drag_quad, a_orbit, a_thrust
    public:: accel_total
    !Объявление общей модели ускорения
    ! Физические параметры
    type(vec3):: thrust_ = vec3(0d0, 0d0, 0d0) ! тяговое ускорение (ракета)
    type(vec3):: g_ = vec3(0d0, -9.81d0, 0d0) ! Вектор ускорения свободного падения
    real(8) :: k1_ = 0d0 !Линейное сопротивление: a = -k1*V
    real(8) :: k2_ = 0d0 !квадратичное сопротивление: a= -k2*V*V
    real(8) :: mu_ = 0d0 ! грав. параметр для орбит: a= -mu*r/|r|^3
    ! Флаги включения компонент ускорения 
    logical :: use_gravity = .false.
    logical :: use_orbit = .false.
    logical :: use_drag_linear = .false.
    logical :: use_drag_quad = .false.
    logical :: use_thrust = .false.
contains
    ! Сеттеры
    subroutine set_gravity(g)
        type(vec3), intent(in):: g
        g_=g
    end subroutine set_gravity
    
    subroutine set_drag_linear(k1)
        real(8), intent(in)::k1
        k1_=k1
    end subroutine set_drag_linear
    
    subroutine set_drag_quad(k2)
        real(8), intent(in)::k2
        k2_=k2
    end subroutine set_drag_quad
    
    subroutine set_orbit_mu(mu)
        real(8), intent(in)::mu
        mu_=mu
    end subroutine set_orbit_mu
    
    subroutine set_thrust(thrust_a)
        type(vec3), intent(in):: thrust_a
        thrust_=thrust_a
    end subroutine set_thrust
    ! Управление флагами
    subroutine disable_all()
        use_gravity = .false.
        use_orbit = .false.
        use_drag_linear = .false.
        use_drag_quad = .false.
        use_thrust = .false.
    end subroutine disable_all
    
    subroutine enable_gravity(on)
        logical, intent(in):: on
        use_gravity = on 
    end subroutine enable_gravity

    subroutine enable_drag_linear(on)
        logical, intent(in):: on
        use_drag_linear = on 
    end subroutine enable_drag_linear

    subroutine enable_drag_quad(on)
        logical, intent(in):: on
        use_drag_quad = on 
    end subroutine enable_drag_quad

    subroutine enable_thrust(on)
        logical, intent(in):: on
        use_thrust = on 
    end subroutine enable_thrust

    subroutine enable_orbit_mu(on)
        logical, intent(in):: on
        use_orbit = on 
    end subroutine enable_orbit_mu

    ! Компоненты ускорения 
    function a_gravity() result(a)
        type(vec3):: a
        a = g_
    end function a_gravity

    function a_drag_linear(v) result(a)
        type(vec3), intent(in):: v
        type(vec3):: a
        a = (-k1_)*v
    end function a_drag_linear

    function a_drag_quad(v) result(a)
        type(vec3), intent(in):: v
        type(vec3):: a
        a = (-k2_*norm(v))*v
    end function a_drag_quad

    function a_orbit(r) result(a)
        type(vec3), intent(in):: r
        type(vec3):: a
        real(8) :: r2, rabs, inv_r3
        r2 = norm2(r)
        if(r2 == 0d0) then 
            a= zero()
        else 
            rabs = sqrt(r2)
            inv_r3= 1d0/(r2*rabs)
            a = (-mu_*inv_r3)*r 
        end if 
    end function a_orbit

    function a_thrust() result(a)
        type(vec3):: a 
        a = thrust_
    end function a_thrust

    ! Итоговая сборка
    function accel_total(t,r,v) result(a)
        real(8), intent(in):: t
        type(vec3), intent(in):: r, v 
        type(vec3):: a 
        a = zero()

        if (use_gravity)        a= a+a_gravity()
        if (use_drag_linear)    a= a+a_drag_linear(v)
        if (use_drag_quad)      a= a+a_drag_quad(v)
        if (use_orbit)          a= a+a_orbit(r)
        if (use_thrust)         a= a+a_thrust() 
    end function accel_total
end module forces_m