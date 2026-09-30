module int_m
    use vector_m
    implicit none
    private
    ! Объявление состояния системы, ускорения как функции, функция шага интегрирования по Эйлеру
    ! И по Рунге-Кутта 4 порядка 
    public:: state, accel_func, step_euler, step_rk4
    type :: state
        type(vec3):: r 
        type(vec3):: v 
    end type
    ! Далее интерфейс функции ускорения (каркас для возможности варьировать)
    abstract interface
        ! на вход подаем время, радиус вектор, скорость ожидая на выходе вектор
        function accel_func(t,r,v) result (a)
            import:: vec3 
            real(8), intent(in):: t
            type(vec3), intent(in)::  r, v
            type(vec3):: a
        end function accel_func
    end interface
contains
        ! Функция правой части системы 
        ! dr/dt = v
        ! dv/dt = a(t,r,v)
        function rhs(t,y,accel) result(dy)
            real(8), intent(in):: t 
            type(state), intent(in):: y
            procedure(accel_func):: accel
            type(state) :: dy 
            dy%r = y%v 
            dy%v = accel(t, y%r , y%v )
        end function rhs
        ! Шаг по методу Эйлера
        subroutine step_euler(t, y, dt, accel)
            real(8), intent(inout):: t 
            real(8), intent(in):: dt
            type(state), intent(inout):: y
            procedure(accel_func) :: accel
            type(state):: k1
            k1 = rhs(t,y,accel)
            y%r= y%r + dt*k1%r 
            y%v= y%v + dt*k1%v

            t= t+dt  
        end subroutine step_euler
        ! Шаг по методу Рунге-Кутта 4 порядка
        subroutine step_rk4(t, y, dt, accel)
            real(8), intent(inout):: t 
            real(8), intent(in):: dt
            type(state), intent(inout):: y
            procedure(accel_func) :: accel
            type(state):: k1, k2, k3, k4, yt

            k1 = rhs(t, y, accel)

            yt%r = y%r + 0.5d0*dt*k1%r
            yt%v = y%v + 0.5d0*dt*k1%v
            k2 = rhs(t + 0.5d0*dt, yt, accel)

            yt%r = y%r + 0.5d0*dt*k2%r
            yt%v = y%v + 0.5d0*dt*k2%v
            k3 = rhs(t + 0.5d0*dt, yt, accel)

            yt%r = y%r + dt*k3%r
            yt%v = y%v + dt*k3%v
            k4 = rhs(t + dt, yt, accel)

            y%r = y%r + (dt/6d0)*(k1%r + 2d0*k2%r + 2d0*k3%r + k4%r)
            y%v = y%v + (dt/6d0)*(k1%v + 2d0*k2%v + 2d0*k3%v + k4%v)

            t = t + dt
        end subroutine step_rk4
end module int_m
