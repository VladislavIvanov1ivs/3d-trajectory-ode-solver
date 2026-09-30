module diff_m
    use vector_m
    implicit none
    private
    ! Объявляем методы дифференцирования (Конечные разности)
    ! d1- первая производная
    ! d2 - вторая производная 
    public:: d1_forward, d1_backward, d1_central
    public:: d2_central
contains
    ! Конечная разность, шаг вперед
    pure function d1_forward(f_now, f_next, dt) result (df)
        type(vec3), intent(in):: f_now, f_next
        real(8), intent(in):: dt
        type(vec3):: df
        df = (f_next-f_now)*(1d0/dt)
    end function d1_forward
    ! Конечная разность шаг назад
    pure function d1_backward(f_now, f_prev, dt) result (df)
        type(vec3), intent(in):: f_now, f_prev
        real(8), intent(in):: dt
        type(vec3):: df
        df = (f_now-f_prev)*(1d0/dt)
    end function d1_backward
    ! Конечная разность центральная схема
    pure function d1_central(f_next, f_prev, dt) result (df)
        type(vec3), intent(in):: f_next, f_prev
        real(8), intent(in):: dt
        type(vec3):: df
        df = (f_next-f_prev)*(0.5d0/dt)
    end function d1_central
    ! Конечная разность для второй производной центральная схема
    pure function d2_central(f_next, f_now, f_prev, dt) result (df)
        type(vec3), intent(in):: f_next, f_prev, f_now
        real(8), intent(in):: dt
        type(vec3):: df
        df = (f_next-2d0*f_now+f_prev)*(1d0/(dt*dt))
    end function d2_central
end module diff_m