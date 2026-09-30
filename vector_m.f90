module vector_m
    implicit none
    private
    ! Vec3 - 3d вектор
    ! v3 - конструктор вектора 
    ! zero - функция обнуления вектора
    ! dot - скалярное произведение
    ! cross - векторное произведение 
    ! norm - длина ветора 
    ! norm2 - квадрат длины вектора 
    public:: vec3, v3, zero, dot, cross, norm, norm2
    ! Вспомогательные функции перехода от угла и модуля к составляющим вектора
    public:: from_angle, from_angle_deg
    ! Фиксации действий над ветором (оператор умножения - это умножение на скаляр)
    public:: operator(+),operator(-),operator(*)

    type:: vec3
        real(8):: x=0d0 , y= 0d0, z= 0d0
    end type vec3
    ! Объявляем объединяющий интерфейс
    ! Те за каждым оператором скрывается функция 
    interface operator(+); module procedure add3; end interface
    interface operator(-); module procedure sub3; end interface
    interface operator(*); module procedure smul3_r, smul3_l; end interface
contains
    ! Конструктор вектора
    pure function v3(x,y,z) result(v)
        real(8), intent(in):: x, y, z
        type(vec3):: v
        v%x=x
        v%y=y
        v%z=z
    end function v3
    ! Обнуление вектора
    pure function zero() result(v)
        type(vec3):: v
        v = v3(0d0,0d0,0d0)
    end function zero
    ! Функция сложения
    pure function add3(a,b) result(c)
        type(vec3), intent(in):: a, b
        type(vec3):: c
        c%x= a%x+b%x
        c%y= a%y+b%y
        c%z= a%z+b%z
    end function add3
    ! Функция вычетания
    pure function sub3(a,b) result(c)
        type(vec3), intent(in):: a, b
        type(vec3):: c
        c%x= a%x-b%x
        c%y= a%y-b%y
        c%z= a%z-b%z
    end function sub3
    !Функция умножения на скаляр слева (если скаляр слева)
    pure function smul3_l(k,a) result(c)
        real(8), intent(in):: k
        type(vec3), intent(in):: a 
        type(vec3):: c 
        c%x=k*a%x
        c%y=k*a%y
        c%z=k*a%z
    end function smul3_l
    !Функция умножения на скаляр справа (если скаляр справа)
    pure function smul3_r(a,k) result(c)
        real(8), intent(in):: k
        type(vec3), intent(in):: a 
        type(vec3):: c 
        c= k*a
    end function smul3_r
    ! Функция скалярного произведения
    pure function dot(a,b) result(d)
        type(vec3), intent(in):: a, b
        real(8):: d
        d= a%x*b%x+a%y*b%y+a%z*b%z
    end function dot
    ! Функция векторного произведения
    pure function cross(a,b) result(c)
        type(vec3), intent(in):: a, b
        type(vec3):: c 
        c%x= a%y*b%z - b%y*a%z
        c%y = a%z*b%x - b%z*a%x
        c%z = a%x*b%y - b%x*a%y
    end function cross
    ! Функция квадрата длины вектора
    pure function norm2(a) result(n2)
        type(vec3), intent(in):: a 
        real(8):: n2
        n2= dot(a,a)
    end function norm2 
    ! Функция длины вектора
    pure function norm(a) result(n)
        type(vec3), intent(in):: a 
        real(8):: n
        n= sqrt(norm2(a))
    end function norm
    ! Функция перехода от модуля к составляющим по углу в радианах
    pure function from_angle(a, alpha_xy, alpha_z) result(v)
        real(8), intent(in):: a, alpha_xy, alpha_z
        type(vec3):: v 
        v= v3(a*cos(alpha_z)*cos(alpha_xy), a*sin(alpha_xy)*cos(alpha_z), a*sin(alpha_z))
    end function from_angle
    ! Функция перехода от модуля к составляющим по углу в градусах
    pure function from_angle_deg(a, alpha_xy_deg, alpha_z_deg) result(v)
        real(8), intent(in):: a, alpha_xy_deg, alpha_z_deg
        type(vec3):: v
        real(8):: alp1, alp2 
        real(8), parameter:: pi= acos(-1.0d0)
        alp1= alpha_xy_deg*pi/180d0
        alp2= alpha_z_deg*pi/180d0
        v= from_angle(a,alp1, alp2)
    end function from_angle_deg
end module vector_m