!  Calculations of forces from stresses
subroutine fl_node
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

    !  1 - 3
    !  |   |
    !  2 - 4
    !
    !  diagonal / :
    !
    !   A:        B:
    !
    !  1---3         1
    !  | /         / |
    !  2         2---3
    !
    !  diagonal \ :
    !
    !   C:        D:
    !
    !  1          1---3
    !  | \         \  |
    !  2---3          2
    !
    !    assemblage of forces is COUNTRE CLOCK-WISE !
    !

boff = 0

drat = dt / dt_elastic
if (drat .lt. 1.) drat = 1.

!$OMP parallel private(i, j, fx, fy, &
!$OMP                  p_est, rosubg, &
!$OMP                  press_norm_l, dlx_l, dly_l, &
!$OMP                  press_norm_r, dlx_r, dly_r, &
!$OMP                  iunknown, rho_water_g, water_depth)
!
!$OMP do
do i = 1,nx
    do j = 1,nz
        if(ynstressbc.eq.0.) then
           force(j,i,1) = 0
           force(j,i,2) = 0
           balance(j,i,1)=0
           balance(j,i,2)=0
        endif
        ! REGULAR PART - forces from stresses
        
        ! Element (j-1,i-1). Triangles B,C,D
        if ( j.ne.1 .and. i.ne.1 ) then
            ! triangle B
            ! side 2-3
            fx = stress0(j-1,i-1,1,2) * (cord(j  ,i  ,2)-cord(j  ,i-1,2)) - &
                 stress0(j-1,i-1,3,2) * (cord(j  ,i  ,1)-cord(j  ,i-1,1))
            fy = stress0(j-1,i-1,3,2) * (cord(j  ,i  ,2)-cord(j  ,i-1,2)) - &
                 stress0(j-1,i-1,2,2) * (cord(j  ,i  ,1)-cord(j  ,i-1,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 3-1
            fx = stress0(j-1,i-1,1,2) * (cord(j-1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i-1,3,2) * (cord(j-1,i  ,1)-cord(j  ,i  ,1))
            fy = stress0(j-1,i-1,3,2) * (cord(j-1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i-1,2,2) * (cord(j-1,i  ,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

            ! triangle C
            ! side 2-3
            fx = stress0(j-1,i-1,1,3) * (cord(j  ,i  ,2)-cord(j  ,i-1,2)) - &
                 stress0(j-1,i-1,3,3) * (cord(j  ,i  ,1)-cord(j  ,i-1,1))
            fy = stress0(j-1,i-1,3,3) * (cord(j  ,i  ,2)-cord(j  ,i-1,2)) - &
                 stress0(j-1,i-1,2,3) * (cord(j  ,i  ,1)-cord(j  ,i-1,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 3-1
            fx = stress0(j-1,i-1,1,3) * (cord(j-1,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i-1,3,3) * (cord(j-1,i-1,1)-cord(j  ,i  ,1))
            fy = stress0(j-1,i-1,3,3) * (cord(j-1,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i-1,2,3) * (cord(j-1,i-1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

            ! triangle D
            ! side 1-2
            fx = stress0(j-1,i-1,1,4) * (cord(j  ,i  ,2)-cord(j-1,i-1,2)) - &
                 stress0(j-1,i-1,3,4) * (cord(j  ,i  ,1)-cord(j-1,i-1,1))
            fy = stress0(j-1,i-1,3,4) * (cord(j  ,i  ,2)-cord(j-1,i-1,2)) - &
                 stress0(j-1,i-1,2,4) * (cord(j  ,i  ,1)-cord(j-1,i-1,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 2-3
            fx = stress0(j-1,i-1,1,4) * (cord(j-1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i-1,3,4) * (cord(j-1,i  ,1)-cord(j  ,i  ,1))
            fy = stress0(j-1,i-1,3,4) * (cord(j-1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i-1,2,4) * (cord(j-1,i  ,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

        endif

        ! Element (j-1,i). Triangles A,B,C.
        if ( j.ne.1 .and. i.ne.nx ) then
            ! triangle A
            ! side 1-2
            fx = stress0(j-1,i  ,1,1) * (cord(j  ,i  ,2)-cord(j-1,i  ,2)) - &
                 stress0(j-1,i  ,3,1) * (cord(j  ,i  ,1)-cord(j-1,i  ,1))
            fy = stress0(j-1,i  ,3,1) * (cord(j  ,i  ,2)-cord(j-1,i  ,2)) - &
                 stress0(j-1,i  ,2,1) * (cord(j  ,i  ,1)-cord(j-1,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 2-3
            fx = stress0(j-1,i  ,1,1) * (cord(j-1,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i  ,3,1) * (cord(j-1,i+1,1)-cord(j  ,i  ,1))
            fy = stress0(j-1,i  ,3,1) * (cord(j-1,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i  ,2,1) * (cord(j-1,i+1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

            ! triangle B
            ! side 1-2
            fx = stress0(j-1,i  ,1,2) * (cord(j  ,i  ,2)-cord(j-1,i+1,2)) - &
                 stress0(j-1,i  ,3,2) * (cord(j  ,i  ,1)-cord(j-1,i+1,1))
            fy = stress0(j-1,i  ,3,2) * (cord(j  ,i  ,2)-cord(j-1,i+1,2)) - &
                 stress0(j-1,i  ,2,2) * (cord(j  ,i  ,1)-cord(j-1,i+1,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 2-3
            fx = stress0(j-1,i  ,1,2) * (cord(j  ,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i  ,3,2) * (cord(j  ,i+1,1)-cord(j  ,i  ,1))
            fy = stress0(j-1,i  ,3,2) * (cord(j  ,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i  ,2,2) * (cord(j  ,i+1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

            ! triangle C
            ! side 1-2
            fx = stress0(j-1,i  ,1,3) * (cord(j  ,i  ,2)-cord(j-1,i  ,2)) - &
                 stress0(j-1,i  ,3,3) * (cord(j  ,i  ,1)-cord(j-1,i  ,1))
            fy = stress0(j-1,i  ,3,3) * (cord(j  ,i  ,2)-cord(j-1,i  ,2)) - &
                 stress0(j-1,i  ,2,3) * (cord(j  ,i  ,1)-cord(j-1,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 2-3
            fx = stress0(j-1,i  ,1,3) * (cord(j  ,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i  ,3,3) * (cord(j  ,i+1,1)-cord(j  ,i  ,1))
            fy = stress0(j-1,i  ,3,3) * (cord(j  ,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(j-1,i  ,2,3) * (cord(j  ,i+1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

        endif
        
        ! Element (j,i-1). Triangles A,B,D
        if ( j.ne.nz .and. i.ne.1 ) then
            ! triangle A
            ! side 2-3
            fx = stress0(j  ,i-1,1,1) * (cord(j  ,i  ,2)-cord(j+1,i-1,2)) - &
                 stress0(j  ,i-1,3,1) * (cord(j  ,i  ,1)-cord(j+1,i-1,1))
            fy = stress0(j  ,i-1,3,1) * (cord(j  ,i  ,2)-cord(j+1,i-1,2)) - &
                 stress0(j  ,i-1,2,1) * (cord(j  ,i  ,1)-cord(j+1,i-1,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 3-1
            fx = stress0(j  ,i-1,1,1) * (cord(j  ,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i-1,3,1) * (cord(j  ,i-1,1)-cord(j  ,i  ,1))
            fy = stress0(j  ,i-1,3,1) * (cord(j  ,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i-1,2,1) * (cord(j  ,i-1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

            ! triangle B
            ! side 1-2
            fx = stress0(j  ,i-1,1,2) * (cord(j+1,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i-1,3,2) * (cord(j+1,i-1,1)-cord(j  ,i  ,1))
            fy = stress0(j  ,i-1,3,2) * (cord(j+1,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i-1,2,2) * (cord(j+1,i-1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 3-1
            fx = stress0(j  ,i-1,1,2) * (cord(j  ,i  ,2)-cord(j+1,i  ,2)) - &
                 stress0(j  ,i-1,3,2) * (cord(j  ,i  ,1)-cord(j+1,i  ,1))
            fy = stress0(j  ,i-1,3,2) * (cord(j  ,i  ,2)-cord(j+1,i  ,2)) - &
                 stress0(j  ,i-1,2,2) * (cord(j  ,i  ,1)-cord(j+1,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

            ! triangle D
            ! side 2-3
            fx = stress0(j  ,i-1,1,4) * (cord(j  ,i  ,2)-cord(j+1,i  ,2)) - &
                 stress0(j  ,i-1,3,4) * (cord(j  ,i  ,1)-cord(j+1,i  ,1))
            fy = stress0(j  ,i-1,3,4) * (cord(j  ,i  ,2)-cord(j+1,i  ,2)) - &
                 stress0(j  ,i-1,2,4) * (cord(j  ,i  ,1)-cord(j+1,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 3-1
            fx = stress0(j  ,i-1,1,4) * (cord(j  ,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i-1,3,4) * (cord(j  ,i-1,1)-cord(j  ,i  ,1))
            fy = stress0(j  ,i-1,3,4) * (cord(j  ,i-1,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i-1,2,4) * (cord(j  ,i-1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

        endif

        ! Element (j,i). Triangles A,C,D
        if ( j.ne.nz .and. i.ne.nx ) then
            ! triangle A
            ! side 1-2
            fx = stress0(j  ,i  ,1,1) * (cord(j+1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i  ,3,1) * (cord(j+1,i  ,1)-cord(j  ,i  ,1))
            fy = stress0(j  ,i  ,3,1) * (cord(j+1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i  ,2,1) * (cord(j+1,i  ,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 3-1
            fx = stress0(j  ,i  ,1,1) * (cord(j  ,i  ,2)-cord(j  ,i+1,2)) - &
                 stress0(j  ,i  ,3,1) * (cord(j  ,i  ,1)-cord(j  ,i+1,1))
            fy = stress0(j  ,i  ,3,1) * (cord(j  ,i  ,2)-cord(j  ,i+1,2)) - &
                 stress0(j  ,i  ,2,1) * (cord(j  ,i  ,1)-cord(j  ,i+1,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

            ! triangle C
            ! side 1-2
            fx = stress0(j  ,i  ,1,3) * (cord(j+1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i  ,3,3) * (cord(j+1,i  ,1)-cord(j  ,i  ,1))
            fy = stress0(j  ,i  ,3,3) * (cord(j+1,i  ,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i  ,2,3) * (cord(j+1,i  ,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 3-1
            fx = stress0(j  ,i  ,1,3) * (cord(j  ,i  ,2)-cord(j+1,i+1,2)) - &
                 stress0(j  ,i  ,3,3) * (cord(j  ,i  ,1)-cord(j+1,i+1,1))
            fy = stress0(j  ,i  ,3,3) * (cord(j  ,i  ,2)-cord(j+1,i+1,2)) - &
                 stress0(j  ,i  ,2,3) * (cord(j  ,i  ,1)-cord(j+1,i+1,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

            ! triangle D
            ! side 1-2
            fx = stress0(j  ,i  ,1,4) * (cord(j+1,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i  ,3,4) * (cord(j+1,i+1,1)-cord(j  ,i  ,1))
            fy = stress0(j  ,i  ,3,4) * (cord(j+1,i+1,2)-cord(j  ,i  ,2)) - &
                 stress0(j  ,i  ,2,4) * (cord(j+1,i+1,1)-cord(j  ,i  ,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)
            ! side 3-1
            fx = stress0(j  ,i  ,1,4) * (cord(j  ,i  ,2)-cord(j  ,i+1,2)) - &
                 stress0(j  ,i  ,3,4) * (cord(j  ,i  ,1)-cord(j  ,i+1,1))
            fy = stress0(j  ,i  ,3,4) * (cord(j  ,i  ,2)-cord(j  ,i+1,2)) - &
                 stress0(j  ,i  ,2,4) * (cord(j  ,i  ,1)-cord(j  ,i+1,1))
            force(j,i,1) = force(j,i,1) - 0.25*fx
            force(j,i,2) = force(j,i,2) - 0.25*fy
            balance(j,i,1) = balance(j,i,1) + 0.25*abs(fx)
            balance(j,i,2) = balance(j,i,2) + 0.25*abs(fy)

        endif

        ! GRAVITY FORCE
        force(j,i,2) = force(j,i,2) - rmass(j,i)*g
        balance(j,i,2) = balance(j,i,2) + abs(rmass(j,i)*g)

  enddo
enddo
!$OMP end do

! BOUNDARY CONDITIONS
if(nyhydro.gt.0) then
    !$OMP do
    do i=1,nx

        ! pressure from water sea on top
        rho_water_g = 1030. * g
        if(i.lt.nx) then
            water_depth = 0.5*(cord(1,i+1,2)+cord(1,i,2))
        else
            water_depth = 0.5*(cord(1,i-1,2)+cord(1,i,2))
        endif

        if (water_depth.lt.0.) then ! No water (above sea level)
            if(i.eq.1) then
                press_norm_l = 0
                dlx_l = 0
                dly_l = 0
                press_norm_r = rho_water_g*((cord(1,i+1,2)+cord(1,i,2))/2.)
                dlx_r = cord(1,i+1,1)-cord(1,i  ,1)
                dly_r = cord(1,i+1,2)-cord(1,i  ,2)
            elseif(i.eq.nx) then
                press_norm_l = rho_water_g*((cord(1,i-1,2)+cord(1,i,2))/2.)
                dlx_l = cord(1,i  ,1)-cord(1,i-1,1)
                dly_l = cord(1,i  ,2)-cord(1,i-1,2)
                press_norm_r = 0
                dlx_r = 0
                dly_r = 0
            else
                press_norm_l = rho_water_g*((cord(1,i-1,2)+cord(1,i,2))/2.)
                dlx_l = cord(1,i  ,1)-cord(1,i-1,1)
                dly_l = cord(1,i  ,2)-cord(1,i-1,2)
                press_norm_r = rho_water_g*((cord(1,i+1,2)+cord(1,i,2))/2.)
                dlx_r = cord(1,i+1,1)-cord(1,i  ,1)
                dly_r = cord(1,i+1,2)-cord(1,i  ,2)
            endif
            force(1,i,1) = force(1,i,1)-0.5*press_norm_l*dly_l-0.5*press_norm_r*dly_r
            force(1,i,2) = force(1,i,2)+0.5*press_norm_l*dlx_l+0.5*press_norm_r*dlx_r
            balance(1,i,1) = 1.0d+17
        endif
    enddo
    !$OMP end do

    !$OMP do
    do i=1,nx

        ! bottom support - Archimed force (normal to the surface, shear component = 0)
        p_est = pisos + 0.5*(den(iphsub)+drosub)*g*(cord(nz,i,2)-rzbo)
        rosubg = g * (den(iphsub)+drosub) * (1-alfa(iphsub)*temp(nz,i)+beta(iphsub)*p_est)

        if(i.eq.1) then
            press_norm_l = 0
            dlx_l = 0
            dly_l = 0

            press_norm_r = pisos-rosubg*((cord(nz,i+1,2)+cord(nz,i,2))/2-rzbo)
            dlx_r = cord(nz,i+1,1)-cord(nz,i  ,1)
            dly_r = cord(nz,i+1,2)-cord(nz,i  ,2)
        elseif(i.eq.nx) then
            press_norm_l = pisos-rosubg*((cord(nz,i-1,2)+cord(nz,i,2))/2-rzbo)
            dlx_l = cord(nz,i  ,1)-cord(nz,i-1,1)
            dly_l = cord(nz,i  ,2)-cord(nz,i-1,2)

            press_norm_r = 0
            dlx_r = 0
            dly_r = 0
        else
            press_norm_l = pisos-rosubg*((cord(nz,i-1,2)+cord(nz,i,2))/2-rzbo)
            dlx_l = cord(nz,i  ,1)-cord(nz,i-1,1)
            dly_l = cord(nz,i  ,2)-cord(nz,i-1,2)

            press_norm_r = pisos-rosubg*((cord(nz,i+1,2)+cord(nz,i,2))/2-rzbo)
            dlx_r = cord(nz,i+1,1)-cord(nz,i  ,1)
            dly_r = cord(nz,i+1,2)-cord(nz,i  ,2)
        endif
            
        force(nz,i,1) = force(nz,i,1)-0.5*press_norm_l*dly_l-0.5*press_norm_r*dly_r
        force(nz,i,2) = force(nz,i,2)+0.5*press_norm_l*dlx_l+0.5*press_norm_r*dlx_r

        balance(nz,i,1) = 1.0d+17
        !write(*,*) i,pisos,force(nz,i,1),force(nz,i,2),press_norm_l,press_norm_r,dlx_l,dlx_r,dly_l,dly_r

    enddo
    !$OMP end do
endif

! Traction bc for side walls
!! LEFT WALL

do j=1,nz

!    rho_mantle_g = 1.e6 !3e8 !3300. * g
    if ((time .lt. 3.85*1e3*3.1536e7)) then
        rho_mantle_g = time*1.e8 / (1*1e6*3.1536e7)
!    else
!        rho_mantle_g = 0 !3.85*1e3*3.1536e7*1.e8 / (1*1e6*3.1536e7)
!    end if    
    
        !if () then ! For every elements in model
        if(j.eq.1) then
            press_norm_u = 0
            dlx_u = 0
            dly_u = 0
            press_norm_d = rho_mantle_g !*((cord(j+1,1,2)+cord(j,1,2))/2.)
            dlx_d = cord(j+1,1,1)-cord(j,1  ,1)
            dly_d = cord(j+1,1,2)-cord(j,1  ,2)
        elseif(j.eq.nz) then
            press_norm_u = rho_mantle_g !*((cord(j-1,1,2)+cord(j,1,2))/2.)
            dlx_u = cord(j,1  ,1)-cord(j-1,1,1)
            dly_u = cord(j,1  ,2)-cord(j-1,1,2)
            press_norm_d = 0
            dlx_d = 0
            dly_d = 0
        else
            press_norm_u = rho_mantle_g !*((cord(j-1,1,2)+cord(j,1,2))/2.)
            dlx_u = cord(j,1  ,1)-cord(j-1,1,1)
            dly_u = cord(j,1  ,2)-cord(j-1,1,2)
            press_norm_d = rho_mantle_g !*((cord(j+1,1,2)+cord(j,1,2))/2.)
            dlx_d = cord(j+1,1,1)-cord(j,1  ,1)
            dly_d = cord(j+1,1,2)-cord(j,1  ,2)
        endif
        !print *, 'left j=', j, force(j,1,1), 0.5*press_norm_u*dly_u, 0.5*press_norm_d*dly_d
!        if ((time .gt. 3.85*1.e3*3.1536e7 .and. time .gt. 3.8500001*1.e3*3.1536e7)) then
!            print *, '1) left j=', j, force(j,1,1), 'time =', time
!            print *, '2) left j=', j, force(j,1,2), 'time =', time
!        endif
        force(j,1,1) = 0.5*press_norm_u*dly_u + 0.5*press_norm_d*dly_d
        force(j,1,2) = - 0.5*press_norm_u*dlx_u - 0.5*press_norm_d*dlx_d
        balance(j,1,1) = 1.0d+17
        !print *, force(j,1,1)
        !endif
    else
        if(j.eq.1) then
            force(j,1,1) = force(j,1,1) - 2748272910.4164391
            force(j,1,2) = force(j,1,2) + 1835361956.0045931
        elseif(j.eq.2) then
            force(j,1,1) = force(j,1,1) + 6116402938.0449286
            force(j,1,2) = force(j,1,2) + 2588859961.2459397
        elseif(j.eq.3) then
            force(j,1,1) = force(j,1,1) + 17111038783.990341
            force(j,1,2) = force(j,1,2) + 3613845491.1236172
        elseif(j.eq.4) then
            force(j,1,1) = force(j,1,1) + 26702263578.216248
            force(j,1,2) = force(j,1,2) + 4365026689.4011660
        elseif(j.eq.5) then
            force(j,1,1) = force(j,1,1) + 35830619899.019089
            force(j,1,2) = force(j,1,2) + 4992805618.8505888
        elseif(j.eq.6) then
            force(j,1,1) = force(j,1,1) + 44680929366.182289
            force(j,1,2) = force(j,1,2) + 5632265763.7193165
        elseif(j.eq.7) then
            force(j,1,1) = force(j,1,1) + 53358315566.769890
            force(j,1,2) = force(j,1,2) + 6241032003.4354467
        elseif(j.eq.8) then
            force(j,1,1) = force(j,1,1) + 62137885676.840935
            force(j,1,2) = force(j,1,2) + 6878487437.8966713
        elseif(j.eq.9) then
            force(j,1,1) = force(j,1,1) + 71204361593.618759
            force(j,1,2) = force(j,1,2) + 7580835918.5424433
        elseif(j.eq.10) then
            force(j,1,1) = force(j,1,1) + 81263856729.717529
            force(j,1,2) = force(j,1,2) + 8334132644.9650574
        elseif(j.eq.11) then
            force(j,1,1) = force(j,1,1) + 93663397180.344116
            force(j,1,2) = force(j,1,2) + 7003342030.0599861
        elseif(j.eq.12) then
            force(j,1,1) = force(j,1,1) + 97447462417.585678
            force(j,1,2) = force(j,1,2) + 9630661994.6603203
        elseif(j.eq.13) then
            force(j,1,1) = force(j,1,1) + 91026839667.313629
            force(j,1,2) = force(j,1,2) + 1641438515.7779603
        elseif(j.eq.14) then
            force(j,1,1) = force(j,1,1) + 98345356917.754395
            force(j,1,2) = force(j,1,2) + 191752490.99300003
        elseif(j.eq.15) then
            force(j,1,1) = force(j,1,1) + 106439163934.14070
            force(j,1,2) = force(j,1,2) + 591707904.12245846
        elseif(j.eq.16) then
            force(j,1,1) = force(j,1,1) + 114519548987.10474
            force(j,1,2) = force(j,1,2) + 222905194.84442663
        elseif(j.eq.17) then
            force(j,1,1) = force(j,1,1) + 122720519242.91457
            force(j,1,2) = force(j,1,2) + 346513061.68616056
        elseif(j.eq.18) then
            force(j,1,1) = force(j,1,1) + 130897916136.60071
            force(j,1,2) = force(j,1,2) + 171584555.76412392
        elseif(j.eq.19) then
            force(j,1,1) = force(j,1,1) + 139134753754.36417
            force(j,1,2) = force(j,1,2) + 146237417.45102406
        elseif(j.eq.20) then
            force(j,1,1) = force(j,1,1) + 147387544263.11316
            force(j,1,2) = force(j,1,2) + 16361320.445938110
        elseif(j.eq.21) then
            force(j,1,1) = force(j,1,1) + 155696663201.50354
            force(j,1,2) = force(j,1,2) - 79670123.641295910
        elseif(j.eq.22) then
            force(j,1,1) = force(j,1,1) + 164036736889.04822
            force(j,1,2) = force(j,1,2) - 188576391.09359694
        elseif(j.eq.23) then
            force(j,1,1) = force(j,1,1) + 172392460846.37064
            force(j,1,2) = force(j,1,2) - 269383907.20983744
        elseif(j.eq.24) then
            force(j,1,1) = force(j,1,1) + 180759244882.09344
            force(j,1,2) = force(j,1,2) - 355545881.95638084
        elseif(j.eq.25) then
            force(j,1,1) = force(j,1,1) + 189177920364.40732
            force(j,1,2) = force(j,1,2) - 445475540.71091461
        elseif(j.eq.26) then
            force(j,1,1) = force(j,1,1) + 197637371946.97766
            force(j,1,2) = force(j,1,2) - 503441064.32550764
        elseif(j.eq.27) then
            force(j,1,1) = force(j,1,1) + 206108126794.22797
            force(j,1,2) = force(j,1,2) - 558405622.41003513
        elseif(j.eq.28) then
            force(j,1,1) = force(j,1,1) + 214609824114.71829
            force(j,1,2) = force(j,1,2) - 631589190.89238358
        elseif(j.eq.29) then
            force(j,1,1) = force(j,1,1) + 223089152678.70657
            force(j,1,2) = force(j,1,2) - 746910579.91093731
        elseif(j.eq.30) then
            force(j,1,1) = force(j,1,1) + 231572135593.86963
            force(j,1,2) = force(j,1,2) - 852219435.41282749
        elseif(j.eq.31) then
            force(j,1,1) = force(j,1,1) + 240084408836.58768
            force(j,1,2) = force(j,1,2) - 1003767979.0382090
        elseif(j.eq.32) then
            force(j,1,1) = force(j,1,1) + 248610325921.29648
            force(j,1,2) = force(j,1,2) - 1110593809.2604723
        elseif(j.eq.33) then
            force(j,1,1) = force(j,1,1) + 257174189806.74487
            force(j,1,2) = force(j,1,2) - 1316683926.5874944
        elseif(j.eq.34) then
            force(j,1,1) = force(j,1,1) + 265745144840.12967
            force(j,1,2) = force(j,1,2) - 1405522978.1025901
        elseif(j.eq.35) then
            force(j,1,1) = force(j,1,1) + 274380069057.13968
            force(j,1,2) = force(j,1,2) - 1714722623.3716435
        elseif(j.eq.36) then
            force(j,1,1) = force(j,1,1) + 282998578074.14081
            force(j,1,2) = force(j,1,2) - 1732695485.2171059
        elseif(j.eq.37) then
            force(j,1,1) = force(j,1,1) + 291735901524.72498
            force(j,1,2) = force(j,1,2) - 2258093168.5952892
        elseif(j.eq.38) then
            force(j,1,1) = force(j,1,1) + 300351829919.42566
            force(j,1,2) = force(j,1,2) - 2032572717.0067153
        elseif(j.eq.39) then
            force(j,1,1) = force(j,1,1) + 309214654833.23981
            force(j,1,2) = force(j,1,2) - 3130813832.2078409
        elseif(j.eq.40) then
            force(j,1,1) = force(j,1,1) + 317094463961.80096
            force(j,1,2) = force(j,1,2) - 1735267442.2393627
        elseif(j.eq.nz) then
            force(j,1,1) = force(j,1,1) + 157840450051.41693
            force(j,1,2) = force(j,1,2) - 4162869132.3101807
        end if
!force(j,1,1) = force(j,1,1) +
!force(j,1,2) = force(j,1,2) +
    balance(j,1,1) = 1.0d+17
    end if
end do

!! RIGHT WALL
!$OMP do
do j=1,nz

! bottom support - Archimed force (normal to the surface, shear component = 0)
!p_est = pisos + 0.5*(den(iphsub)+drosub)*g*(cord(nz,i,2)-rzbo)
!rosubg = g * (den(iphsub)+drosub) * (1-alfa(iphsub)*temp(nz,i)+beta(iphsub)*p_est)
!        rho_mantle_g = 1.0e6 !3e8 !3300. * g
    if ((time .lt. 3.85*1.e3*3.1536e7)) then
        rho_mantle_g = time*1.e8 / (1*1e6*3.1536e7)
        !else
        !    rho_mantle_g = 0 !3.85*1.e3*3.1536e7*1.e8 / (1*1e6*3.1536e7)
        !end if

        if(j.eq.1) then
            press_norm_u = 0
            dlx_u = 0
            dly_u = 0
            press_norm_d = rho_mantle_g !*((cord(j+1,nx,2)+cord(j,nx,2))/2.)
            dlx_d = cord(j+1,nx,1)-cord(j,nx  ,1)
            dly_d = cord(j+1,nx,2)-cord(j,nx  ,2)
        elseif(j.eq.nz) then
            press_norm_u = rho_mantle_g !*((cord(j-1,nx,2)+cord(j,nx,2))/2.)
            dlx_u = cord(j,nx  ,1)-cord(j-1,nx,1)
            dly_u = cord(j,nx  ,2)-cord(j-1,nx,2)
            press_norm_d = 0
            dlx_d = 0
            dly_d = 0
        else
            press_norm_u = rho_mantle_g !*((cord(j-1,nx,2)+cord(j,nx,2))/2.)
            dlx_u = cord(j,nx  ,1)-cord(j-1,nx,1)
            dly_u = cord(j,nx  ,2)-cord(j-1,nx,2)
            press_norm_d = rho_mantle_g !*((cord(j+1,nx,2)+cord(j,nx,2))/2.)
            dlx_d = cord(j+1,nx,1)-cord(j,nx  ,1)
            dly_d = cord(j+1,nx,2)-cord(j,nx  ,2)
        endif
        !print *, 'right j=', j, force(j,nx,1), 0.5*press_norm_u*dly_u, 0.5*press_norm_d*dly_d
        force(j,nx,1) = -0.5*press_norm_u*dly_u-0.5*press_norm_d*dly_d
        force(j,nx,2) = 0.5*press_norm_u*dlx_u+0.5*press_norm_d*dlx_d
        balance(j,nx,1) = 1.0d+17
        !print *, force(j,nx,1)
    !write(*,*) i,pisos,force(nz,i,1),force(nz,i,2),press_norm_l,press_norm_r,dlx_l,dlx_r,dly_l,dly_r
    else
        if(j.eq.1) then
            force(j,1,1) = force(j,1,1) + 2748272910.4164391
            force(j,1,2) = force(j,1,2) - 1835361956.0045931
        elseif(j.eq.2) then
            force(j,1,1) = force(j,1,1) - 6116402938.0449286
            force(j,1,2) = force(j,1,2) - 2588859961.2459397
        elseif(j.eq.3) then
            force(j,1,1) = force(j,1,1) - 17111038783.990341                 
            force(j,1,2) = force(j,1,2) - 3613845491.1236172
        elseif(j.eq.4) then
            force(j,1,1) = force(j,1,1) - 26702263578.216248
            force(j,1,2) = force(j,1,2) - 4365026689.4011660
        elseif(j.eq.5) then
            force(j,1,1) = force(j,1,1) - 35830619899.019089
            force(j,1,2) = force(j,1,2) - 4992805618.8505888
        elseif(j.eq.6) then
            force(j,1,1) = force(j,1,1) - 44680929366.182289
            force(j,1,2) = force(j,1,2) - 5632265763.7193165
        elseif(j.eq.7) then
            force(j,1,1) = force(j,1,1) - 53358315566.769890
            force(j,1,2) = force(j,1,2) - 6241032003.4354467
        elseif(j.eq.8) then
            force(j,1,1) = force(j,1,1) - 62137885676.840935
            force(j,1,2) = force(j,1,2) - 6878487437.896671
        elseif(j.eq.9) then
            force(j,1,1) = force(j,1,1) - 71204361593.618759
            force(j,1,2) = force(j,1,2) - 7580835918.5424433
        elseif(j.eq.10) then
            force(j,1,1) = force(j,1,1) - 81263856729.717529
            force(j,1,2) = force(j,1,2) - 8334132644.9650574
        elseif(j.eq.11) then
            force(j,1,1) = force(j,1,1) - 93663397180.344116
            force(j,1,2) = force(j,1,2) - 7003342030.0599861
        elseif(j.eq.12) then
            force(j,1,1) = force(j,1,1) - 97447462417.585678
            force(j,1,2) = force(j,1,2) - 9630661994.6603203
        elseif(j.eq.13) then
            force(j,1,1) = force(j,1,1) - 91026839667.313629
            force(j,1,2) = force(j,1,2) - 1641438515.7779603
        elseif(j.eq.14) then
            force(j,1,1) = force(j,1,1) - 98345356917.754395
            force(j,1,2) = force(j,1,2) - 191752490.99300003
        elseif(j.eq.15) then
            force(j,1,1) = force(j,1,1) - 106439163934.14070
            force(j,1,2) = force(j,1,2) - 591707904.12245846
        elseif(j.eq.16) then
            force(j,1,1) = force(j,1,1) - 114519548987.10474
            force(j,1,2) = force(j,1,2) - 222905194.84442663
        elseif(j.eq.17) then
            force(j,1,1) = force(j,1,1) - 122720519242.91457
            force(j,1,2) = force(j,1,2) - 346513061.68616056
        elseif(j.eq.18) then
            force(j,1,1) = force(j,1,1) - 130897916136.60071
            force(j,1,2) = force(j,1,2) - 171584555.76412392
        elseif(j.eq.19) then
            force(j,1,1) = force(j,1,1) - 139134753754.36417
            force(j,1,2) = force(j,1,2) - 146237417.45102406
        elseif(j.eq.20) then
            force(j,1,1) = force(j,1,1) - 147387544263.11316
            force(j,1,2) = force(j,1,2) - 16361320.445938110
        elseif(j.eq.21) then
            force(j,1,1) = force(j,1,1) - 155696663201.50354
            force(j,1,2) = force(j,1,2) + 79670123.641295910
        elseif(j.eq.22) then
            force(j,1,1) = force(j,1,1) - 164036736889.04822
            force(j,1,2) = force(j,1,2) + 188576391.09359694
        elseif(j.eq.23) then
            force(j,1,1) = force(j,1,1) - 172392460846.37064
            force(j,1,2) = force(j,1,2) + 269383907.20983744
        elseif(j.eq.24) then
            force(j,1,1) = force(j,1,1) - 180759244882.09344
            force(j,1,2) = force(j,1,2) + 355545881.95638084
        elseif(j.eq.25) then
            force(j,1,1) = force(j,1,1) - 189177920364.40732
            force(j,1,2) = force(j,1,2) + 445475540.71091461
        elseif(j.eq.26) then
            force(j,1,1) = force(j,1,1) - 197637371946.97766
            force(j,1,2) = force(j,1,2) + 503441064.32550764
        elseif(j.eq.27) then
            force(j,1,1) = force(j,1,1) - 206108126794.22797
            force(j,1,2) = force(j,1,2) + 558405622.41003513
        elseif(j.eq.28) then
            force(j,1,1) = force(j,1,1) - 214609824114.71829
            force(j,1,2) = force(j,1,2) + 631589190.89238358
        elseif(j.eq.29) then
            force(j,1,1) = force(j,1,1) - 223089152678.70657
            force(j,1,2) = force(j,1,2) + 746910579.91093731
        elseif(j.eq.30) then
            force(j,1,1) = force(j,1,1) - 231572135593.86963
            force(j,1,2) = force(j,1,2) + 852219435.41282749
        elseif(j.eq.31) then
            force(j,1,1) = force(j,1,1) - 240084408836.58768
            force(j,1,2) = force(j,1,2) + 1003767979.0382090
        elseif(j.eq.32) then
            force(j,1,1) = force(j,1,1) - 248610325921.29648
            force(j,1,2) = force(j,1,2) + 1110593809.2604723
        elseif(j.eq.33) then
            force(j,1,1) = force(j,1,1) - 257174189806.74487
            force(j,1,2) = force(j,1,2) + 1316683926.5874944
        elseif(j.eq.34) then
            force(j,1,1) = force(j,1,1) - 265745144840.12967
            force(j,1,2) = force(j,1,2) + 1405522978.1025901
        elseif(j.eq.35) then
            force(j,1,1) = force(j,1,1) - 274380069057.13968
            force(j,1,2) = force(j,1,2) + 1714722623.3716435
        elseif(j.eq.36) then
            force(j,1,1) = force(j,1,1) - 282998578074.14081
            force(j,1,2) = force(j,1,2) + 1732695485.2171059
        elseif(j.eq.37) then
            force(j,1,1) = force(j,1,1) - 291735901524.72498
            force(j,1,2) = force(j,1,2) + 2258093168.5952892
        elseif(j.eq.38) then
            force(j,1,1) = force(j,1,1) - 300351829919.42566
            force(j,1,2) = force(j,1,2) + 2032572717.0067153
        elseif(j.eq.39) then
            force(j,1,1) = force(j,1,1) - 309214654833.23981
            force(j,1,2) = force(j,1,2) + 3130813832.2078409
        elseif(j.eq.40) then
            force(j,1,1) = force(j,1,1) - 317094463961.80096
            force(j,1,2) = force(j,1,2) + 1735267442.2393627
        elseif(j.eq.nz) then
            force(j,1,1) = force(j,1,1) - 157840450051.41693
            force(j,1,2) = force(j,1,2) + 4162869132.3101807
        end if
    end if
    balance(j,1,1) = 1.0d+17
end do            
!$OMP end do
!endif



!$OMP do reduction(max:boff)
do i=1,nx
    do j=1,nz

        ! BALANCE-OFF
        if( iand(ncod(j,i,1),1).eq.1 .or. j.le.n_boff_cutoff ) then
            balance(j,i,1) = 0
        else
            balance(j,i,1) = abs(force(j,i,1)) / (balance(j,i,1) + 1.0d-9)
        endif

        if( iand(ncod(j,i,2),2).eq.2 .or. j.le.n_boff_cutoff ) then
            balance(j,i,2) = 0
        else
            balance(j,i,2) = abs(force(j,i,2)) / (balance(j,i,2) + 1.0d-9)
        endif

        ! DAMPING
        if( iand(ncod(j,i,1),1).ne.1 .and. abs(vel(j,i,1)).gt.1.0d-13 ) then
            force(j,i,1) = force(j,i,1) - demf*sign(force(j,i,1),vel(j,i,1))
        endif

        if( iand(ncod(j,i,2),2).ne.2 .and. abs(vel(j,i,2)).gt.1.0d-13 ) then
            force(j,i,2) = force(j,i,2) - demf*sign(force(j,i,2),vel(j,i,2))
        endif

        ! VELOCITIES FROM FORCES
        if( ncod(j,i,1) .eq. 1 ) then
            vel(j,i,1) = bc(j,i,1) 
        else
            vel(j,i,1) = vel(j,i,1) + dt*force(j,i,1)/(amass(j,i)*drat*drat)
        endif
        if( ncod(j,i,2) .eq. 1 ) then
            vel(j,i,2) = bc(j,i,2)
        else
            vel(j,i,2) = vel(j,i,2) + dt*force(j,i,2)/(amass(j,i)*drat*drat)
        endif
        ! MAX balance-off
        boff = max(boff,balance(j,i,1))
        boff = max(boff,balance(j,i,2))

    end do
end do
!$OMP end do
!$OMP end parallel
! Prestress to form the topo when density differences are present WITHOUT PUSHING OR PULLING!
if (i_prestress.eq.1.and.time.lt.600.e3*sec_year) then
    do k = 1,2
        do i = 1, nx
            vel(nz,i,k) = 0.
        enddo
        do j = 1, nz
            vel(j,1,k) = 0.
            vel(j,nx,k) = 0.
        enddo
    enddo
endif
return
end subroutine fl_node
