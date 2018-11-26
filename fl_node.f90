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
!$OMP do
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
            force(j,1,1) = force(j,1,1) - 2772511879.1606851
            force(j,1,2) = force(j,1,2) + 1851107371.9701347
        elseif(j.eq.2) then
            force(j,1,1) = force(j,1,1) + 6147369066.5476007
            force(j,1,2) = force(j,1,2) + 2618599619.2106161
        elseif(j.eq.3) then
            force(j,1,1) = force(j,1,1) + 17150625223.323690
            force(j,1,2) = force(j,1,2) + 3647156761.7931528
        elseif(j.eq.4) then
            force(j,1,1) = force(j,1,1) + 26731113653.777958
            force(j,1,2) = force(j,1,2) + 4396740221.0536909
        elseif(j.eq.5) then
            force(j,1,1) = force(j,1,1) + 35865872439.688438
            force(j,1,2) = force(j,1,2) + 5021429434.3409786
        elseif(j.eq.6) then
            force(j,1,1) = force(j,1,1) + 44730837149.789131
            force(j,1,2) = force(j,1,2) + 5663493635.7465458
        elseif(j.eq.7) then
            force(j,1,1) = force(j,1,1) + 53397359248.944855
            force(j,1,2) = force(j,1,2) + 6282380903.8103409
        elseif(j.eq.8) then
            force(j,1,1) = force(j,1,1) + 62150368442.934250
            force(j,1,2) = force(j,1,2) + 6912610827.9923925
        elseif(j.eq.9) then
            force(j,1,1) = force(j,1,1) + 71217263617.187027
            force(j,1,2) = force(j,1,2) + 7608653457.5902281
        elseif(j.eq.10) then
            force(j,1,1) = force(j,1,1) + 81275792349.752380
            force(j,1,2) = force(j,1,2) + 8357471523.2074280
        elseif(j.eq.11) then
            force(j,1,1) = force(j,1,1) + 93672545752.985382
            force(j,1,2) = force(j,1,2) + 7025564499.6591930
        elseif(j.eq.12) then
            force(j,1,1) = force(j,1,1) + 97463648910.896515
            force(j,1,2) = force(j,1,2) + 9639521856.7541275
        elseif(j.eq.13) then
            force(j,1,1) = force(j,1,1) + 91043637396.485153
            force(j,1,2) = force(j,1,2) + 1644097760.1884775
        elseif(j.eq.14) then
            force(j,1,1) = force(j,1,1) + 98362236819.562576
            force(j,1,2) = force(j,1,2) + 194190868.87531996
        elseif(j.eq.15) then
            force(j,1,1) = force(j,1,1) + 106455940801.41135
            force(j,1,2) = force(j,1,2) + 594740167.14225292
        elseif(j.wq.16) then    
            force(j,1,1) = force(j,1,1) + 114535820366.86194
            force(j,1,2) = force(j,1,2) + 226148080.39129925
        elseif(j.eq.17) then
            force(j,1,1) = force(j,1,1) + 122736902212.65063
            force(j,1,2) = force(j,1,2) + 351491455.52102232
        elseif(j.eq.18) then
            force(j,1,1) = force(j,1,1) + 130914291840.54945
            force(j,1,2) = force(j,1,2) + 176772845.18415451
        elseif(j.eq.19) then
            force(j,1,1) = force(j,1,1) + 139149937354.01999
            force(j,1,2) = force(j,1,2) + 150698120.64690018
        elseif(j.eq.20) then
            force(j,1,1) = force(j,1,1) + 147401539907.14032
            force(j,1,2) = force(j,1,2) + 19618379.133758545
        elseif(j.eq.21) then
            force(j,1,1) = force(j,1,1) + 155709091661.09506 
            force(j,1,2) = force(j,1,2) - 76738594.925148964
        elseif(j.eq.22) then
            force(j,1,1) = force(j,1,1) + 164046451117.19232
            force(j,1,2) = force(j,1,2) - 185849507.55278587
        elseif(j.eq.23) then
            force(j,1,1) = force(j,1,1) + 172400540410.13101
            force(j,1,2) = force(j,1,2) - 268515692.66992140
        elseif(j.eq.24) then
            force(j,1,1) = force(j,1,1) + 180766329317.52972
            force(j,1,2) = force(j,1,2) - 353956758.14543915
        elseif(j.eq.25) then
            force(j,1,1) = force(j,1,1) + 189183494633.77896
            force(j,1,2) = force(j,1,2) - 444742589.74883461
        elseif(j.eq.26) then
            force(j,1,1) = force(j,1,1) + 197642112456.51712
            force(j,1,2) = force(j,1,2) - 502732317.71121788
        elseif(j.eq.27) then
            force(j,1,1) = force(j,1,1) + 206111904342.74066
            force(j,1,2) = force(j,1,2) - 557888893.39622593
        elseif(j.eq.28) then
            force(j,1,1) = force(j,1,1) + 214613038077.78442
            force(j,1,2) = force(j,1,2) - 631045736.93895864
        elseif(j.eq.29) then
            force(j,1,1) = force(j,1,1) + 223091928068.71719
            force(j,1,2) = force(j,1,2) - 746424484.28039837
        elseif(j.eq.30) then
            force(j,1,1) = force(j,1,1) + 231574502477.35065
            force(j,1,2) = force(j,1,2) - 851761335.24914074
        elseif(j.eq.31) then
            force(j,1,1) = force(j,1,1) + 240086392165.11655
            force(j,1,2) = force(j,1,2) - 1003326259.0377941
        elseif(j.eq.32) then
            force(j,1,1) = force(j,1,1) + 248611935286.99850  
            force(j,1,2) = force(j,1,2) - 1110181271.7651095
        elseif(j.eq.33) then
            force(j,1,1) = force(j,1,1) + 257175435132.34756
            force(j,1,2) = force(j,1,2) - 1316312941.6318150
        elseif(j.eq.34) then
            force(j,1,1) = force(j,1,1) + 265746045694.70230
            force(j,1,2) = force(j,1,2) - 1405216781.2217941
        elseif(j.eq.35) then
            force(j,1,1) = force(j,1,1) + 274380671896.55719
            force(j,1,2) = force(j,1,2) - 1714497826.3099451
        elseif(j.eq.36) then
            force(j,1,1) = force(j,1,1) + 282998955044.36511
            force(j,1,2) = force(j,1,2) - 1732552148.5807509
        elseif(j.eq.37) then
            force(j,1,1) = force(j,1,1) + 291736127087.70844
            force(j,1,2) = force(j,1,2) - 2258011612.5836325
        elseif(j.eq.38) then
            force(j,1,1) = force(j,1,1) + 300351962482.39581
            force(j,1,2) = force(j,1,2) - 2032538047.9748321
        elseif(j.eq.39) then
            force(j,1,1) = force(j,1,1) + 309214730616.05273
            force(j,1,2) = force(j,1,2) - 3130803951.6576819
        elseif(j.eq.40) then
            force(j,1,1) = force(j,1,1) + 317094497208.67175
            force(j,1,2) = force(j,1,2) - 1735279304.1086807
        elseif(j.eq.nz) then
            force(j,1,1) = force(j,1,1) + 157840443523.02090
            force(j,1,2) = force(j,1,2) - 4162873390.3673401
        end if
!force(j,1,1) = force(j,1,1) +
!force(j,1,2) = force(j,1,2) +
    balance(j,1,1) = 1.0d+17
    end if
end do
!$OMP end do

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
            force(j,1,1) = force(j,1,1) + 3919368836.2816157
            force(j,1,2) = force(j,1,2) + 2207273655.1352592
        elseif(j.eq.2) then
            force(j,1,1) = force(j,1,1) - 4362941711.4482355
            force(j,1,2) = force(j,1,2) + 2722041939.2233171
        elseif(j.eq.3) then
            force(j,1,1) = force(j,1,1) - 15544765363.970612
            force(j,1,2) = force(j,1,2) + 3648446572.4409709
        elseif(j.eq.4) then
            force(j,1,1) = force(j,1,1) - 25337705937.219574
            force(j,1,2) = force(j,1,2) + 4386500329.2463531
        elseif(j.eq.5) then
            force(j,1,1) = force(j,1,1) - 34354865535.907543
            force(j,1,2) = force(j,1,2) + 5055818367.5982513
        elseif(j.eq.6) then
            force(j,1,1) = force(j,1,1) - 43108041180.903374
            force(j,1,2) = force(j,1,2) + 5617751395.5300255
        elseif(j.eq.7) then
            force(j,1,1) = force(j,1,1) - 51864793780.329102
            force(j,1,2) = force(j,1,2) + 6170390041.0386391
        elseif(j.eq.8) then
            force(j,1,1) = force(j,1,1) - 60696960786.750519
            force(j,1,2) = force(j,1,2) + 6825094155.1005373
        elseif(j.eq.9) then
            force(j,1,1) = force(j,1,1) - 69775784805.669022
            force(j,1,2) = force(j,1,2) + 7454434813.7655754
        elseif(j.eq.10) then
            force(j,1,1) = force(j,1,1) - 80052457776.723038
            force(j,1,2) = force(j,1,2) + 8180152177.0099554
        elseif(j.eq.11) then
            force(j,1,1) = force(j,1,1) - 92523815015.974976 
            force(j,1,2) = force(j,1,2) + 6998205903.1676826
        elseif(j.eq.12) then
            force(j,1,1) = force(j,1,1) - 96910681473.188889
            force(j,1,2) = force(j,1,2) + 9358340644.3453693
        elseif(j.eq.13) then
            force(j,1,1) = force(j,1,1) - 90941533896.161743
            force(j,1,2) = force(j,1,2) + 1627185064.7878022
        elseif(j.eq.14) then
            force(j,1,1) = force(j,1,1) - 98274745002.433594
            force(j,1,2) = force(j,1,2) + 199581699.74137735
        elseif(j.eq.15) then
            force(j,1,1) = force(j,1,1) - 106367756682.05446
            force(j,1,2) = force(j,1,2) + 583089167.93898201
        elseif(j.eq.16) then
            force(j,1,1) = force(j,1,1) - 114442586979.70863
            force(j,1,2) = force(j,1,2) + 228930512.34689569
        elseif(j.eq.17) then
            force(j,1,1) = force(j,1,1) - 122643312517.22057
            force(j,1,2) = force(j,1,2) + 336140835.00807238
        elseif(j.eq.18) then
            force(j,1,1) = force(j,1,1) - 130827983101.38046
            force(j,1,2) = force(j,1,2) + 161017544.21415186
        elseif(j.eq.19) then
            force(j,1,1) = force(j,1,1) - 139070694053.38589
            force(j,1,2) = force(j,1,2) + 125743002.27758646
        elseif(j.eq.20) then
            force(j,1,1) = force(j,1,1) - 147326648099.60791
            force(j,1,2) = force(j,1,2) - 8164114.7779426575
        elseif(j.eq.21) then
            force(j,1,1) = force(j,1,1) - 155632540306.09311 
            force(j,1,2) = force(j,1,2) - 110732340.45475149
        elseif(j.eq.22) then
            force(j,1,1) = force(j,1,1) - 163962826619.91904
            force(j,1,2) = force(j,1,2) - 221794069.01076221
        elseif(j.eq.23) then
            force(j,1,1) = force(j,1,1) - 172313695654.79288
            force(j,1,2) = force(j,1,2) - 310499162.52732563
        elseif(j.eq.24) then
            force(j,1,1) = force(j,1,1) - 180688256926.06458
            force(j,1,2) = force(j,1,2) - 388190100.23139668
        elseif(j.eq.25) then
            force(j,1,1) = force(j,1,1) - 189094465561.68448
            force(j,1,2) = force(j,1,2) - 455728553.56906319
        elseif(j.eq.26) then
            force(j,1,1) = force(j,1,1) - 197534040052.32297
            force(j,1,2) = force(j,1,2) - 512237179.93121862
        elseif(j.eq.27) then
            force(j,1,1) = force(j,1,1) - 206006287602.67157
            force(j,1,2) = force(j,1,2) - 562870273.89819956
        elseif(j.eq.28) then
            force(j,1,1) = force(j,1,1) - 214531196007.73810
            force(j,1,2) = force(j,1,2) - 629102974.47695065
        elseif(j.eq.29) then
            force(j,1,1) = force(j,1,1) - 222998652926.62656
            force(j,1,2) = force(j,1,2) - 746748373.96940136
        elseif(j.eq.30) then
            force(j,1,1) = force(j,1,1) - 231470495045.09955
            force(j,1,2) = force(j,1,2) - 856951166.78671741
        elseif(j.eq.31) then
            force(j,1,1) = force(j,1,1) - 239975003980.63934
            force(j,1,2) = force(j,1,2) - 1019335034.8649273
        elseif(j.eq.32) then
            force(j,1,1) = force(j,1,1) - 248502630156.26144
            force(j,1,2) = force(j,1,2) - 1130508198.1725221
        elseif(j.eq.33) then
            force(j,1,1) = force(j,1,1) - 257074891020.60251
            force(j,1,2) = force(j,1,2) - 1339377959.9246883
        elseif(j.eq.34) then
            force(j,1,1) = force(j,1,1) - 265658034094.64590
            force(j,1,2) = force(j,1,2) - 1427157801.7170124
        elseif(j.eq.35) then
            force(j,1,1) = force(j,1,1) - 274307711716.32956
            force(j,1,2) = force(j,1,2) - 1735814244.3973713
        elseif(j.eq.36) then
            force(j,1,1) = force(j,1,1) - 282941155379.53113
            force(j,1,2) = force(j,1,2) - 1749420566.1836853
        elseif(j.eq.37) then
            force(j,1,1) = force(j,1,1) - 291692348134.44769
            force(j,1,2) = force(j,1,2) - 2273369123.8404708
        elseif(j.eq.38) then
            force(j,1,1) = force(j,1,1) - 300319559314.79816
            force(j,1,2) = force(j,1,2) - 2043391963.2981815
        elseif(j.eq.39) then
            force(j,1,1) = force(j,1,1) - 309191664413.28156
            force(j,1,2) = force(j,1,2) - 3141653229.0652809
        elseif(j.eq.40) then
            force(j,1,1) = force(j,1,1) - 317077029155.09100
            force(j,1,2) = force(j,1,2) - 1742905002.7869158
        elseif(j.eq.nz) then
            force(j,1,1) = force(j,1,1) - 157828505360.35623
            force(j,1,2) = force(j,1,2) - 4172846153.0714111
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
