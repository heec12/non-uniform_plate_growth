!----------------------------------------------------------------
! Inertial Masses and time steps - elastic and maxwell
!---------------------------------------------------------------

!    j,i           j,i+1
!     1-------------3
!     ! \  1     4/ !
!     !  1\     / 4 !
!     !     \ /     !
!     !     / \  2  !
!     !  3/    2\   !
!     ! / 3       \ !
!     2-------------4
!   j+1,i        j+1,i+1

subroutine dt_mass

use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'


real*8, parameter :: c1d12 = 1./12.

! minimal propagation distance
dlmin = dlmin_prop()
!if(time/sec_year.le.500000.) strain_inert = 1.e-5
if (idt_scale .eq. 0) then
    ! find dt below
elseif (idt_scale.eq.1) then 
    ! choosing dt_elastic from sup.dat (non-automatic) 
    dt_elastic = dt_scale
elseif (idt_scale.eq.2) then
    ! choosing dt_elastic from tolerance conditions (automatic)
    if( vbc .eq. 0 ) then
        ! let it be 1 cm/year
        dt_elastic = dlmin*frac*strain_inert/(1./sec_year/100)
    else
!write(*,'(5F45.20)')dt_elastic,dlmin,frac,strain_inert,vbc
        dt_elastic = dlmin*frac*strain_inert/vbc
!write(*,'(F45.20)')dt_elastic
    endif
endif
vel_max = 0.
do k = 1,2
do i = 1,nx
do j = 1,nz
   vel_max = max(vel_max,abs(vel(j,i,k)))
enddo
enddo
enddo

if (idt_scale .eq. 0) then
    amass = rmass
else
    amass = 0
end if

dtmax_therm = 1.e+28
dt_maxwell = 1.e+28



do 1 i = 1,nx-1
    do 1 j = 1,nz-1

        iph     = iphase(j,i)
        pwave   = rl(iph) + 0.6666*rm(iph)  
        dens    = den(iph)
        vel_sound = dlmin*frac/dt_elastic  
        rho_inert = pwave/(vel_sound*vel_sound)  
        if (i_rey.eq.1.and.vel_max.gt.0.) then
            rho_inert2 = (xReyn*v_min)/(vel_max*abs(rzbo))
!           write(*,*) rho_inert, rho_inert2,vel_max
            if (rho_inert.gt.rho_inert2) rho_inert = rho_inert2
        endif
        ! Find the inert. density for given geometry, elas_mod and dt_scale
        ! idt_scale = 0 (dt = frac*dx_min * sqrt(dens/pwave) )
        ! idt_scale = 1 (dt is taken from sup.dat: dt = dt_scale)
        ! idt_scale = 2 (dt = frac*dx_min * tolerance/Vbc_max)
        if (idt_scale.gt.0) then 

            ! Distribution 1/3 of the inertial mass of each element to the nodes 
            ! am3=c1d12*area(j,i,ii)*pwave*(dlmax*dt_scale/frac)**2
            ! 1/12 = 1/3 * 1/2 (2 meshes) * 1/2 (1/area_num = 2 area_real) 
            am3=c1d12*rho_inert/area(j,i,1)
            amass(j  ,i  ) = amass(j  ,i  ) + am3
            amass(j+1,i  ) = amass(j+1,i  ) + am3
            amass(j  ,i+1) = amass(j  ,i+1) + am3

            am3=c1d12*rho_inert/area(j,i,2)
            amass(j+1,i  ) = amass(j+1,i  ) + am3
            amass(j+1,i+1) = amass(j+1,i+1) + am3
            amass(j  ,i+1) = amass(j  ,i+1) + am3

            am3=c1d12*rho_inert/area(j,i,3)
            amass(j  ,i  ) = amass(j  ,i  ) + am3
            amass(j+1,i  ) = amass(j+1,i  ) + am3
            amass(j+1,i+1) = amass(j+1,i+1) + am3

            am3=c1d12*rho_inert/area(j,i,4)
            amass(j  ,i  ) = amass(j  ,i  ) + am3
            amass(j+1,i+1) = amass(j+1,i+1) + am3
            amass(j  ,i+1) = amass(j  ,i+1) + am3

        else  ! idt_scale=0
            ! Find the dtime for given geometry, density and elas_mod
            dte = frac*dlmin*sqrt(dens/pwave)
            dt_elastic = min(dt_elastic,dte)
        endif 

        ! Find the maximum THERMAL time step from Stability Criterion
        ! dtmax = dxmin^2/diffusivity = dx^2/(lyamda/cp*dens)  
        diff = Eff_conduct(j,i)/den(iph)/Eff_cp(j,i)
        dtt = dlmin*dlmin/diff
        dtmax_therm =min (dtmax_therm,dtt)

        ! Calculate maxwell time step
        if (ivis_present .eq. 1) then
            !dt_m =visn(j,i)/rm(iph)*fracm
            if( (irheol(iph).eq.3 .OR. irheol(iph).eq.12) .AND. rm(iph).lt.1.e+11 ) then
                visc_cut = 1.e+10
                if( v_min .lt. visc_cut ) then
                    rmu = rm(iph) * v_min/visc_cut
                else
                    rmu = rm(iph)
                endif
                dt_m =v_min/rmu * fracm
                dt_maxwell = min (dt_m,dt_maxwell)
            endif
        endif

1 continue 

return
end


!==============================================================
!   Global (in the whole domain) calculation of the shortest side - 
!   i.e. minimal propagation distance
!   dlmin = Area/Dmax for each triangle

function dlmin_prop()
use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

dlmin = 1.e+28

do 1 i = 1,nx-1
    do 1 j = 1,nz-1

        ! side 1-2 (triangles 1 and 3)
        dl = sqrt( (cord(j+1,i  ,1)-cord(j  ,i  ,1))**2 + (cord(j+1,i  ,2)-cord(j  ,i  ,2))**2 )
        dlm = 1./(area(j,i,1)*dl)
        if( dlm.lt.dlmin ) dlmin = dlm
        dlm = 1./(area(j,i,3)*dl)
        if( dlm.lt.dlmin ) dlmin = dlm

        ! side 2-4 (triangles 2 and 3)
        dl = sqrt( (cord(j+1,i+1,1)-cord(j+1,i  ,1))**2 + (cord(j+1,i+1,2)-cord(j+1,i  ,2))**2 )
        dlm = 1./(area(j,i,2)*dl)
        if( dlm.lt.dlmin ) dlmin = dlm
        dlm = 1./(area(j,i,3)*dl)
        if( dlm.lt.dlmin ) dlmin = dlm

        ! side 4-3 (triangles 2 and 4)
        dl = sqrt( (cord(j+1,i+1,1)-cord(j  ,i+1,1))**2 + (cord(j+1,i+1,2)-cord(j  ,i+1,2))**2 )
        dlm = 1./(area(j,i,2)*dl)
        if( dlm.lt.dlmin ) dlmin = dlm
        dlm = 1./(area(j,i,4)*dl)
        if( dlm.lt.dlmin ) dlmin = dlm

        ! side 3-1 (triangles 1 and 4)
        dl = sqrt( (cord(j  ,i+1,1)-cord(j  ,i  ,1))**2 + (cord(j  ,i+1,2)-cord(j  ,i  ,2))**2 )
        dlm = 1./(area(j,i,1)*dl)
        if( dlm.lt.dlmin ) dlmin = dlm
        dlm = 1./(area(j,i,4)*dl)
        if( dlm.lt.dlmin ) dlmin = dlm

        ! diagonal 1-4 (triangles 3 and 4)
        dl = sqrt( (cord(j+1,i+1,1)-cord(j  ,i  ,1))**2 + (cord(j+1,i+1,2)-cord(j  ,i  ,2))**2 )
        dlm = 1./(area(j,i,3)*dl)
        if( dlm.lt.dlmin ) dlmin = dlm
        dlm = 1./(area(j,i,4)*dl)
        if( dlm.lt.dlmin ) dlmin = dlm

        ! diagonal 2-3 (triangles 1 and 2)
        dl = sqrt( (cord(j+1,i  ,1)-cord(j  ,i+1,1))**2 + (cord(j+1,i  ,2)-cord(j  ,i+1,2))**2 )
        dlm = 1./(area(j,i,1)*dl)
        if( dlm.lt.dlmin ) dlmin = dlm
        dlm = 1./(area(j,i,2)*dl)
        if( dlm.lt.dlmin ) dlmin = dlm

1 continue

dlmin_prop = dlmin

return
end
