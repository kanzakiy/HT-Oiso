
!**************************************************************************************************************************************
program HT1
! 1st attempt for hydrothermal convection 
! working? 10/10/2019
! from HT1_v4uw.f90 
! modifying the left boundary condition from constant temp to flux (Iyer et al. 2010)
! using sparse matrix solver instead of blas 
! from v5; adding flx check
! from v6; removing blas calc option to save memory
! from v7; removing change in q from heat calculation (delta.q = 0 from mass balance)
! from v8: trying irregular grid 
implicit none

integer(kind=4),parameter :: nx = 200, ny = 320
real(kind=8),parameter :: xmax = 30d3   !  m
real(kind=8),parameter :: ymax = 12d3    !  m 
real(kind=8),parameter :: yr2sec = 60d0*60d0*24d0*365.25d0 ! sec/yr
real(kind=8) :: rl = 1d8 ! m ridge length
real(kind=8) dx(nx), dy(ny), dt, time  ! m,m, sec, sec
real(kind=8) xe(nx+1),ye(ny+1),xm(nx),ym(ny)
real(kind=8) temp(nx,ny),pres(nx,ny), visc(nx,ny), cp(nx,ny), rho(nx,ny)
real(kind=8) vx(nx+1,ny), vy(nx,ny+1), perm(nx,ny), poro(nx,ny)
real(kind=8) vmx(nx,ny), vmy(nx,ny), tempk(nx,ny), dv(nx,ny), qvx(nx,ny), qvy(nx,ny)
real(kind=8) tempk_pre(nx,ny),ds(nx,ny),permt_x
real(kind=8) :: kappa =  3d0  ! Heat conductivity of rock  [W m-1 K-1] = [J s-1 m-1 K-1]
real(kind=8) :: cpm =  1d3  ! Heat capacity of rock [J kg-1 K-1]
real(kind=8) :: grav = -9.8d0  ! gravity [m s-2]
real(kind=8) :: rhom = 3d3 ! rock density [kg m-3]
real(kind=8) :: tempm = 1200d0 ! C intrusion temp
real(kind=8) :: temps = 2d0 ! C seawater temp
real(kind=8) :: w = 3d-2/yr2sec ! m/s spreading rate
! real(kind=8) :: w = 9d-2/yr2sec ! m/s spreading rate
! real(kind=8) :: w = 1d-2/yr2sec ! m/s spreading rate
real(kind=8) :: presw = 25d6 !  Pa assuming 2.5 km depth 
! real(kind=8) :: presw = 50d6 !  Pa assuming 5.0 km depth 
! real(kind=8) :: presw = 10d6 !  Pa assuming 1.0 km depth 
real(kind=8) :: permb = 10d0**(-16.8d0)  ! base permeability [m2] 
real(kind=8) :: permt = 10d0**(-11.8d0)  ! top permeability [m2] 
real(kind=8) :: zhalf = 300d0 ! scale depth for half decrease of permeability (see Fisher 1998) 
real(kind=8) :: zscale = 300d0 ! scale depth for half decrease of permeability (see Fisher 1998) 
integer(kind = 4) ix,iy
real(kind=8) flxt(nx,ny), flxadv(nx,ny), flxcnd(nx,ny), flxres(nx,ny), flxb(nx,ny)
! matrix solver
integer(kind = 4) nmx 
real(kind=8), allocatable :: amx(:,:), ymx(:), emx(:)
integer(kind = 4), allocatable :: ipiv(:)
integer(kind = 4) infobls 
! when using sparse matrix solver 
integer(kind=4) n, nnz
integer(kind=4), allocatable :: ai(:), ap(:) 
real(kind=8), allocatable :: ax(:), bx(:) 
real(kind=8) :: control(20)
integer(kind=4) i
real(kind=8) info(90)
integer(kind=8) numeric
integer(kind=4) status
integer(kind=8) symbolic
integer(kind=4) sys
real(kind=8), allocatable :: kai(:)
!
real(kind=8) pmpa, ps, rhol, rhov, rho_start,  a, cjth, cjtt, cp2, cv, dpdr, dpdt, g, h, s, u, eta2, rhodum, tk 
real(kind=8),parameter :: t_crit = 647.1260000001D+00
real(kind=8),parameter :: temp_crit = 600d0
real(kind=8),parameter :: y_crit = 6d3  
real(kind=8) gascon  ! a function in steam library 
integer(kind=4) row,col, it 
real(kind=8) rhof, cpf, viscf, rhoout, cpout, viscout, presout
integer(kind=4) cd(nx,ny),cd2(nx,ny),cnt, cnt2, ixp, ixn, iyp, iyn, tmpint(5),tmpint2(5)
real(kind=8) tmprl(5)  
character*100 workdir
logical :: flgskip, flgend
logical :: flgesc = .false.
! logical :: flgesc = .true.
real(kind=8) beta
! 
character*10 dumchr(3)
integer dumint(8)

call date_and_time(dumchr(1),dumchr(2),dumchr(3),dumint)

workdir = '../ht-oiso_output/'
workdir = trim(adjustl(workdir))//'perm_expexp_-16_8-11_8_zh300_spx1_200x320_irr'
workdir = trim(adjustl(workdir))//'-'//trim(adjustl(dumchr(1)))
call system ('mkdir -p '//trim(adjustl(workdir)))
workdir = trim(adjustl(workdir))//'/'


! open(unit=100,file=trim(adjustl(workdir))//'comment.txt',action='write',status='unknown')
! write(100,*) 'write some comments'
! close(100)
 
!  initial conditions 
xe(1) = 0d0
ye(1) = 0d0
beta = 1.5d0
call make_grid(  &
    nx,beta,xmax  &! input 
    ,dx         &! output
    )
call make_grid(  &
    ny,beta,ymax  &! input 
    ,dy         &! output
    )
! dx = xmax/nx
! dy = ymax/ny
! print *, dx
! print *, sum(dx)
! print *
! print *, dy
! print *, sum(dy)
! stop
do ix=1,nx
    do iy=1,ny
        ds(ix,iy)=dx(ix)*dy(iy)
    enddo
enddo 
do ix = 2,nx+1
    xe(ix) = xe(ix-1) + dx(ix-1)
enddo
do iy = 2,ny+1
    ye(iy) = ye(iy-1) + dy(iy-1)
enddo
xm(:) = 0.5d0*(xe(1:nx) + xe(2:nx+1))
ym(:) = 0.5d0*(ye(1:ny) + ye(2:ny+1))

open(unit=100,file=trim(adjustl(workdir))//'x_edges.txt',action='write',status='unknown')
do ix=1,nx+1
    write(100,*) xe(ix)
enddo 
close(100)
open(unit=100,file=trim(adjustl(workdir))//'y_edges.txt',action='write',status='unknown')
do iy=1,ny+1
    write(100,*) ye(iy)
enddo 
close(100)
open(unit=100,file=trim(adjustl(workdir))//'x.txt',action='write',status='unknown')
do ix=1,nx
    write(100,*) xm(ix)
enddo 
close(100)
open(unit=100,file=trim(adjustl(workdir))//'y.txt',action='write',status='unknown')
do iy=1,ny
    write(100,*) ym(iy)
enddo 
close(100)
! stop

perm = 5d-14 ! m2
poro = 0.05d0 
do iy = 1, ny
    do ix = 1,nx
    pres(ix,iy) = presw + ym(iy)*1d-3*10d6  ! 1km 10MPa
    temp(ix,iy) = temps + ym(iy)*0.1d0     ! 0.03 K/m geothermal gradient 
    enddo
    ! exponential decrease
    perm(:,iy) = (permt - permb)*exp(-ym(iy)/zhalf) + permb
    ! stepwise 
    if (ym(iy) < zhalf) perm(:,iy) = permt
    if (ym(iy)>= zhalf) perm(:,iy) = permb
    ! log linear
    ! perm(:,iy) = 10d0**(max(log10(permb),(log10(permb)-log10(permt))*ym(iy)/zhalf+log10(permt)))
    ! 10**exp
    ! do ix =1,nx
        ! if (xm(ix)<=30d3) permt_x = permt 
        ! if (xm(ix)>30d3) permt_x = 10d0**max(log10(permt)+2d0*(xm(ix)-30d3)/(xmax-30d3),log10(permt))
        ! perm(ix,iy) = 10d0**((log10(permt_x)-log10(permb))*exp(-ym(iy)/zhalf)+log10(permb))
    ! enddo 
    perm(:,iy) = 10d0**((log10(permt)-log10(permb))*exp(-ym(iy)/zhalf)+log10(permb))
    ! 10**tanh
    ! perm(:,iy) = 10d0**((log10(permt)-log10(permb))*0.5d0*(1d0-tanh(ym(iy)-zhalf/zscale))+log10(permb))
enddo

! pres = presw
! temp = 0d0

temp = temps

ps = 20000.0D+00

pmpa = presw/1d6
tk = temps + 273.15D+00

rhol = 0.0D+00
rhov = 0.0D+00

if ( tk < t_crit ) then
  call psat ( tk, ps, rhol, rhov )
end if

if ( pmpa > ps ) then
  rho_start = rhol
else
  rho_start = pmpa / ( gascon() * tk )
end if

call dense ( pmpa, tk, rho_start, rhodum, dpdr )
call therm ( tk, rhodum, a, cjth, cjtt, cp2, cv, dpdr, dpdt, g, h, pmpa, s, u )
call viscosity ( tk, rhodum, eta2 )

rhof = rhodum*1d3  ! kg/m3
cpf = cp2           !  [J kg-1 K-1]
viscf = eta2*1d-6  ! Pa s   (eta is  MPa sec) dynamic viscosity 


pmpa = ((ym(ny)+dy(ny))*1d-3*10d6)/1d6
tk = tempm + 273.15D+00

rhol = 0.0D+00
rhov = 0.0D+00

if ( tk < t_crit ) then
  call psat ( tk, ps, rhol, rhov )
end if

if ( pmpa > ps ) then
  rho_start = rhol
else
  rho_start = pmpa / ( gascon() * tk )
end if

call dense ( pmpa, tk, rho_start, rhodum, dpdr )
call therm ( tk, rhodum, a, cjth, cjtt, cp2, cv, dpdr, dpdt, g, h, pmpa, s, u )
call viscosity ( tk, rhodum, eta2 )
! stop
dt = 1d12  ! sec  (~3000 yr)

rhoout = rhodum*1d3  ! kg/m3
cpout = cp2           !  [J kg-1 K-1]
viscout = eta2*1d-6  ! [Pa s]
it = 0
do while (it < 1000)
! if (it > 500 .and.presw> 25d6) presw =  25d6
flgend = .false.

print *,'it',it
do iy = 1,ny
    do ix = 1, nx
        flgskip = .false.
        if (it/=0 .and. pres(ix,iy)<0d0 ) then 
            if (cd(ix,iy)/=0) then 
                print*,"!!! NEGATIVE PRESSURE !!!",pres(ix,iy)
                open(unit=100,file=trim(adjustl(workdir))//'log.txt',action='write',status='unknown')
                write(100,*)"!!! NEGATIVE PRESSURE !!!",it,pres(ix,iy)
                close(100)
                flgend=.true.
                pres(ix,iy) = 0.1d6
                exit
                ! stop
            else 
                flgskip = .true.
            endif
        endif
        if (it/=0 .and. temp(ix,iy)< 0d0 ) then
            print *, "!!! NEGATIVE TEMPERATURE !!!",temp(ix,iy)
            open(unit=100,file=trim(adjustl(workdir))//'log.txt',action='write',status='unknown')
            write(100,*)"!!! NEGATIVE TEMPERATURE !!!",it,temp(ix,iy)
            close(100)
            ! flgend=.true.
            temp(ix,iy) = temps
            ! exit 
            ! stop
        endif
        if (flgskip .or. flgend) cycle
        ! if (temp(ix,iy) > temp_crit) perm(ix,iy) = 1d-100
        ! if (ym(iy) > y_crit) perm(ix,iy) = 1d-100
        pmpa = pres(ix,iy)/1d6
        tk = temp(ix,iy) + 273.15D+00
        
        rhol = 0.0D+00
        rhov = 0.0D+00
        
        if ( tk < t_crit ) then
          call psat ( tk, ps, rhol, rhov )
        end if

        if ( pmpa > ps ) then
          rho_start = rhol
        else
          rho_start = pmpa / ( gascon() * tk )
        end if
        call dense ( pmpa, tk, rho_start, rhodum, dpdr )
        call therm ( tk, rhodum, a, cjth, cjtt, cp2, cv, dpdr, dpdt, g, h, pmpa, s, u )
        call viscosity ( tk, rhodum, eta2 )
        
        rho(ix,iy) = rhodum*1d3  ! [kg m-3]
        cp(ix,iy) = cp2           !  [J kg-1 K-1]
        visc(ix,iy) = eta2*1d-6  ! [Pa s]  ! dynamic viscosity 
    enddo
enddo

if(flgend) exit

visc = visc/rho ! converting dynamic viscosity to kinematic viscosity [m2 s-1]

cd = 0
cd2 = 0
cnt = 100
do while(cnt>0) !  classifying pixel less than critical temperature as 1 recorded on cd2 

cnt = 0
do iy = 1,ny
    do ix = 1,nx
        ixp = ix + 1
        ixn = ix - 1
        iyp = iy + 1
        iyn = iy - 1
        if (ixp > nx) ixp = nx
        if (ixn < 1) ixn = 1
        if (iyp > ny) iyp = ny 
        if (iyn < 1) iyn = 1
        
        if (temp(ix,iy)>=temp_crit.or.ym(iy)>= y_crit) cycle
        if (cd2(ix,iy) == 1) cycle
        
        if (iy == 1) then 
            cd2(ix,iy) = 1
            cnt = cnt + 1
            cycle 
        else 
            if (cd2(ix,iyn) == 1 .or. &
                cd2(ix,iyp) == 1 .or. &
                cd2(ixp,iy) == 1 .or. &
                cd2(ixn,iy) == 1) then
                cd2(ix,iy) = 1
                cnt = cnt + 1
                cycle 
            endif 
        endif
    enddo
enddo
enddo

cnt = 0
cnt2 = 0
do iy = 1,ny
    do ix = 1, nx
        ixp = ix + 1
        ixn = ix - 1
        iyp = iy + 1
        iyn = iy - 1
        if (ixp > nx) ixp = nx
        if (ixn < 1) ixn = 1
        if (iyp > ny) iyp = ny 
        if (iyn < 1) iyn = 1
        
        if (cd2(ix,iy) == 0) cycle 
        if (cd2(ix,iy) == 1) then 
            cnt = cnt + 1
            cd(ix,iy) = cnt
            
            cnt2 = cnt2 + 1 
            if (cd2(ixp,iy) == 1 .and. cd(ixp,iy) /= cnt) cnt2 = cnt2 + 1
            if (cd2(ixn,iy) == 1 .and. cd(ixn,iy) /= cnt) cnt2 = cnt2 + 1
            if (cd2(ix,iyp) == 1 .and. cd(ix,iyp) /= cnt) cnt2 = cnt2 + 1
            if (cd2(ix,iyn) == 1 .and. cd(ix,iyn) /= cnt) cnt2 = cnt2 + 1
        endif 
    enddo
enddo
print*,cnt   

nmx = cnt
n = nmx
nnz = cnt2

if (allocated(ai)) deallocate(ai)
if (allocated(ap)) deallocate(ap)
if (allocated(ax)) deallocate(ax)
if (allocated(bx)) deallocate(bx)
if (allocated(kai)) deallocate(kai)

allocate(ai(nnz))
allocate(ap(n+1))
allocate(ax(nnz))
allocate(bx(n))
allocate(kai(n))

ai = 0
ap = 0
ax = 0d0
bx = 0d0
kai = 0d0

cnt = 0
cnt2 = 0

ap(1) = 0


do iy = 1, ny
    do ix = 1, nx 
        if (cd2(ix,iy)==0) cycle
        row = cd(ix,iy)
        
        tmpint =0
        tmprl = 0d0
        cnt2 = 1
        
        tmpint(1) = row
        
        if (iy == 1) then 
            if (ix == 1) then 
                bx(row) = & 
                    + perm(ix,iy)/((visc(ix,iy)+viscf)*0.5d0)*(grav*(rho(ix,iy)+rhof)*0.5d0 & 
                    - presw/dy(iy))/dy(iy)
                tmprl(1) = &
                    + perm(ix,iy)/((visc(ix,iy)+viscf)*0.5d0)*(1d0/dy(iy))/dy(iy)
                
                if (cd2(ix+1,iy) == 1) then 
                    col = cd(ix+1,iy) 
                    tmprl(1) = tmprl(1) &
                        -(perm(ix,iy)+perm(ix+1,iy))/(visc(ix,iy)+visc(ix+1,iy))*(-1d0)  &
                            /(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix)
                    
                    cnt2 = cnt2+1
                    tmprl(cnt2) = (perm(ix,iy)+perm(ix+1,iy))/(visc(ix,iy)+visc(ix+1,iy))*(-1d0)  &
                        /(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix+1)
                    tmpint(cnt2) = col
                    
                endif
                if (cd2(ix,iy+1) == 1) then 
                    col = cd(ix,iy+1) 
                    tmprl(1) = tmprl(1)  &
                        -(perm(ix,iy)+perm(ix,iy+1))/(visc(ix,iy)+visc(ix,iy+1))*(-1d0)  &
                            /(0.5d0*(dy(iy)+dy(iy+1)))/dy(iy)
                    bx(row) = bx(row) &
                        -(perm(ix,iy)+perm(ix,iy+1))/(visc(ix,iy)+visc(ix,iy+1))*grav*(rho(ix,iy+1)+rho(ix,iy))*0.5d0  &
                            /dy(iy)
                        
                    cnt2 = cnt2 + 1
                    tmprl(cnt2)  = (perm(ix,iy)+perm(ix,iy+1))/(visc(ix,iy)+visc(ix,iy+1))*(-1d0)  &
                        /(0.5d0*(dy(iy)+dy(iy+1)))/dy(iy+1)
                    tmpint(cnt2) = col

                endif
                
            else if (ix == nx) then 
                bx(row) = & 
                    + perm(ix,iy)/((visc(ix,iy)+viscf)*0.5d0)*(grav*(rho(ix,iy)+rhof)*0.5d0 & 
                    - presw/dy(iy))/dy(iy)
                tmprl(1) = &
                    + perm(ix,iy)/((visc(ix,iy)+viscf)*0.5d0)*(1d0/dy(iy))/dy(iy)
                
                if (cd2(ix-1,iy) == 1) then 
                    col = cd(ix-1,iy) 
                    tmprl(1) = tmprl(1)  &
                        + (perm(ix,iy)+perm(ix-1,iy))/(visc(ix,iy)+visc(ix-1,iy))*(1d0)  &
                            /(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix)
                    
                    cnt2 = cnt2 + 1
                    tmprl(cnt2) = -(perm(ix,iy)+perm(ix-1,iy))/(visc(ix,iy)+visc(ix-1,iy))*(1d0)  &
                        /(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix-1)
                    tmpint(cnt2) = col 
                endif
                if (cd2(ix,iy+1) == 1) then 
                    col = cd(ix,iy+1) 
                    tmprl(1) = tmprl(1)  &
                        -(perm(ix,iy)+perm(ix,iy+1))/(visc(ix,iy)+visc(ix,iy+1))*(-1d0)  &
                            /(0.5d0*(dy(iy)+dy(iy+1)))/dy(iy)
                    bx(row) = bx(row)  &
                        -(perm(ix,iy)+perm(ix,iy+1))/(visc(ix,iy)+visc(ix,iy+1))*grav*(rho(ix,iy)+rho(ix,iy+1))*0.5d0  &
                            /dy(iy)
                    
                    cnt2 = cnt2 + 1
                    tmprl(cnt2) = (perm(ix,iy)+perm(ix,iy+1))/(visc(ix,iy)+visc(ix,iy+1))*(-1d0)  &
                        /(0.5d0*(dy(iy)+dy(iy+1)))/dy(iy+1)
                    tmpint(cnt2) = col 
                    
                endif
            else 
                bx(row) = & 
                    + perm(ix,iy)/((visc(ix,iy)+viscf)*0.5d0)*(grav*(rho(ix,iy)+rhof)*0.5d0 & 
                    - presw/dy(iy))/dy(iy)
                tmprl(1) = &
                    + perm(ix,iy)/((visc(ix,iy)+viscf)*0.5d0)*(1d0/dy(iy))/dy(iy)
                
                if (cd2(ix-1,iy) == 1) then 
                    col = cd(ix-1,iy) 
                    tmprl(1) = tmprl(1) &
                        + (perm(ix,iy)+perm(ix-1,iy))/(visc(ix,iy)+visc(ix-1,iy))*(1d0)  &
                            /(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix)
                    
                    cnt2 = cnt2 + 1
                    tmprl(cnt2) = -(perm(ix,iy)+perm(ix-1,iy))/(visc(ix,iy)+visc(ix-1,iy))*(1d0)  &
                        /(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix-1)
                    tmpint(cnt2) = col 
                endif
                if (cd2(ix+1,iy) == 1) then 
                    col = cd(ix+1,iy) 
                    tmprl(1) = tmprl(1) &
                        -(perm(ix,iy)+perm(ix+1,iy))/(visc(ix,iy)+visc(ix+1,iy))*(-1d0)  &
                            /(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix)
                    
                    cnt2 = cnt2 + 1
                    tmprl(cnt2) = (perm(ix,iy)+perm(ix+1,iy))/(visc(ix,iy)+visc(ix+1,iy))*(-1d0)  &
                        /(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix+1)
                    tmpint(cnt2) = col 
                endif
                if (cd2(ix,iy+1) == 1) then 
                    col = cd(ix,iy+1) 
                    tmprl(1) = tmprl(1) &
                        -(perm(ix,iy)+perm(ix,iy+1))/(visc(ix,iy)+visc(ix,iy+1))*(-1d0)  &
                            /(0.5d0*(dy(iy)+dy(iy+1)))/dy(iy)
                    bx(row) = bx(row) &
                        -(perm(ix,iy)+perm(ix,iy+1))/(visc(ix,iy)+visc(ix,iy+1))*grav*(rho(ix,iy)+rho(ix,iy+1))*0.5d0  &
                            /dy(iy)
                    
                    cnt2 = cnt2 + 1
                    tmprl(cnt2) = (perm(ix,iy)+perm(ix,iy+1))/(visc(ix,iy)+visc(ix,iy+1))*(-1d0)  &
                        /(0.5d0*(dy(iy)+dy(iy+1)))/dy(iy+1)
                    tmpint(cnt2) = col 
                endif
            endif 
        else if (iy == ny) then 
            if (ix == 1) then 
                bx(row) = 0d0
                if (cd2(ix+1,iy) == 1) then 
                    col = cd(ix+1,iy) 
                    tmprl(1) = tmprl(1) &
                        -(perm(ix,iy)+perm(ix+1,iy))/(visc(ix,iy)+visc(ix+1,iy))*(-1d0)/(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix)
                    
                    cnt2 = cnt2 + 1
                    tmprl(cnt2) = (perm(ix,iy)+perm(ix+1,iy))/(visc(ix,iy)+visc(ix+1,iy))*(-1d0)/(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix+1)
                    tmpint(cnt2) = col 
                endif
                if (cd2(ix,iy-1) == 1) then 
                    col = cd(ix,iy-1) 
                    tmprl(1) = tmprl(1) &
                        + (perm(ix,iy)+perm(ix,iy-1))/(visc(ix,iy)+visc(ix,iy-1))*(1d0)/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)
                    bx(row) = bx(row) &
                        + (perm(ix,iy)+perm(ix,iy-1))/(visc(ix,iy)+visc(ix,iy-1))*grav*(rho(ix,iy)+rho(ix,iy-1))*0.5d0/dy(iy)
                    
                    cnt2 = cnt2 + 1
                    tmprl(cnt2) = -(perm(ix,iy)+perm(ix,iy-1))/(visc(ix,iy)+visc(ix,iy-1))*(1d0)/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy-1)
                    tmpint(cnt2) = col 
                endif
            else if (ix == nx) then 
                bx(row) = 0d0
                if (cd2(ix-1,iy) == 1) then 
                    col = cd(ix-1,iy) 
                    tmprl(1) = tmprl(1) &
                        + (perm(ix,iy)+perm(ix-1,iy))/(visc(ix,iy)+visc(ix-1,iy))*(1d0)/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix)
                    
                    cnt2 = cnt2 + 1
                    tmprl(cnt2) = -(perm(ix,iy)+perm(ix-1,iy))/(visc(ix,iy)+visc(ix-1,iy))*(1d0)/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix-1)
                    tmpint(cnt2) = col 
                endif
                if (cd2(ix,iy-1) == 1) then 
                    col = cd(ix,iy-1) 
                    tmprl(1) = tmprl(1) &
                        + (perm(ix,iy)+perm(ix,iy-1))/(visc(ix,iy)+visc(ix,iy-1))*(1d0)/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)
                    bx(row) = bx(row) &
                        + (perm(ix,iy)+perm(ix,iy-1))/(visc(ix,iy)+visc(ix,iy-1))*grav*(rho(ix,iy)+rho(ix,iy-1))*0.5d0/dy(iy)
                    
                    cnt2 = cnt2 + 1
                    tmprl(cnt2) = (perm(ix,iy)+perm(ix,iy-1))/(visc(ix,iy)+visc(ix,iy-1))*(-1d0)/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy-1)
                    tmpint(cnt2) = col 
                endif
            else 
                bx(row) = 0d0
                if (cd2(ix-1,iy) == 1) then 
                    col = cd(ix-1,iy) 
                    tmprl(1) = tmprl(1) &
                        + (perm(ix,iy)+perm(ix-1,iy))/(visc(ix,iy)+visc(ix-1,iy))*(1d0)/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix)
                    
                    cnt2 = cnt2 + 1
                    tmprl(cnt2) = -(perm(ix,iy)+perm(ix-1,iy))/(visc(ix,iy)+visc(ix-1,iy))*(1d0)/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix-1)
                    tmpint(cnt2) = col 
                endif
                if (cd2(ix+1,iy) == 1) then 
                    col = cd(ix+1,iy) 
                    tmprl(1) = tmprl(1) &
                        -(perm(ix,iy)+perm(ix+1,iy))/(visc(ix,iy)+visc(ix+1,iy))*(-1d0)/(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix)
                    
                    cnt2 = cnt2 + 1
                    tmprl(cnt2) = (perm(ix,iy)+perm(ix+1,iy))/(visc(ix,iy)+visc(ix+1,iy))*(-1d0)/(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix+1)
                    tmpint(cnt2) = col 
                endif
                if (cd2(ix,iy-1) == 1) then 
                    col = cd(ix,iy-1) 
                    tmprl(1) = tmprl(1) & 
                        + (perm(ix,iy)+perm(ix,iy-1))/(visc(ix,iy)+visc(ix,iy-1))*(1d0)/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)
                    bx(row) = bx(row) &
                        + (perm(ix,iy)+perm(ix,iy-1))/(visc(ix,iy)+visc(ix,iy-1))*grav*(rho(ix,iy)+rho(ix,iy-1))*0.5d0/dy(iy)
                    
                    cnt2 = cnt2 + 1
                    tmprl(cnt2) = -(perm(ix,iy)+perm(ix,iy-1))/(visc(ix,iy)+visc(ix,iy-1))*(1d0)/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy-1)
                    tmpint(cnt2) = col 
                endif
            endif 
        else 
            if (ix == 1) then 
                bx(row) = 0d0
                if (cd2(ix+1,iy) == 1) then 
                    col = cd(ix+1,iy) 
                    tmprl(1) = tmprl(1) &
                        -(perm(ix,iy)+perm(ix+1,iy))/(visc(ix,iy)+visc(ix+1,iy))*(-1d0)/(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix)
                    
                    cnt2 = cnt2 + 1
                    tmprl(cnt2) = (perm(ix,iy)+perm(ix+1,iy))/(visc(ix,iy)+visc(ix+1,iy))*(-1d0)/(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix+1)
                    tmpint(cnt2) = col 
                endif
                if (cd2(ix,iy-1) == 1) then 
                    col = cd(ix,iy-1) 
                    tmprl(1) = tmprl(1) & 
                        + (perm(ix,iy)+perm(ix,iy-1))/(visc(ix,iy)+visc(ix,iy-1))*(1d0)/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)
                    bx(row) = bx(row) &
                        + (perm(ix,iy)+perm(ix,iy-1))/(visc(ix,iy)+visc(ix,iy-1))*grav*(rho(ix,iy)+rho(ix,iy-1))*0.5d0/dy(iy)
                    
                    cnt2 = cnt2 + 1
                    tmprl(cnt2) = -(perm(ix,iy)+perm(ix,iy-1))/(visc(ix,iy)+visc(ix,iy-1))*(1d0)/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy-1)
                    tmpint(cnt2) = col 
                endif
                if (cd2(ix,iy+1) == 1) then 
                    col = cd(ix,iy+1) 
                    tmprl(1) = tmprl(1) &
                        -(perm(ix,iy)+perm(ix,iy+1))/(visc(ix,iy)+visc(Ix,iy+1))*(-1d0)/(0.5d0*(dy(iy)+dy(iy+1)))/dy(iy)
                    bx(row) = bx(row) &
                        - (perm(ix,iy)+perm(ix,iy+1))/(visc(ix,iy)+visc(ix,iy+1))*grav*(rho(ix,iy)+rho(ix,iy+1))*0.5d0/dy(iy)
                    
                    cnt2 = cnt2 + 1
                    tmprl(cnt2) = (perm(ix,iy)+perm(ix,iy+1))/(visc(ix,iy)+visc(Ix,iy+1))*(-1d0)/(0.5d0*(dy(iy)+dy(iy+1)))/dy(iy+1)
                    tmpint(cnt2) = col 
                endif
            else if (ix == nx) then 
                bx(row) = 0d0
                if (cd2(ix-1,iy) == 1) then 
                    col = cd(ix-1,iy) 
                    tmprl(1) = tmprl(1) &
                        +(perm(ix,iy)+perm(ix-1,iy))/(visc(ix,iy)+visc(ix-1,iy))*(1d0)/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix)
                    
                    cnt2 = cnt2 + 1
                    tmprl(cnt2) = -(perm(ix,iy)+perm(ix-1,iy))/(visc(ix,iy)+visc(ix-1,iy))*(1d0)/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix-1)
                    tmpint(cnt2) = col 
                endif
                if (cd2(ix,iy-1) == 1) then 
                    col = cd(ix,iy-1) 
                    tmprl(1) = tmprl(1) &
                        + (perm(ix,iy)+perm(ix,iy-1))/(visc(ix,iy)+visc(ix,iy-1))*(1d0)/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)
                    bx(row) = bx(row) &
                        + (perm(ix,iy)+perm(ix,iy-1))/(visc(ix,iy)+visc(ix,iy-1))*grav*(rho(ix,iy)+rho(ix,iy-1))*0.5d0/dy(iy)
                    
                    cnt2 = cnt2 + 1
                    tmprl(cnt2) = -(perm(ix,iy)+perm(ix,iy-1))/(visc(ix,iy)+visc(ix,iy-1))*(1d0)/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy-1)
                    tmpint(cnt2) = col 
                endif
                if (cd2(ix,iy+1) == 1) then 
                    col = cd(ix,iy+1) 
                    tmprl(1) = tmprl(1) &
                        -(perm(ix,iy)+perm(ix,iy+1))/(visc(ix,iy)+visc(ix,iy+1))*(-1d0)/(0.5d0*(dy(iy)+dy(iy+1)))/dy(iy)
                    bx(row) = bx(row) &
                        - (perm(ix,iy)+perm(ix,iy+1))/(visc(ix,iy)+visc(ix,iy+1))*grav*(rho(ix,iy)+rho(ix,iy+1))*0.5d0/dy(iy)
                    
                    cnt2 = cnt2 + 1
                    tmprl(cnt2) = (perm(ix,iy)+perm(ix,iy+1))/(visc(ix,iy)+visc(ix,iy+1))*(-1d0)/(0.5d0*(dy(iy)+dy(iy+1)))/dy(iy+1)
                    tmpint(cnt2) = col 
                endif
            else 
                bx(row) = 0d0
                if (cd2(ix-1,iy) == 1) then 
                    col = cd(ix-1,iy) 
                    tmprl(1) = (perm(ix,iy)+perm(ix-1,iy))/(visc(ix,iy)+visc(ix-1,iy))*(1d0)/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix)
                    
                    cnt2 = cnt2 + 1
                    tmprl(cnt2) = -(perm(ix,iy)+perm(ix-1,iy))/(visc(ix,iy)+visc(ix-1,iy))*(1d0)/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix-1)
                    tmpint(cnt2) = col 
                endif
                if (cd2(ix+1,iy) == 1) then 
                    col = cd(ix+1,iy) 
                    tmprl(1) = tmprl(1) &
                        -(perm(ix,iy)+perm(ix+1,iy))/(visc(ix,iy)+visc(ix+1,iy))*(-1d0)/(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix)
                    
                    cnt2 = cnt2 + 1
                    tmprl(cnt2) = (perm(ix,iy)+perm(ix+1,iy))/(visc(ix,iy)+visc(ix+1,iy))*(-1d0)/(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix+1)
                    tmpint(cnt2) = col 
                endif
                if (cd2(ix,iy-1) == 1) then 
                    col = cd(ix,iy-1) 
                    tmprl(1) = tmprl(1) &
                        + (perm(ix,iy)+perm(ix,iy-1))/(visc(ix,iy)+visc(ix,iy-1))*(1d0)/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)
                    bx(row) = bx(row) &
                        + (perm(ix,iy)+perm(ix,iy-1))/(visc(ix,iy)+visc(ix,iy-1))*grav*(rho(ix,iy)+rho(ix,iy-1))*0.5d0/dy(iy)
                    
                    cnt2 = cnt2 + 1
                    tmprl(cnt2) = -(perm(ix,iy)+perm(ix,iy-1))/(visc(ix,iy)+visc(ix,iy-1))*(1d0)/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy-1)
                    tmpint(cnt2) = col 
                endif
                if (cd2(ix,iy+1) == 1) then 
                    col = cd(ix,iy+1) 
                    tmprl(1) = tmprl(1) &
                        -(perm(ix,iy)+perm(ix,iy+1))/(visc(ix,iy)+visc(ix,iy+1))*(-1d0)/(0.5d0*(dy(iy)+dy(iy+1)))/dy(iy)
                    bx(row) = bx(row) &
                        - (perm(ix,iy)+perm(ix,iy+1))/(visc(ix,iy)+visc(ix,iy+1))*grav*(rho(ix,iy)+rho(ix,iy+1))*0.5d0/dy(iy)
                    
                    cnt2 = cnt2 + 1
                    tmprl(cnt2) = (perm(ix,iy)+perm(ix,iy+1))/(visc(ix,iy)+visc(ix,iy+1))*(-1d0)/(0.5d0*(dy(iy)+dy(iy+1)))/dy(iy+1)
                    tmpint(cnt2) = col 
                endif
            endif 
        endif
        
        ! print*,cnt2
        if (cnt2 < 1 .or. cnt2 > 5) then 
            print *,"FATAL ERROR"
            stop
        endif
        tmpint2 = 0 
        ! print*,'here'
        call heapsort2(cnt2, tmpint(1:cnt2),tmpint2(1:cnt2))
        do i = 1,cnt2 
            ai(cnt+i)=tmpint(i) - 1
            ax(cnt+i)=tmprl(tmpint2(i))
        enddo
        cnt = cnt + cnt2 
        ap(row+1) = cnt
        
    enddo
enddo

bx = - bx


! call DGESV(nmx,int(1),amx,nmx,IPIV,ymx,nmx,infobls) 


! solving matrix with UMFPACK (following is pasted from umfpack_simple.f90)
! Set the default control parameters.
call umf4def( control )
! From the matrix data, create the symbolic factorization information.
call umf4sym ( n, n, ap, ai, ax, symbolic, control, info )
if ( info(1) < 0.0D+00 ) then
    write ( *, * ) ''
    write ( *, *) 'UMFPACK_SIMPLE - Fatal error!'
    write ( *, * ) '  UMF4SYM returns INFO(1) = ', info(1)
    stop 1
end if
! From the symbolic factorization information, carry out the numeric factorization.
call umf4num ( ap, ai, ax, symbolic, numeric, control, info )
if ( info(1) < 0.0D+00 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'UMFPACK_SIMPLE - Fatal error!'
    write ( *, '(a,g14.6)' ) '  UMF4NUM returns INFO(1) = ', info(1)
    stop 1
end if
!  Free the memory associated with the symbolic factorization.
call umf4fsym ( symbolic )
! Solve the linear system.
sys = 0
call umf4sol ( sys, kai, bx, numeric, control, info )
if ( info(1) < 0.0D+00 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'UMFPACK_SIMPLE - Fatal error!'
    write ( *, '(a,g14.6)' ) '  UMF4SOL returns INFO(1) = ', info(1)
    stop 1
end if
! Free the memory associated with the numeric factorization.
call umf4fnum ( numeric )
!  Print the solution.
write ( *, * ) ''
write ( *, * ) '  Computed solution'
write ( *, * ) ''


do iy=1,ny
    do ix = 1,nx
        if (cd2(ix,iy) == 0) cycle
        row = cd(ix,iy)
        pres(ix,iy) = kai(row)
    enddo
enddo

do iy = 1, ny
    do ix = 1, nx
        if (iy/=1) then 
            if (cd2(ix,iy)==0 .and.cd2(ix,iy-1)==1) then 
                pres(ix,iy) = pres(ix,iy-1)-0.5d0*(dy(iy)+dy(iy-1))*grav*(rho(ix,iy)+rho(ix,iy-1))*0.5d0
            endif
        endif 
        if (iy/=ny) then 
            if (cd2(ix,iy)==0 .and.cd2(ix,iy+1)==1) then 
                pres(ix,iy) = pres(ix,iy+1)-0.5d0*(dy(iy)+dy(iy+1))*grav*(rho(ix,iy)+rho(ix,iy+1))*0.5d0
            endif
        endif 
        
        if (ix/=1) then  
            if (cd2(ix,iy)==0 .and.cd2(ix-1,iy)==1) then 
                pres(ix,iy) = pres(ix-1,iy)
            endif
        endif 
        if (ix/=nx) then  
            if (cd2(ix,iy)==0 .and.cd2(ix+1,iy)==1) then 
                pres(ix,iy) = pres(ix+1,iy)
            endif
        endif 
    enddo
enddo    

! stop
! print *, pres/1d6
vx = 0d0
vy = 0d0

do ix = 2, nx
    vx(ix,:) = -(perm(ix,:)+perm(ix-1,:))*0.5d0/((visc(ix,:)+visc(ix-1,:))*0.5d0)*(pres(ix,:)-pres(ix-1,:))  &
        /(0.5d0*(dx(ix)+dx(ix-1)))
enddo

vy(:,1) = -(perm(:,1)+perm(:,1))*0.5d0/((visc(:,1)+viscf)*0.5d0)*((pres(:,1)-presw)/dy(1) + grav*(rho(:,1)+rhof)*0.5d0)
do iy = 2, ny
    vy(:,iy) = -(perm(:,iy)+perm(:,iy-1))*0.5d0/((visc(:,iy)+visc(:,iy-1))*0.5d0)*((pres(:,iy)-pres(:,iy-1))  &
        /(0.5d0*(dy(iy)+dy(iy-1)))  &
        +grav*(rho(:,iy)+rho(:,iy-1))*0.5d0)
enddo

do iy = 1, ny
    do ix = 1, nx
        if (cd2(ix,iy)==0) then 
            vx(ix,iy) = 0d0
            vx(ix+1,iy) = 0d0
            vy(ix,iy) = 0d0
            vy(ix,iy+1) = 0d0
        endif
        ! dv(ix,iy) = (vx(ix+1,iy)-vx(ix,iy))/dx + (vy(ix,iy+1)-vy(ix,iy))/dy
    enddo
enddo    

do iy = 1, ny
    do ix = 1, nx
        dv(ix,iy) = (vx(ix+1,iy)-vx(ix,iy))/dx(ix) + (vy(ix,iy+1)-vy(ix,iy))/dy(iy)
    enddo
enddo 


qvx = 0.0d0
do iy = 1, ny
	if (iy == 1) cycle
	qvx(1,iy) = qvx(1,iy-1) + 0.5d0*(dy(iy-1)+dy(iy))*(vx(1,iy-1)+vx(2,iy)+vx(1,iy)+vx(2,iy-1))/4.0d0
	do ix = 1, nx
		if (ix == 1) cycle
		qvx(ix,iy) = qvx(ix-1,iy) - 0.5d0*(dx(ix-1)+dx(ix))*(vy(ix-1,iy) + vy(ix-1, iy+1) &
                 +vy(ix,iy) + vy(ix, iy+1))/4.0d0
	end do
end do 
qvy = 0.0d0
do ix = 1, nx
	if (ix == 1) cycle
	qvy(ix,1) = qvy(ix-1,1) - 0.5d0*(dx(ix-1)+dx(ix))*(vy(ix-1,1) + vy(ix, 2)+vy(ix,1) + vy(ix-1, 2))/4.0d0
	do iy = 1, ny
		if (iy == 1) cycle
		qvy(ix,iy) = qvy(ix,iy-1) + 0.5d0*(dy(iy-1)+dy(iy))*(vx(ix,iy-1) + vx(ix+1,iy-1)  &
                     +vx(ix,iy) + vx(ix+1,iy))/4.0d0
	end do 
end do

vmx = 0.5d0*(vx(1:nx,:)+vx(2:nx+1,:))
vmy = 0.5d0*(vy(:,1:ny)+vy(:,2:ny+1))

do iy = 1, ny
    do ix = 1, nx
        if (cd2(ix,iy)==0) then 
            vmx(ix,iy) = 0d0
            vmy(ix,iy) = 0d0
        endif
    enddo
enddo    

! returning kinematic to dynamic 
visc = visc*rho

! net flow (should be negligible for mass balance) 
print '("net flow in kg/yr:",E10.3)',sum(vy(:,1)*dx(:)*rl*yr2sec)
print '("net flow in kg/yr:",E10.3)',sum(vmy(:,1)*dx(:)*rl*yr2sec)
print '("exchanged flow in kg/yr:",E10.3)',sum(abs(vy(:,1))*dx(:)*rl*yr2sec)
print '("exchanged flow in kg/yr:",E10.3)',sum(abs(vmy(:,1))*dx(:)*rl*yr2sec)


! pause

tempk = temp + 273.15d0

nmx = nx*ny

n = nmx
nnz = 3*4 + 4*((nx-2)*2+(ny-2)*2) + 5*(nx-2)*(ny-2)

if (allocated(ai)) deallocate(ai)
if (allocated(ap)) deallocate(ap)
if (allocated(ax)) deallocate(ax)
if (allocated(bx)) deallocate(bx)
if (allocated(kai)) deallocate(kai)

allocate(ai(nnz))
allocate(ap(n+1))
allocate(ax(nnz))
allocate(bx(n))
allocate(kai(n))

ai = 0
ap = 0
ax = 0d0
bx = 0d0
kai = 0d0

cnt = 0
cnt2 = 0

ap(1) = 0
do iy = 1,ny
    do ix=1,nx
        row = (iy-1)*nx + ix
        
        tmpint =0
        tmprl = 0d0
        cnt2 = 1
        
        tmpint(1) = row
        
        if (iy == 1) then 
            if (ix==1) then 
                bx(row) = (poro(ix,iy)*rho(ix,iy)*cp(ix,iy)+(1d0-poro(ix,iy))*rhom*cpm)*(-tempk(ix,iy))/dt &
                    ! + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*cp(ix,iy)*(-(tempm + 273.15d0))/dx &
                    ! - kappa*(tempm + 273.15d0)/dx/dx  &
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)*(-(temps + 273.15d0))/dy(iy) &
                    - kappa*(temps + 273.15d0)/dy(iy)/dy(iy)  &
                    - rhom*w*(tempm+273.15d0)*cpm/dx(ix) 
                tmprl(1) = (poro(ix,iy)*rho(ix,iy)*cp(ix,iy)+(1d0-poro(ix,iy))*rhom*cpm)*(1d0)/dt &
                    ! + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(cp(ix,iy)-cpm)/dx*(1d0)  &
                    + (vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*(cp(ix+1,iy)-cp(ix,iy))/dx(ix)*(1d0)  &
                    ! + (vmx(ix,iy) - 0d0)*cp(ix,iy)/dx*(1d0)  &
                    ! + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*cp(ix,iy)/dx*(1d0)  &
                    + (vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*cp(ix,iy)/dx(ix)*(-1d0)  &
                    ! - kappa*(-2d0)/dx/dx   &
                    - kappa*(-1d0/(0.5d0*(dx(ix)+dx(ix+1))))/dx(ix)   &
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(cp(ix,iy)-cpf)/dy(iy)*(1d0) & 
                    + (vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*(cp(ix,iy+1)-cp(ix,iy))/dy(iy)*(1d0) & 
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)/dy(iy)*(1d0)  &
                    + (vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)/dy(iy)*(-1d0)  &
                    - kappa*(-1d0/dy(iy)-1d0/(0.5d0*(dy(iy)+dy(iy+1))))/dy(iy)     &
                    - rhom*w*(-1d0)*cpm/dx(ix) 
                tmprl(2) =  (vmx(ix+1,iy)+abs(vmx(ix+1,iy)))*0.50d0*cp(ix+1,iy)*(-1d0)/dx(ix+1)   &
                    - kappa*(1d0)/(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix+1)  
                tmprl(3) = (vmy(ix,iy+1)+abs(vmy(ix,iy+1)))*0.5d0*cp(ix,iy+1)/dy(iy+1)*(-1d0)  &
                    - kappa*(1d0)/(0.5d0*(dy(iy)+dy(iy+1)))/dy(iy+1) 
                     
                tmpint(2) = row+1
                tmpint(3) = row+nx
                cnt2 = cnt2 + 2
                
            else if (ix==nx) then 
                bx(row) = (poro(ix,iy)*rho(ix,iy)*cp(ix,iy)+(1d0-poro(ix,iy))*rhom*cpm)*(-tempk(ix,iy))/dt &
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)*(-(temps + 273.15d0))/dy(iy) &
                    - kappa*(temps + 273.15d0)/dy(iy)/dy(iy)  
                tmprl(1) = (poro(ix,iy)*rho(ix,iy)*cp(ix,iy)+(1d0-poro(ix,iy))*rhom*cpm)*(1d0)/dt &
                    + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(cp(ix,iy)-cp(ix-1,iy))/dx(ix)*(1d0)  &
                    ! + (vmx(ix,iy)-vmx(ix-1,iy))*cp(ix,iy)/dx*(1d0)  &
                    + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*cp(ix,iy)/dx(ix)*(1d0)  &
                    - kappa*(-1d0)/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix)   &
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(cp(ix,iy)-cpf)/dy(iy)*(1d0) & 
                    + (vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*(cp(ix,iy+1)-cp(ix,iy))/dy(iy)*(1d0) & 
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)/dy(iy)*(1d0) & 
                    + (vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)/dy(iy)*(-1d0) & 
                    - kappa*(-1d0/dy(iy)-1d0/(0.5d0*(dy(iy)+dy(iy+1))))/dy(iy)   
                tmprl(2) =  (vmx(ix-1,iy)-abs(vmx(ix-1,iy)))*0.5d0*cp(ix-1,iy)/dx(ix-1)*(1d0) &
                    - kappa*(1d0)/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix-1)  
                tmprl(3) = (vmy(ix,iy+1)+abs(vmy(ix,iy+1)))*0.5d0*cp(ix,iy+1)/dy(iy+1)*(-1d0)  &
                    - kappa*(1d0)/(0.5d0*(dy(iy)+dy(iy+1)))/dy(iy+1) 
                    
                tmpint(2) = row-1
                tmpint(3) = row+nx
                cnt2 = cnt2 + 2
                
            else 
                bx(row) = (poro(ix,iy)*rho(ix,iy)*cp(ix,iy)+(1d0-poro(ix,iy))*rhom*cpm)*(-tempk(ix,iy))/dt &
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)*(-(temps + 273.15d0))/dy(iy) &
                    - kappa*(temps + 273.15d0)/dy(iy)/dy(iy)  
                tmprl(1) = (poro(ix,iy)*rho(ix,iy)*cp(ix,iy)+(1d0-poro(ix,iy))*rhom*cpm)*(1d0)/dt &
                    + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(cp(ix,iy)-cp(ix-1,iy))/dx(ix)*(1d0)  &
                    + (vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*(cp(ix+1,iy)-cp(ix,iy))/dx(ix)*(1d0)  &
                    ! + (vmx(ix,iy)-vmx(ix-1,iy))*cp(ix,iy)/dx*(1d0) & 
                    + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*cp(ix,iy)/dx(ix)*(1d0)  &
                    + (vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*cp(ix,iy)/dx(ix)*(-1d0)  &
                    - kappa*(-1d0/(0.5d0*(dx(ix)+dx(ix-1)))-1d0/(0.5d0*(dx(ix)+dx(ix+1))))/dx(ix)   &
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(cp(ix,iy)-cpf)/dy(iy)*(1d0)  &
                    + (vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*(cp(ix,iy+1)-cp(ix,iy))/dy(iy)*(1d0)  &
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)/dy(iy)*(1d0)  &
                    + (vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)/dy(iy)*(-1d0)  &
                    - kappa*(-1d0/dy(iy)-1d0/(0.5d0*(dy(iy)+dy(iy+1))))/dy(iy)   
                tmprl(2) = (vmx(ix-1,iy)-abs(vmx(ix-1,iy)))*0.5d0*cp(ix-1,iy)/dx(ix-1)*(1d0)   &
                    - kappa*(1d0)/(0.5d0*(dx(ix)+dx(Ix-1)))/dx(ix-1)  
                tmprl(3) = (vmx(ix+1,iy)+abs(vmx(ix+1,iy)))*0.5d0*cp(ix+1,iy)/dx(ix+1)*(-1d0)   &
                    - kappa*(1d0)/(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix+1)  
                tmprl(4) = (vmy(ix,iy+1)+abs(vmy(ix,iy+1)))*0.5d0*cp(ix,iy+1)/dy(iy+1)*(-1d0)  &
                    - kappa*(1d0)/(0.5d0*(dy(iy)+dy(iy+1)))/dy(iy+1)  
                    
                tmpint(2) = row-1
                tmpint(3) = row+1
                tmpint(4) = row+nx
                cnt2 = cnt2 + 3
                
            endif
        else if (iy==ny) then 
            if (ix==1) then 
                bx(row) = (poro(ix,iy)*rho(ix,iy)*cp(ix,iy)+(1d0-poro(ix,iy))*rhom*cpm)*(-tempk(ix,iy))/dt &
                    ! + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*cp(ix,iy)*(-(tempm + 273.15d0))/dx &
                    ! - kappa*(tempm + 273.15d0)/dx/dx  &
                    - rhom*w*(tempm+273.15d0)*cpm/dx(ix) 
                tmprl(1) = (poro(ix,iy)*rho(ix,iy)*cp(ix,iy)+(1d0-poro(ix,iy))*rhom*cpm)*(1d0)/dt &
                    ! + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(cp(ix,iy)-cpm)/dx*(1d0)  &
                    + (vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*(cp(ix+1,iy)-cp(ix,iy))/dx(ix)*(1d0)  &
                    ! + (vmx(ix,iy)-0d0)*cp(ix,iy)/dx*(1d0)  &
                    ! + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*cp(ix,iy)/dx*(1d0)  &
                    + (vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*cp(ix,iy)/dx(ix)*(-1d0)  &
                    ! - kappa*(-2d0)/dx/dx   &
                    - kappa*(-1d0)/(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix)   &
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(cp(ix,iy)-cp(ix,iy-1))/dy(iy)*(1d0) &
                    ! + (vmy(ix,iy)-vmy(ix,iy-1))*cp(ix,iy)/dy*(1d0)  &
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)/dy(iy)*(1d0) &
                    - kappa*(-1d0)/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)     &
                    - rhom*w*(-1d0)*cpm/dx(ix) 
                tmprl(2) = (vmx(ix+1,iy)+abs(vmx(ix+1,iy)))*0.5d0*cp(ix+1,iy)/dx(ix+1)*(-1d0)   &
                    - kappa*(1d0)/(0.5d0*(dx(ix)+dx(Ix+1)))/dx(ix+1)  
                tmprl(3) =  (vmy(ix,iy-1)-abs(vmy(ix,iy-1)))*0.5d0*cp(ix,iy-1)/dy(iy-1)*(1d0) &
                    - kappa*(1d0)/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy-1) 
                    
                tmpint(2) = row+1
                tmpint(3) = row-nx
                cnt2 = cnt2 + 2
                
            else if (ix==nx) then 
                bx(row) = (poro(ix,iy)*rho(ix,iy)*cp(ix,iy)+(1d0-poro(ix,iy))*rhom*cpm)*(-tempk(ix,iy))/dt 
                tmprl(1) = (poro(ix,iy)*rho(ix,iy)*cp(ix,iy)+(1d0-poro(ix,iy))*rhom*cpm)*(1d0)/dt &
                    + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(cp(ix,iy)-cp(ix-1,iy))/dx(ix)*(1d0)  &
                    ! + (vmx(ix,iy)-vmx(ix-1,iy))*cp(ix,iy)/dx*(1d0)  &
                    + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*cp(ix,iy)/dx(ix)*(1d0)  &
                    - kappa*(-1d0)/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix)   &
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(cp(ix,iy)-cp(ix,iy-1))/dy(iy)*(1d0) & 
                    ! + (vmy(ix,iy)-vmy(ix,iy-1))*cp(ix,iy)/dy*(1d0) & 
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)/dy(iy)*(1d0) & 
                    - kappa*(-1d0)/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)   
                tmprl(2) =  (vmx(ix-1,iy)-abs(vmx(ix-1,iy)))*0.5d0*cp(ix-1,iy)/dx(ix-1)*(1d0) &
                    - kappa*(1d0)/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix-1)  
                tmprl(3) = (vmy(ix,iy-1)-abs(vmy(ix,iy-1)))*0.5d0*cp(ix,iy-1)/dy(iy-1)*(1d0)  &
                    - kappa*(1d0)/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy-1)  
                    
                tmpint(2) = row-1
                tmpint(3) = row-nx
                cnt2 = cnt2 + 2
                
            else 
                bx(row) = (poro(ix,iy)*rho(ix,iy)*cp(ix,iy)+(1d0-poro(ix,iy))*rhom*cpm)*(-tempk(ix,iy))/dt 
                tmprl(1) = (poro(ix,iy)*rho(ix,iy)*cp(ix,iy)+(1d0-poro(ix,iy))*rhom*cpm)*(1d0)/dt &
                    + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(cp(ix,iy)-cp(ix-1,iy))/dx(ix)*(1d0)  &
                    + (vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*(cp(ix+1,iy)-cp(ix,iy))/dx(ix)*(1d0)  &
                    ! + (vmx(ix,iy)-vmx(ix-1,iy))*cp(ix,iy)/dx*(1d0) & 
                    + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*cp(ix,iy)/dx(ix)*(1d0)  &
                    + (vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*cp(ix,iy)/dx(ix)*(-1d0)  &
                    - kappa*(-1d0/(0.5d0*(dx(ix)+dx(ix-1)))-1d0/(0.5d0*(dx(ix)+dx(ix+1))))/dx(ix)   &
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(cp(ix,iy)-cp(ix,iy-1))/dy(iy)*(1d0) & 
                    ! + (vmy(ix,iy)-vmy(ix,iy-1))*cp(ix,iy)/dy*(1d0) & 
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)/dy(iy)*(1d0) & 
                    - kappa*(-1d0)/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)   
                tmprl(2) = (vmx(ix-1,iy)-abs(vmx(ix-1,iy)))*0.5d0*cp(ix-1,iy)/dx(ix-1)*(1d0)   &
                    - kappa*(1d0)/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix-1)  
                tmprl(3) = (vmx(ix+1,iy)+abs(vmx(ix+1,iy)))*0.5d0*cp(ix+1,iy)/dx(ix+1)*(-1d0)   &
                    - kappa*(1d0)/(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix+1)  
                tmprl(4) = (vmy(ix,iy-1)-abs(vmy(ix,iy-1)))*0.5d0*cp(ix,iy-1)/dy(iy-1)*(1d0)  &
                    - kappa*(1d0)/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy-1) 
                    
                tmpint(2) = row-1
                tmpint(3) = row+1
                tmpint(4) = row-nx
                cnt2 = cnt2 + 3
                
            endif
        else 
            if (ix==1) then 
                bx(row) = (poro(ix,iy)*rho(ix,iy)*cp(ix,iy)+(1d0-poro(ix,iy))*rhom*cpm)*(-tempk(ix,iy))/dt &
                    ! + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*cp(ix,iy)*(-(tempm + 273.15d0))/dx &
                    ! - kappa*(tempm + 273.15d0)/dx/dx  
                    - rhom*w*(tempm+273.15d0)*cpm/dx(ix) 
                tmprl(1) = (poro(ix,iy)*rho(ix,iy)*cp(ix,iy)+(1d0-poro(ix,iy))*rhom*cpm)*(1d0)/dt &
                    ! + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(cp(ix,iy)-cpm)/dx*(1d0)  &
                    + (vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*(cp(ix+1,iy)-cp(ix,iy))/dx(ix)*(1d0)  &
                    ! + (vmx(ix,iy)-0d0)*cp(ix,iy)/dx*(1d0) & 
                    ! + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*cp(ix,iy)/dx*(1d0)  &
                    + (vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*cp(ix,iy)/dx(ix)*(-1d0)  &
                    ! - kappa*(-2d0)/dx/dx   &
                    - kappa*(-1d0)/(0.5d0*(dx(ix)+dx(Ix+1)))/dx(ix)   &
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(cp(ix,iy)-cp(ix,iy-1))/dy(iy)*(1d0) & 
                    + (vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*(cp(ix,iy+1)-cp(ix,iy))/dy(iy)*(1d0) & 
                    ! + (vmy(ix,iy)-vmy(ix,iy-1))*cp(ix,iy)/dy*(1d0) & 
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)/dy(iy)*(1d0) & 
                    + (vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)/dy(iy)*(-1d0) & 
                    - kappa*(-1d0/(0.5d0*(dy(iy)+dy(iy-1)))-1d0/(0.5d0*(dy(iy)+dy(iy+1))))/dy(iy)   &
                    - rhom*w*(-1d0)*cpm/dx(ix) 
                tmprl(2) = (vmx(ix+1,iy)+abs(vmx(ix+1,iy)))*0.5d0*cp(ix+1,iy)/dx(ix+1)*(-1d0)   &
                    - kappa*(1d0)/(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix+1)  
                tmprl(3) =  (vmy(ix,iy-1)-abs(vmy(ix,iy-1)))*0.5d0*cp(ix,iy-1)/dy(iy-1)*(1d0)  &
                    - kappa*(1d0)/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy-1) 
                tmprl(4) =  (vmy(ix,iy+1)+abs(vmy(ix,iy+1)))*0.5d0*cp(ix,iy+1)/dy(iy+1)*(1d0)   &
                    - kappa*(1d0)/(0.5d0*(dy(iy)+dy(iy+1)))/dy(iy+1) 
                    
                tmpint(2) = row+1
                tmpint(3) = row-nx
                tmpint(4) = row+nx
                cnt2 = cnt2 + 3
                
            else if (ix==nx) then 
                bx(row) = (poro(ix,iy)*rho(ix,iy)*cp(ix,iy)+(1d0-poro(ix,iy))*rhom*cpm)*(-tempk(ix,iy))/dt 
                tmprl(1) = (poro(ix,iy)*rho(ix,iy)*cp(ix,iy)+(1d0-poro(ix,iy))*rhom*cpm)*(1d0)/dt &
                    + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(cp(ix,iy)-cp(ix-1,iy))/dx(ix)*(1d0)  &
                    ! + (vmx(ix,iy)-vmx(ix-1,iy))*cp(ix,iy)/dx*(1d0) & 
                    + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*cp(ix,iy)/dx(ix)*(1d0)  &
                    - kappa*(-1d0)/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix)   &
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(cp(ix,iy)-cp(ix,iy-1))/dy(iy)*(1d0) & 
                    + (vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*(cp(ix,iy+1)-cp(ix,iy))/dy(iy)*(1d0) & 
                    ! + (vmy(ix,iy)-vmy(ix,iy-1))*cp(ix,iy)/dy*(1d0) & 
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)/dy(iy)*(1d0) & 
                    + (vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)/dy(iy)*(-1d0) & 
                    - kappa*(-1d0/(0.5d0*(dy(iy)+dy(iy-1)))-1d0/(0.5d0*(dy(iy)+dy(iy+1))))/dy(iy)   
                tmprl(2) = (vmx(ix-1,iy)-abs(vmx(ix-1,iy)))*0.5d0*cp(ix-1,iy)/dx(ix-1)*(1d0)  &
                    - kappa*(1d0)/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix-1)  
                tmprl(3) = (vmy(ix,iy-1)-abs(vmy(ix,iy-1)))*0.5d0*cp(ix,iy-1)/dy(iy-1)*(1d0) &
                    - kappa*(1d0)/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy-1) 
                tmprl(4) = (vmy(ix,iy+1)+abs(vmy(ix,iy+1)))*0.5d0*cp(ix,iy+1)/dy(iy+1)*(-1d0) &
                    - kappa*(1d0)/(0.5d0*(dy(iy)+dy(iy+1)))/dy(iy+1) 
                    
                tmpint(2) = row-1
                tmpint(3) = row-nx 
                tmpint(4) = row+nx
                cnt2 = cnt2 + 3
                
            else 
                bx(row) = (poro(ix,iy)*rho(ix,iy)*cp(ix,iy)+(1d0-poro(ix,iy))*rhom*cpm)*(-tempk(ix,iy))/dt 
                tmprl(1) = (poro(ix,iy)*rho(ix,iy)*cp(ix,iy)+(1d0-poro(ix,iy))*rhom*cpm)*(1d0)/dt &
                    + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(cp(ix,iy)-cp(ix-1,iy))/dx(ix)*(1d0)  &
                    + (vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*(cp(ix+1,iy)-cp(ix,iy))/dx(ix)*(1d0)  &
                    ! + (vmx(ix,iy)-vmx(ix-1,iy))*cp(ix,iy)/dx*(1d0) & 
                    + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*cp(ix,iy)/dx(ix)*(1d0)  &
                    + (vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*cp(ix,iy)/dx(ix)*(-1d0)  &
                    - kappa*(-1d0/(0.5d0*(dx(ix)+dx(Ix-1)))-1d0/(0.5d0*(dx(ix)+dx(Ix+1))))/dx(ix)   &
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(cp(ix,iy)-cp(ix,iy-1))/dy(iy)*(1d0) & 
                    + (vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*(cp(ix,iy+1)-cp(ix,iy))/dy(iy)*(1d0) & 
                    ! + (vmy(ix,iy)-vmy(ix,iy-1))*cp(ix,iy)/dy*(1d0) & 
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)/dy(iy)*(1d0) & 
                    + (vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)/dy(iy)*(-1d0) & 
                    - kappa*(-1d0/(0.5d0*(dy(iy)+dy(iy-1)))-1d0/(0.5d0*(dy(iy)+dy(iy+1))))/dy(iy)   
                tmprl(3) = (vmx(ix+1,iy)+abs(vmx(ix+1,iy)))*0.5d0*cp(ix+1,iy)/dx(ix+1)*(-1d0)   &
                    - kappa*(1d0)/(0.5d0*(dx(ix)+dx(Ix+1)))/dx(ix+1)  
                tmprl(2) = (vmx(ix-1,iy)-abs(vmx(ix-1,iy)))*0.5d0*cp(ix-1,iy)/dx(ix-1)*(1d0)   &
                    - kappa*(1d0)/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix-1)  
                tmprl(4) = (vmy(ix,iy-1)-abs(vmy(ix,iy-1)))*0.5d0*cp(ix,iy-1)/dy(iy-1)*(1d0) &
                    - kappa*(1d0)/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy-1) 
                tmprl(5) = (vmy(ix,iy+1)+abs(vmy(ix,iy+1)))*0.5d0*cp(ix,iy+1)/dy(iy+1)*(-1d0) &
                    - kappa*(1d0)/(0.5d0*(dy(iy)+dy(iy+1)))/dy(iy+1) 
                    
                tmpint(2) = row-1
                tmpint(3) = row+1
                tmpint(4) = row-nx
                tmpint(5) = row+nx
                cnt2 = cnt2 + 4
                
            endif
        endif
        
        ! print*,cnt2
        if (cnt2 < 1 .or. cnt2 > 5) then 
            print *,"FATAL ERROR"
            stop
        endif
        tmpint2 = 0 
        ! print*,'here'
        call heapsort2(cnt2, tmpint(1:cnt2),tmpint2(1:cnt2))
        do i = 1,cnt2 
            ai(cnt+i)=tmpint(i) - 1
            ax(cnt+i)=tmprl(tmpint2(i))
        enddo
        cnt = cnt + cnt2 
        ap(row+1) = cnt
        
    enddo
enddo

bx = - bx

! call DGESV(nmx,int(1),amx,nmx,IPIV,ymx,nmx,infobls) 


! solving matrix with UMFPACK (following is pasted from umfpack_simple.f90)
! Set the default control parameters.
call umf4def( control )
! From the matrix data, create the symbolic factorization information.
call umf4sym ( n, n, ap, ai, ax, symbolic, control, info )
if ( info(1) < 0.0D+00 ) then
    write ( *, * ) ''
    write ( *, *) 'UMFPACK_SIMPLE - Fatal error!'
    write ( *, * ) '  UMF4SYM returns INFO(1) = ', info(1)
    stop 1
end if
! From the symbolic factorization information, carry out the numeric factorization.
call umf4num ( ap, ai, ax, symbolic, numeric, control, info )
if ( info(1) < 0.0D+00 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'UMFPACK_SIMPLE - Fatal error!'
    write ( *, '(a,g14.6)' ) '  UMF4NUM returns INFO(1) = ', info(1)
    stop 1
end if
!  Free the memory associated with the symbolic factorization.
call umf4fsym ( symbolic )
! Solve the linear system.
sys = 0
call umf4sol ( sys, kai, bx, numeric, control, info )
if ( info(1) < 0.0D+00 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'UMFPACK_SIMPLE - Fatal error!'
    write ( *, '(a,g14.6)' ) '  UMF4SOL returns INFO(1) = ', info(1)
    stop 1
end if
! Free the memory associated with the numeric factorization.
call umf4fnum ( numeric )
!  Print the solution.
write ( *, * ) ''
write ( *, * ) '  Computed solution'
write ( *, * ) ''

tempk_pre = tempk

do iy=1,ny
    do ix = 1,nx
        row = (iy-1)*nx + ix
        tempk(ix,iy) = kai(row)
    enddo
enddo

temp = tempk - 273.15d0

it = it + 1

flxt = 0d0
flxadv = 0d0 
flxcnd = 0d0 
flxres =0d0
flxb = 0d0

do iy = 1,ny
    do ix=1,nx
        
        if (iy == 1) then 
            if (ix==1) then 
                flxt(ix,iy) = flxt(ix,iy) &
                    + (poro(ix,iy)*rho(ix,iy)*cp(ix,iy)+(1d0-poro(ix,iy))*rhom*cpm)*(tempk(ix,iy)-tempk_pre(ix,iy))/dt
                    
                flxadv(ix,iy) = flxadv(ix,iy) & 
                    ! + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(cp(ix,iy)-cpm)*tempk(ix,iy)/dx  &
                    + (vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*(cp(ix+1,iy)-cp(ix,iy))*tempk(ix,iy)/dx(ix)  &
                    ! + (vmx(ix,iy) - 0d0)*cp(ix,iy)*tempk(ix,iy)/dx  &
                    + (vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*cp(ix,iy)*(tempk(ix+1,iy)-tempk(ix,iy))/dx(ix)  &
                    ! 
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(cp(ix,iy)-cpf)*tempk(ix,iy)/dy(iy) & 
                    + (vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*(cp(ix,iy+1)-cp(ix,iy))*tempk(ix,iy)/dy(iy) & 
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)*(tempk(ix,iy)-(temps+273.15d0))/dy(iy)  &
                    + (vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)*(tempk(ix,iy+1)-tempk(ix,iy))/dy(iy)  
                    
                flxcnd(ix,iy) = flxcnd(ix,iy) & 
                    - kappa*(tempk(ix+1,iy)-tempk(ix,iy))/(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix)  &
                    - kappa*((tempk(ix,iy+1)-tempk(ix,iy))/(0.5d0*(dy(iy)+dy(iy+1)))  &
                        -(tempk(ix,iy)-(temps+273.15d0))/dy(iy))/dy(iy)
                
                flxb(ix,iy) = flxb(ix,iy) &
                    - rhom*w*((tempm+273.15d0)-tempk(ix,iy))*cpm/dx(ix) 
                
            else if (ix==nx) then 
                flxt(ix,iy) = flxt(ix,iy) &
                    + (poro(ix,iy)*rho(ix,iy)*cp(ix,iy)+(1d0-poro(ix,iy))*rhom*cpm)*(tempk(ix,iy)-tempk_pre(ix,iy))/dt
                    
                flxadv(ix,iy) = flxadv(ix,iy) & 
                    + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(cp(ix,iy)-cp(ix-1,iy))*tempk(ix,iy)/dx(ix)  &
                    ! + (vmx(ix,iy) - vmx(ix-1,iy))*cp(ix,iy)*tempk(ix,iy)/dx  &
                    + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*cp(ix,iy)*(tempk(ix,iy)-tempk(ix-1,iy))/dx(ix)  &
                    !
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(cp(ix,iy)-cpf)*tempk(ix,iy)/dy(iy) & 
                    + (vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*(cp(ix,iy+1)-cp(ix,iy))*tempk(ix,iy)/dy(iy) & 
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)*(tempk(ix,iy)-(temps+273.15d0))/dy(iy)  &
                    + (vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)*(tempk(ix,iy+1)-tempk(ix,iy))/dy(iy)  
                    
                flxcnd(ix,iy) = flxcnd(ix,iy) & 
                    - kappa*(tempk(ix-1,iy)-tempk(ix,iy))/(0.5d0*(dx(ix)+dx(Ix-1)))/dx(ix)  &
                    - kappa*((tempk(ix,iy+1)-tempk(ix,iy))/(0.5d0*(dy(iy)+dy(iy+1)))  &
                        -(tempk(ix,iy)-(temps+273.15d0))/dy(iy))/dy(iy)
                
            else 
                flxt(ix,iy) = flxt(ix,iy) &
                    + (poro(ix,iy)*rho(ix,iy)*cp(ix,iy)+(1d0-poro(ix,iy))*rhom*cpm)*(tempk(ix,iy)-tempk_pre(ix,iy))/dt
                    
                flxadv(ix,iy) = flxadv(ix,iy) & 
                    + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(cp(ix,iy)-cp(ix-1,iy))*tempk(ix,iy)/dx(ix)  &
                    + (vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*(cp(ix+1,iy)-cp(ix,iy))*tempk(ix,iy)/dx(ix)  &
                    ! + (vmx(ix,iy) - vmx(ix-1,iy))*cp(ix,iy)*tempk(ix,iy)/dx  &
                    + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*cp(ix,iy)*(tempk(ix,iy)-tempk(ix-1,iy))/dx(ix)  &
                    + (vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*cp(ix,iy)*(tempk(ix+1,iy)-tempk(ix,iy))/dx(ix)  &
                    !
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(cp(ix,iy)-cpf)*tempk(ix,iy)/dy(iy) & 
                    + (vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*(cp(ix,iy+1)-cp(ix,iy))*tempk(ix,iy)/dy(iy) & 
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)*(tempk(ix,iy)-(temps+273.15d0))/dy(iy)  &
                    + (vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)*(tempk(ix,iy+1)-tempk(ix,iy))/dy(iy)  
                    
                flxcnd(ix,iy) = flxcnd(ix,iy) & 
                    - kappa*((tempk(ix+1,iy)-tempk(ix,iy))/(0.5d0*(dx(ix)+dx(Ix+1)))  &
                        -(tempk(ix,iy)-tempk(ix-1,iy))/(0.5d0*(dx(ix)+dx(Ix-1))))/dx(ix)  &
                    - kappa*((tempk(ix,iy+1)-tempk(ix,iy))/(0.5d0*(dy(iy)+dy(iy+1)))  &
                        -(tempk(ix,iy)-(temps+273.15d0))/dy(iy))/dy(iy)
                
            endif
        else if (iy==ny) then 
            if (ix==1) then 
                flxt(ix,iy) = flxt(ix,iy) &
                    + (poro(ix,iy)*rho(ix,iy)*cp(ix,iy)+(1d0-poro(ix,iy))*rhom*cpm)*(tempk(ix,iy)-tempk_pre(ix,iy))/dt
                    
                flxadv(ix,iy) = flxadv(ix,iy) & 
                    ! + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(cp(ix,iy)-cpm)*tempk(ix,iy)/dx  &
                    + (vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*(cp(ix+1,iy)-cp(ix,iy))*tempk(ix,iy)/dx(ix)  &
                    ! + (vmx(ix,iy) - 0d0)*cp(ix,iy)*tempk(ix,iy)/dx  &
                    + (vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*cp(ix,iy)*(tempk(ix+1,iy)-tempk(ix,iy))/dx(ix)  &
                    !
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(cp(ix,iy)-cp(ix,iy-1))*tempk(ix,iy)/dy(iy) & 
                    ! + (vmy(ix,iy)-vmy(ix,iy-1))*cp(ix,iy)*tempk(ix,iy)/dy  &
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)*(tempk(ix,iy)-tempk(ix,iy-1))/dy(iy)  
                    
                flxcnd(ix,iy) = flxcnd(ix,iy) & 
                    - kappa*(tempk(ix+1,iy)-1d0*tempk(ix,iy))/(0.5d0*(dx(ix)+dx(Ix+1)))/dx(ix)  &
                    - kappa*(tempk(ix,iy-1)-1d0*tempk(ix,iy))/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)
                    
                flxb(ix,iy) = flxb(ix,iy) & 
                    - rhom*w*((tempm+273.15d0)-tempk(ix,iy))*cpm/dx(ix)
                
            else if (ix==nx) then 
                flxt(ix,iy) = flxt(ix,iy) &
                    + (poro(ix,iy)*rho(ix,iy)*cp(ix,iy)+(1d0-poro(ix,iy))*rhom*cpm)*(tempk(ix,iy)-tempk_pre(ix,iy))/dt
                    
                flxadv(ix,iy) = flxadv(ix,iy) & 
                    + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(cp(ix,iy)-cp(ix-1,iy))*tempk(ix,iy)/dx(ix)  &
                    ! + (vmx(ix,iy) - vmx(ix-1,iy))*cp(ix,iy)*tempk(ix,iy)/dx  &
                    + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*cp(ix,iy)*(tempk(ix,iy)-tempk(ix-1,iy))/dx(ix)  &
                    !
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(cp(ix,iy)-cp(ix,iy-1))*tempk(ix,iy)/dy(iy) & 
                    ! + (vmy(ix,iy)-vmy(ix,iy-1))*cp(ix,iy)*tempk(ix,iy)/dy  &
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)*(tempk(ix,iy)-tempk(ix,iy-1))/dy(iy)  
                    
                flxcnd(ix,iy) = flxcnd(ix,iy) & 
                    - kappa*(tempk(ix-1,iy)-1d0*tempk(ix,iy))/(0.5d0*(dx(ix)+dx(Ix-1)))/dx(ix)  &
                    - kappa*(tempk(ix,iy-1)-1d0*tempk(ix,iy))/(0.5d0*(dy(iy)+dy(Iy-1)))/dy(iy)
                
            else 
                flxt(ix,iy) = flxt(ix,iy) &
                    + (poro(ix,iy)*rho(ix,iy)*cp(ix,iy)+(1d0-poro(ix,iy))*rhom*cpm)*(tempk(ix,iy)-tempk_pre(ix,iy))/dt
                    
                flxadv(ix,iy) = flxadv(ix,iy) & 
                    + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(cp(ix,iy)-cp(ix-1,iy))*tempk(ix,iy)/dx(ix)  &
                    + (vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*(cp(ix+1,iy)-cp(ix,iy))*tempk(ix,iy)/dx(ix)  &
                    ! + (vmx(ix,iy) - vmx(ix-1,iy))*cp(ix,iy)*tempk(ix,iy)/dx  &
                    + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*cp(ix,iy)*(tempk(ix,iy)-tempk(ix-1,iy))/dx(ix)  &
                    + (vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*cp(ix,iy)*(tempk(ix+1,iy)-tempk(ix,iy))/dx(ix)  &
                    !
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(cp(ix,iy)-cp(ix,iy-1))*tempk(ix,iy)/dy(iy) & 
                    ! + (vmy(ix,iy)-vmy(ix,iy-1))*cp(ix,iy)*tempk(ix,iy)/dy  &
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)*(tempk(ix,iy)-tempk(ix,iy-1))/dy(iy)  
                    
                flxcnd(ix,iy) = flxcnd(ix,iy) & 
                    - kappa*((tempk(ix+1,iy)-tempk(ix,iy))/(0.5d0*(dx(ix)+dx(Ix+1))) &
                        -(tempk(ix,iy)-tempk(ix-1,iy))/(0.5d0*(dx(ix)+dx(Ix-1))))/dx(ix)  &
                    - kappa*(tempk(ix,iy-1)-1d0*tempk(ix,iy))/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)
                
            endif
        else 
            if (ix==1) then 
                flxt(ix,iy) = flxt(ix,iy) &
                    + (poro(ix,iy)*rho(ix,iy)*cp(ix,iy)+(1d0-poro(ix,iy))*rhom*cpm)*(tempk(ix,iy)-tempk_pre(ix,iy))/dt
                    
                flxadv(ix,iy) = flxadv(ix,iy) & 
                    ! + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(cp(ix,iy)-cpm)*tempk(ix,iy)/dx  &
                    + (vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*(cp(ix+1,iy)-cp(ix,iy))*tempk(ix,iy)/dx(ix)  &
                    ! + (vmx(ix,iy) - 0d0)*cp(ix,iy)*tempk(ix,iy)/dx  &
                    + (vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*cp(ix,iy)*(tempk(ix+1,iy)-tempk(ix,iy))/dx(ix)  &
                    !
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(cp(ix,iy)-cp(ix,iy-1))*tempk(ix,iy)/dy(iy) & 
                    + (vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*(cp(ix,iy+1)-cp(ix,iy))*tempk(ix,iy)/dy(iy) & 
                    ! + (vmy(ix,iy)-vmy(ix,iy-1))*cp(ix,iy)*tempk(ix,iy)/dy  &
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)*(tempk(ix,iy)-tempk(ix,iy-1))/dy(iy)  &
                    + (vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)*(tempk(ix,iy+1)-tempk(ix,iy))/dy(iy)  
                    
                flxcnd(ix,iy) = flxcnd(ix,iy) & 
                    - kappa*(tempk(ix+1,iy)-1d0*tempk(ix,iy))/(0.5d0*(dx(ix)+dx(Ix+1)))/dx(ix)  &
                    - kappa*((tempk(ix,iy+1)-tempk(ix,iy))/(0.5d0*(dy(iy)+dy(iy+1)))  &
                        -(tempk(ix,iy)-tempk(ix,iy-1))/(0.5d0*(dy(iy)+dy(iy-1))))/dy(iy)
                    
                flxb(ix,iy) = flxb(ix,iy) & 
                    - rhom*w*((tempm+273.15d0)-tempk(ix,iy))*cpm/dx(ix)
                
            else if (ix==nx) then 
                flxt(ix,iy) = flxt(ix,iy) &
                    + (poro(ix,iy)*rho(ix,iy)*cp(ix,iy)+(1d0-poro(ix,iy))*rhom*cpm)*(tempk(ix,iy)-tempk_pre(ix,iy))/dt
                    
                flxadv(ix,iy) = flxadv(ix,iy) & 
                    + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(cp(ix,iy)-cp(ix-1,iy))*tempk(ix,iy)/dx(ix)  &
                    ! + (vmx(ix,iy) - vmx(ix-1,iy))*cp(ix,iy)*tempk(ix,iy)/dx  &
                    + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*cp(ix,iy)*(tempk(ix,iy)-tempk(ix-1,iy))/dx(ix)  &
                    !
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(cp(ix,iy)-cp(ix,iy-1))*tempk(ix,iy)/dy(iy) & 
                    + (vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*(cp(ix,iy+1)-cp(ix,iy))*tempk(ix,iy)/dy(iy) & 
                    ! + (vmy(ix,iy)-vmy(ix,iy-1))*cp(ix,iy)*tempk(ix,iy)/dy  &
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)*(tempk(ix,iy)-tempk(ix,iy-1))/dy(iy)  &
                    + (vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)*(tempk(ix,iy+1)-tempk(ix,iy))/dy(iy)  
                    
                flxcnd(ix,iy) = flxcnd(ix,iy) & 
                    - kappa*(tempk(ix-1,iy)-1d0*tempk(ix,iy))/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix)  &
                    - kappa*((tempk(ix,iy+1)-tempk(ix,iy))/(0.5d0*(dy(iy)+dy(iy+1)))  &
                        -(tempk(ix,iy)-tempk(ix,iy-1))/(0.5d0*(dy(iy)+dy(iy-1))))/dy(iy)
                
            else 
                flxt(ix,iy) = flxt(ix,iy) &
                    + (poro(ix,iy)*rho(ix,iy)*cp(ix,iy)+(1d0-poro(ix,iy))*rhom*cpm)*(tempk(ix,iy)-tempk_pre(ix,iy))/dt
                    
                flxadv(ix,iy) = flxadv(ix,iy) & 
                    + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(cp(ix,iy)-cp(ix-1,iy))*tempk(ix,iy)/dx(ix)  &
                    + (vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*(cp(ix+1,iy)-cp(ix,iy))*tempk(ix,iy)/dx(ix)  &
                    ! + (vmx(ix,iy) - vmx(ix-1,iy))*cp(ix,iy)*tempk(ix,iy)/dx  &
                    + (vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*cp(ix,iy)*(tempk(ix,iy)-tempk(ix-1,iy))/dx(ix)  &
                    + (vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*cp(ix,iy)*(tempk(ix+1,iy)-tempk(ix,iy))/dx(ix)  &
                    !
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(cp(ix,iy)-cp(ix,iy-1))*tempk(ix,iy)/dy(iy) & 
                    + (vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*(cp(ix,iy+1)-cp(ix,iy))*tempk(ix,iy)/dy(iy) & 
                    ! + (vmy(ix,iy)-vmy(ix,iy-1))*cp(ix,iy)*tempk(ix,iy)/dy  &
                    + (vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)*(tempk(ix,iy)-tempk(ix,iy-1))/dy(iy)  &
                    + (vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*cp(ix,iy)*(tempk(ix,iy+1)-tempk(ix,iy))/dy(iy)  
                    
                flxcnd(ix,iy) = flxcnd(ix,iy) & 
                    - kappa*((tempk(ix+1,iy)-tempk(ix,iy))/(0.5d0*(dx(Ix)+dx(Ix+1)))  &
                        -(tempk(ix,iy)-tempk(ix-1,iy))/(0.5d0*(dx(ix)+dx(ix-1))))/dx(ix)  &
                    - kappa*((tempk(ix,iy+1)-tempk(ix,iy))/(0.5d0*(dy(iy)+dy(iy+1)))  &
                        -(tempk(ix,iy)-tempk(ix,iy-1))/(0.5d0*(dy(iy)+dy(iy-1))))/dy(iy)
                
            endif
        endif
        
    enddo
enddo

! flxt = flxt*dx*dy*rl ! [J sec-1 = W]
! flxadv = flxadv*dx*dy*rl 
! flxcnd = flxcnd*dx*dy*rl
! flxb = flxb*dx*dy*rl
flxt = flxt*ds*rl ! [J sec-1 = W]
flxadv = flxadv*ds*rl 
flxcnd = flxcnd*ds*rl
flxb = flxb*ds*rl

flxres = flxt+flxadv+flxcnd+flxb

print*,'flxt, flxadv, flxcnd, flxb, res'
print'(5E9.2)',sum(flxt), sum(flxadv), sum(flxcnd), sum(flxb), sum(flxres)

if (flgesc) then 
    if (abs(sum(flxt)/sum(flxb)) < 1d-5) then 
        open(unit=100,file=trim(adjustl(workdir))//'log.txt',action='write',status='unknown')
        write(100,*)"!!! Enough iteration? Tentative end !!!",it,sum(flxt), sum(flxadv), sum(flxcnd), sum(flxb)
        close(100)
        exit 
    endif
endif

enddo 

open(unit=100,file=trim(adjustl(workdir))//'rho.txt',action='write',status='unknown')
open(unit=200,file=trim(adjustl(workdir))//'visc.txt',action='write',status='unknown')
open(unit=300,file=trim(adjustl(workdir))//'cp.txt',action='write',status='unknown')
open(unit=400,file=trim(adjustl(workdir))//'perm.txt',action='write',status='unknown')
open(unit=500,file=trim(adjustl(workdir))//'temp.txt',action='write',status='unknown')
open(unit=600,file=trim(adjustl(workdir))//'vmx.txt',action='write',status='unknown')
open(unit=700,file=trim(adjustl(workdir))//'vmy.txt',action='write',status='unknown')
open(unit=800,file=trim(adjustl(workdir))//'Dv.txt',action='write',status='unknown')
open(unit=900,file=trim(adjustl(workdir))//'qy.txt',action='write',status='unknown')
open(unit=990,file=trim(adjustl(workdir))//'qx.txt',action='write',status='unknown')
open(unit=980,file=trim(adjustl(workdir))//'pres.txt',action='write',status='unknown')
open(unit=970,file=trim(adjustl(workdir))//'cd.txt',action='write',status='unknown')
open(unit=960,file=trim(adjustl(workdir))//'cd2.txt',action='write',status='unknown')
do iy=1,ny
    write(100,*) (rho(ix,iy),ix=1,nx) 
    write(200,*) (visc(ix,iy),ix=1,nx) 
    write(300,*) (cp(ix,iy),ix=1,nx) 
    write(400,*) (perm(ix,iy),ix=1,nx) 
    write(500,*) (temp(ix,iy),ix=1,nx) 
    write(600,*) (vmx(ix,iy),ix=1,nx) 
    write(700,*) (vmy(ix,iy),ix=1,nx) 
    write(800,*) (dv(ix,iy),ix=1,nx) 
    write(900,*) (qvy(ix,iy),ix=1,nx)
    write(990,*) (qvx(ix,iy),ix=1,nx) 
    write(980,*) (pres(ix,iy)/1d6,ix=1,nx) 
    write(970,*) (cd(ix,iy),ix=1,nx) 
    write(960,*) (cd2(ix,iy),ix=1,nx) 
enddo
close(100)
close(200)
close(300)
close(400)
close(500)
close(600)
close(700)
close(800)
close(900)
close(990)
close(980)
close(970)
close(960)


open(unit=100,file=trim(adjustl(workdir))//'heatflx.txt',action='write',status='unknown')
write(100,*)'flxt, flxadv, flxcnd, flxb, flxres'
write(100,*) sum(flxt), sum(flxadv), sum(flxcnd), sum(flxb), sum(flxres)
close(100)


open(unit=100,file=trim(adjustl(workdir))//'wflx.txt',action='write',status='unknown')
write(100,*)'total exchange (marker), total exchange (cell), net (m), net (c) [kg/yr]'
write(100,*) 0.5d0*sum(abs(vy(:,1))*dx(:)*rl*yr2sec), 0.5d0*sum(abs(vmy(:,1))*dx(:)*rl*yr2sec) &
    , sum(vy(:,1)*dx(:)*rl*yr2sec), sum(vmy(:,1)*dx(:)*rl*yr2sec)
close(100)

End program
!**************************************************************************************************************************************


!**************************************************************************************************************************************
subroutine heapsort2(n,array,turn)
!!!  from http://slpr.sakura.ne.jp/qp/sortf90/
  implicit none
  integer,intent(in)::n  ! size of array 
  integer,intent(out)::turn(1:n)  ! order of input array to make sorted array
  integer,intent(inout)::array(1:n)  ! sorted array in increasing order
 
  integer::i,k,j,l,m
  integer::t
 
  if(n.le.0)then
     write(6,*)"Error, at heapsort"; stop
  endif
  if(n.eq.1)return

  do i=1,N
     turn(i)=i
  enddo

  l=n/2+1
  k=n
  do while(k.ne.1)
     if(l.gt.1)then
        l=l-1
        t=array(l)
        m=turn(l)
     else
        t=array(k)
        m=turn(k)
        array(k)=array(1)
        turn(k)=turn(1)
        k=k-1
        if(k.eq.1) then
           array(1)=t
           turn(1)=m
           exit
        endif
     endif
     i=l
     j=l+l
     do while(j.le.k)
        if(j.lt.k)then
           if(array(j).lt.array(j+1))j=j+1
        endif
        if (t.lt.array(j))then
           array(i)=array(j)
           turn(i)=turn(j)
           i=j
           j=j+j
        else
           j=k+1
        endif
     enddo
     array(i)=t
     turn(i)=m
  enddo

  return
end subroutine heapsort2
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine make_grid(  &
    nz,beta,ztot  &! input 
    ,dz         &! output
    )

implicit none 

integer(kind=4),intent(in)::nz
real(kind=8),intent(in)::beta,ztot
real(kind=8),intent(out)::dz(nz)

! local variables 
integer(kind=4)::iz
real(kind=8)::eta(nz),alpha,z(nz)


do iz=1,nz
    eta(iz) = 1d0*iz/(1d0*nz)
enddo

z = ztot*((beta+1d0)-(beta-1d0)*(((beta+1d0)/(beta-1d0))**(1d0-eta)))  &
    /((((beta+1d0)/(beta-1d0))**(1d0-eta))+1d0)

dz = z/sum(z)*ztot

! do iz=1,nz  ! depth is defined at the middle of individual layers 
    ! if (iz==1) z(iz)=dz(iz)*0.5d0  
    ! if (iz/=1) z(iz) = z(iz-1)+dz(iz-1)*0.5d0 + 0.5d0*dz(iz)
! enddo

! open(unit=250,file='test.txt',action = 'write',status='replace')
! do iz=1,nz
    ! write(250,*) eta(iz),dz(iz),z(iz)
! enddo 
! close(250)
! stop
        
endsubroutine make_grid
!**************************************************************************************************************************************