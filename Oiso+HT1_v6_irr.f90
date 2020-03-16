
!**************************************************************************************************************************************
program Oiso_HT1
! oxygen isotope calculation using 2D data for flow, temperature and density
! looking working 
! modified from v3 to include the change in rho for diffusion term
! modified from v4 to include the temp dependence of diffusion coefficient
! remove non-sparse matrix solver to save the memory
! from v5 to remove porosity before q
! from v6 trying irregular grid 
implicit none

integer(kind=4),parameter :: nx = 200, ny = 320
real(kind=8),parameter :: xmax = 30d3   !  m
real(kind=8),parameter :: ymax = 12d3    !  m 
real(kind=8),parameter :: yr2sec = 60d0*60d0*24d0*365.25d0 ! sec/yr
real(kind=8) :: rl = 1d8 ! m ridge length
real(kind=8) dx(nx), dy(ny), dt, time  ! m,m, yr, yr
real(kind=8) xe(nx+1),ye(ny+1),xm(nx),ym(ny),ds(nx,ny)
real(kind=8) temp(nx,ny), rho(nx,ny)  ! need to be read 
real(kind=8) vmx(nx,ny), vmy(nx,ny) ! need to be read 
real(kind=8) visc(nx,ny), vabs(nx,ny), wr(nx,ny), omega(nx,ny)
real(kind=8) poro(nx,ny), rhob(nx,ny), rxn(nx,ny) 
real(kind=8) drxn_dfr(nx,ny),drxn_dfp(nx,ny) 
real(kind=8) :: rhom = 3d3 ! rock density [kg m-3]
real(kind=8) :: w = 3d-2 ! m/yr spreading rate 
! real(kind=8) :: w = 9d-2 ! m/yr spreading rate 
! real(kind=8) :: w = 1d-2 ! m/yr spreading rate 
! real(kind=8) :: w = 3d-1 ! m/yr spreading rate 
real(kind=8) :: temps = 2d0 ! C seawater temp
real(kind=8) :: poroi = 0.05d0 !  
real(kind=8) :: disp = 10d0 !  dispersivity [m]
integer(kind = 4) ix,iy, isw
! o isotope 
real(kind=8) fr(nx,ny), fp(nx,ny), kex(nx,ny), alfa(nx,ny), dif(nx,ny)
real(kind=8) frx(nx,ny), fpx(nx,ny), o18r(nx,ny),o18p(nx,ny)
real(kind=8) fsw, fri
real(kind=8) :: difi = 1d-9*yr2sec  ! m2 yr-1
real(kind=8) :: o18swi = -0d0  ! o/oo
real(kind=8), parameter :: rsmow = 0.0020052d0 
real(kind=8), parameter :: rg = 8.3144d-3 !  kJ mol^-1 K^-1
real(kind=8), parameter :: E = 50d0 ! kJ
real(kind=8), parameter :: kref = 10d0**(-8.5d0) ! mol-1 kg yr-1
! real(kind=8), parameter :: kref = 10d0**(-9.5d0) ! mol-1 kg yr-1
real(kind=8), parameter :: tempref = 5d0 ! C
real(kind=8), parameter :: o18ri = 5.7d0 ! o/oo of solid rock
real(kind=8), parameter :: ms = 0.5d0/16.0d0*1.0d3 ! mol kg-1
real(kind=8), parameter :: mw = 16d0/18d0/16.0d0*1.0d3  !  mol kg-1
real(kind=8), parameter :: tempk_0 = 273.15d0  !  K
real(kind=8) :: beta = 0.876002291d0 ! anddesite (Zhao and Zheng, 2003)
real(kind=8) rflxadv(nx,ny), rflxrxn(nx,ny), rflxt(nx,ny), pflxadv(nx,ny), pflxrxn(nx,ny), pflxdif(nx,ny), pflxt(nx,ny) 
real(kind=8) rflxrxnpls(nx,ny), rflxrxnmns(nx,ny)
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
integer(kind=4) row,col, it, itr
real(kind=8) error
real(kind=8),parameter :: tol = 1d-12
integer(kind=4) cnt, cnt2, ixp, ixn, iyp, iyn, tmpint(6),tmpint2(6)
real(kind=8) tmprl(6)  
character*100 workdir
character*2 intsw
! functions
real(kind=8) :: f2d, r2d, d2f, d2r
real(kind=8) beta_grid 

workdir = 'C:/Users/YK/Desktop/HT_res/'//  &
    'perm_expexp_-16_8-11_3_zh500_spx1_200x320_irr-20191030'  &
    //'/'
 
!  initial conditions 

poro = poroi

xe(1) = 0d0
ye(1) = 0d0
beta_grid = 1.5d0
call make_grid(  &
    nx,beta_grid,xmax  &! input 
    ,dx         &! output
    )
call make_grid(  &
    ny,beta_grid,ymax  &! input 
    ,dy         &! output
    )
! dx = xmax/nx
! dy = ymax/ny
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

open(unit=100,file=trim(adjustl(workdir))//'x_edges_iso.txt',action='write',status='unknown')
do ix=1,nx+1
    write(100,*) xe(ix)
enddo 
close(100)
open(unit=100,file=trim(adjustl(workdir))//'y_edges_iso.txt',action='write',status='unknown')
do iy=1,ny+1
    write(100,*) ye(iy)
enddo 
close(100)

! print*, ym

open(unit=100,file=trim(adjustl(workdir))//'rho.txt',action='read',status='old')
open(unit=500,file=trim(adjustl(workdir))//'temp.txt',action='read',status='old')
open(unit=600,file=trim(adjustl(workdir))//'vmx.txt',action='read',status='old')
open(unit=700,file=trim(adjustl(workdir))//'vmy.txt',action='read',status='old')
open(unit=800,file=trim(adjustl(workdir))//'visc.txt',action='read',status='old')
do iy=1,ny
    read(100,*) (rho(ix,iy),ix=1,nx)  
    read(500,*) (temp(ix,iy),ix=1,nx) 
    read(600,*) (vmx(ix,iy),ix=1,nx) 
    read(700,*) (vmy(ix,iy),ix=1,nx) 
    read(800,*) (visc(ix,iy),ix=1,nx) ! dynamic viscosity [Pa s]
enddo
close(100)
close(500)
close(600)
close(700)
close(800)

open(unit=100,file=trim(adjustl(workdir))//'temp_chk.txt',action='write',status='replace')
do iy=1,ny
    write(100,*) (temp(ix,iy),ix=1,nx) 
enddo
close(100)

vmx = vmx/rho*yr2sec/poro
vmy = vmy/rho*yr2sec/poro

vabs = sqrt(vmx*vmx+vmy*vmy)

! temp = 2d0
! vmx = 0d0

dif = difi
! modified Stokes-Einstein (Krynicki et al., 1978)
dif = 6.9d-15*(tempk_0+temp)/visc  ! water diffusion [m2 s-1]
dif = dif/sqrt(20d0/18d0)  ! correction for H218O (Harris et al., 1979)

rhob = (1d0-poro)*rhom + poro*rho


wr = poro*rho*mw*vabs/((1d0-poro)*rhom*ms*w)

alfa = (6.673d6/(tempk_0+temp)**2.0d0+10.398d3/(tempk_0+temp)-4.78d0)*exp((1.0d0-beta)/(rg*(tempk_0+temp)))*beta  &
    -(2.194d6/(tempk_0+temp)**2.0d0+15.163d3/(tempk_0+temp)-4.72d0)+1.767d0*(2.0d0*beta-1.0d0)

open(unit=100,file=trim(adjustl(workdir))//'alfa.txt',action='write',status='replace')
open(unit=200,file=trim(adjustl(workdir))//'dif.txt',action='write',status='replace')
open(unit=300,file=trim(adjustl(workdir))//'vabs.txt',action='write',status='replace')
open(unit=400,file=trim(adjustl(workdir))//'wr.txt',action='write',status='replace')
do iy=1,ny
    write(100,*) (alfa(ix,iy),ix=1,nx) 
    write(200,*) (dif(ix,iy),ix=1,nx) 
    write(300,*) (vabs(ix,iy),ix=1,nx) 
    write(400,*) (wr(ix,iy),ix=1,nx) 
enddo
close(100)
close(200)
close(300)
close(400)

alfa = exp(alfa/1d3) ! changing from o/oo to fractionation factor

dif = dif*yr2sec   ! converting sec to yr 

dif = dif*poro**1.4d0 + disp*vabs  !  assuming homogeneous dispersivity 

kex = kref*exp(-E*(1.0d0/(tempk_0+temp)-1.0d0/(tempk_0+tempref))/rg)

open(unit=100,file=trim(adjustl(workdir))//'k.txt',action='write',status='replace')
do iy=1,ny
    write(100,*) (log10(kex(ix,iy)),ix=1,nx) 
enddo
close(100)

do isw = 1,13, 2
! do isw = 1,1

o18swi = -1d0*(isw-1)
write(intsw,'(I2.2)') int(abs(o18swi))

fsw = d2f(o18swi)
fri = d2f(o18ri)

print*,fri

! fr=fri
! fp=fsw

frx = fr
fpx = fp

dt = 1d100 ! ~3000 yr

nmx = nx*ny*2
n = nmx
nnz = 3*4 + 4*((nx-2)*2+(ny-2)*2) + 5*(nx-2)*(ny-2) + nx*ny &
    + ny + 2*ny*(nx-1) + nx*ny

if (allocated(ai)) deallocate(ai)
if (allocated(ap)) deallocate(ap)
if (allocated(ax)) deallocate(ax)
if (allocated(bx)) deallocate(bx)
if (allocated(kai)) deallocate(kai)
if (allocated(emx)) deallocate(emx)

allocate(ai(nnz))
allocate(ap(n+1))
allocate(ax(nnz))
allocate(bx(n))
allocate(kai(n))
allocate(emx(n))

time = 0d0
do it = 1,1
print*,'time', it,time
itr = 0
error = 1d4
do while (error>tol)

ai = 0
ap = 0
ax = 0d0
bx = 0d0
kai = 0d0
emx = 0d0

cnt = 0
cnt2 = 0

ap(1) = 0

rxn = frx*(1d0-fpx)-alfa*(1d0-frx)*fpx
drxn_dfr = (1d0-fpx)-alfa*(-1d0)*fpx
drxn_dfp = frx*(-1d0)-alfa*(1d0-frx)

do iy=1,ny
    do ix = 1,nx
        ! rock phase 
        row = 2*(ix-1) + 1 + (iy-1)*nx*2
        
        tmpint =0
        tmprl = 0d0
        cnt2 = 1
        
        tmpint(1) = row
        
        if (ix==1) then 
            bx(row)= &
                ! + (1d0-poro(ix,iy))*rhom*ms*(frx(ix,iy)-fr(ix,iy))/dt &
                + (1d0-poro(ix,iy))*rhom*ms*w*(frx(ix,iy)-fri)/dx(ix)   &
                + rhob(ix,iy)*ms*mw*kex(ix,iy)*rxn(ix,iy) 
            tmprl(1) = &
                ! + (1d0-poro(ix,iy))*rhom*ms*(1d0)/dt &
                + (1d0-poro(ix,iy))*rhom*ms*w*(1d0)/dx(ix)   &
                + rhob(ix,iy)*ms*mw*kex(ix,iy)*drxn_dfr(ix,iy)  
            tmprl(2) =    &
                - rhob(ix,iy)*ms*mw*kex(ix,iy)*drxn_dfr(ix,iy)
            tmprl(3) =    &
                + (1d0-poro(ix+1,iy))*rhom*ms*w*(-1d0)/dx(ix+1)
            
            tmpint(2) = row+1
            tmpint(3) = row+2
            cnt2 = cnt2 + 2
            
        else if (ix==nx) then 
            bx(row)= &
                ! + (1d0-poro(ix,iy))*rhom*ms*(frx(ix,iy)-fr(ix,iy))/dt &
                + (1d0-poro(ix,iy))*rhom*ms*w*(frx(ix,iy)-frx(ix-1,iy))/dx(ix)   &
                + rhob(ix,iy)*ms*mw*kex(ix,iy)*rxn(ix,iy) 
            tmprl(1) = &
                ! + (1d0-poro(ix,iy))*rhom*ms*(1d0)/dt &
                + (1d0-poro(ix,iy))*rhom*ms*w*(1d0)/dx(ix)   &
                + rhob(ix,iy)*ms*mw*kex(ix,iy)*drxn_dfr(ix,iy)  
            tmprl(2) =    &
                - rhob(ix,iy)*ms*mw*kex(ix,iy)*drxn_dfr(ix,iy)
            
            tmpint(2) = row+1
            cnt2 = cnt2 + 1
            
        else 
            bx(row)= &
                ! + (1d0-poro(ix,iy))*rhom*ms*(frx(ix,iy)-fr(ix,iy))/dt &
                + (1d0-poro(ix,iy))*rhom*ms*w*(frx(ix,iy)-frx(ix-1,iy))/dx(ix)   &
                + rhob(ix,iy)*ms*mw*kex(ix,iy)*rxn(ix,iy) 
            tmprl(1) = &
                ! + (1d0-poro(ix,iy))*rhom*ms*(1d0)/dt &
                + (1d0-poro(ix,iy))*rhom*ms*w*(1d0)/dx(ix)   &
                + rhob(ix,iy)*ms*mw*kex(ix,iy)*drxn_dfr(ix,iy)  
            tmprl(2) =    &
                - rhob(ix,iy)*ms*mw*kex(ix,iy)*drxn_dfr(ix,iy)
            tmprl(3) =  &
                + (1d0-poro(ix+1,iy))*rhom*ms*w*(-1d0)/dx(ix+1) 
            
            tmpint(2) = row+1
            tmpint(3) = row+2
            cnt2 = cnt2 + 2
            
        endif
        
        ! print*,cnt2
        if (cnt2 < 1 .or. cnt2 > 6) then 
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
        
        !  water phase ---
        row = row + 1
        
        tmpint =0
        tmprl = 0d0
        cnt2 = 1
        
        tmpint(1) = row
        
        if (iy == 1) then 
            if (ix==1) then 
                bx(row)= &
                    ! + poro(ix,iy)*rho(ix,iy)*mw*(fpx(ix,iy)-fp(ix,iy))/dt   &
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*(fpx(ix+1,iy)-fpx(ix,iy))/dx(ix)  &
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(fpx(ix,iy)-fsw)/dy(iy)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*(fpx(ix,iy+1)-fpx(ix,iy))/dy(iy)  & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(fpx(ix+1,iy)-1d0*fpx(ix,iy))/(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix)  & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*((fpx(ix,iy+1)-fpx(ix,iy))/(0.5d0*(dy(iy)+dy(iy+1)))  &
                        -(fpx(ix,iy)-fsw)/dy(iy))/dy(iy)  & 
                    - rhob(ix,iy)*ms*mw*kex(ix,iy)*rxn(ix,iy) 
                tmprl(1) = &
                    ! + poro(ix,iy)*rho(ix,iy)*mw*(1d0)/dt   &
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(1d0)/dx(ix)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*(-1d0)/dx(ix)  &
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(1d0)/dy(iy)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*(-1d0)/dy(iy)  & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(-1d0)/(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix)  & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(-1d0/(0.5d0*(dy(iy)+dy(iy+1)))-1d0/dy(iy))/dy(iy)  & 
                    - rhob(ix,iy)*ms*mw*kex(ix,iy)*drxn_dfp(ix,iy) 
                tmprl(2) =  & 
                    + rhob(ix,iy)*ms*mw*kex(ix,iy)*drxn_dfp(ix,iy)
                tmprl(3) =    &
                    + poro(ix+1,iy)*rho(ix+1,iy)*mw*(vmx(ix+1,iy)+abs(vmx(ix+1,iy)))*0.5d0*(-1d0)/dx(ix+1)  & 
                    - poro(ix+1,iy)*rho(ix+1,iy)*mw*dif(ix+1,iy)*(1d0)/(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix+1)  &
                    - mw*(poro(ix+1,iy)*rho(ix+1,iy)*dif(ix+1,iy)-poro(ix,iy)*rho(ix,iy)*dif(ix,iy))*(-1d0)  &
                        /(0.5d0*(dx(ix)+dx(Ix+1)))/dx(ix+1)
                tmprl(4) = &
                    + poro(ix,iy+1)*rho(ix,iy+1)*mw*(vmy(ix,iy+1)+abs(vmy(ix,iy+1)))*0.5d0*(-1d0)/dy(iy+1)  &
                    - poro(ix,iy+1)*rho(ix,iy+1)*mw*dif(ix,iy+1)*(1d0)/(0.5d0*(dy(iy)+dy(iy+1)))/dy(iy+1)  &
                    - mw*(poro(ix,iy+1)*rho(ix,iy+1)*dif(ix,iy+1)-poro(ix,iy)*rho(ix,iy)*dif(ix,iy))*(-1d0)  &
                        /(0.5d0*(dy(iy)+dy(iy+1)))/dy(iy+1) 
                
                tmpint(2) = row-1 
                tmpint(3) = row+2
                tmpint(4) = row+2*nx
                cnt2 = cnt2 + 3  
                
            else if (ix==nx) then 
                bx(row)= &
                    ! + poro(ix,iy)*rho(ix,iy)*mw*(fpx(ix,iy)-fp(ix,iy))/dt   &
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(fpx(ix,iy)-fpx(ix-1,iy))/dx(ix)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(fpx(ix,iy)-fsw)/dy(iy)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*(fpx(ix,iy+1)-fpx(ix,iy))/dy(iy)  & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(fpx(ix-1,iy)-1d0*fpx(ix,iy))/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix-1,iy)*rho(ix-1,iy)*dif(ix-1,iy)) &
                        *(fpx(ix,iy)-fpx(ix-1,iy))/(0.5d0*(dx(ix)+dx(Ix-1)))/dx(ix)  &
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*((fpx(ix,iy+1)-fpx(ix,iy))/(0.5d0*(dy(iy)+dy(iy+1)))  &
                        -(fpx(ix,iy)-fsw)/dy(iy))/dy(iy)  &  
                    - rhob(ix,iy)*ms*mw*kex(ix,iy)*rxn(ix,iy) 
                tmprl(1) = &
                    ! + poro(ix,iy)*rho(ix,iy)*mw*(1d0)/dt   &
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(1d0)/dx(ix)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(1d0)/dy(iy)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*(-1d0)/dy(iy)  & 
                    - poro(ix,iy)*rho(ix,iy)*mw*(-1d0)/(0.5d0*(dx(ix)+dx(Ix-1)))/dx(ix)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix-1,iy)*rho(ix-1,iy)*dif(ix-1,iy))*(1d0)  &
                        /(0.5d0*(dx(ix)+dx(Ix-1)))/dx(ix)  &
                    - poro(ix,iy)*rho(ix,iy)*mw*(-1d0/(0.5d0*(dy(iy)+dy(iy+1)))-1d0/dy(iy))/dy(iy)  & 
                    - rhob(ix,iy)*ms*mw*kex(ix,iy)*drxn_dfp(ix,iy) 
                tmprl(2) =  & 
                    + rhob(ix,iy)*ms*mw*kex(ix,iy)*drxn_dfp(ix,iy)
                tmprl(3) =    &
                    + poro(ix-1,iy)*rho(ix-1,iy)*mw*(vmx(ix-1,iy)-abs(vmx(ix-1,iy)))*0.5d0*(1d0)/dx(ix-1)  & 
                    - poro(ix-1,iy)*rho(ix-1,iy)*mw*dif(ix-1,iy)*(1d0)/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix-1)  
                tmprl(4) = &
                    + poro(ix,iy+1)*rho(ix,iy+1)*mw*(vmy(ix,iy+1)+abs(vmy(ix,iy+1)))*0.5d0*(-1d0)/dy(iy+1)  &
                    - poro(ix,iy+1)*rho(ix,iy+1)*mw*dif(ix,iy+1)*(1d0)/(0.5d0*(dy(iy)+dy(iy+1)))/dy(iy+1)  &
                    - mw*(poro(ix,iy+1)*rho(ix,iy+1)*dif(ix,iy+1)-poro(ix,iy)*rho(ix,iy)*dif(ix,iy))*(-1d0)  &
                        /(0.5d0*(dy(iy)+dy(iy+1)))/dy(iy+1)  
                
                tmpint(2) = row-1
                tmpint(3) = row-2
                tmpint(4) = row+2*nx
                cnt2 = cnt2 + 3  
                
            else 
                bx(row)= &
                    ! + poro(ix,iy)*rho(ix,iy)*mw*(fpx(ix,iy)-fp(ix,iy))/dt   &
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(fpx(ix,iy)-fpx(ix-1,iy))/dx(ix)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*(fpx(ix+1,iy)-fpx(ix,iy))/dx(ix)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(fpx(ix,iy)-fsw)/dy(iy)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*(fpx(ix,iy+1)-fpx(ix,iy))/dy(iy)  & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*((fpx(ix+1,iy)-fpx(ix,iy))/(0.5d0*(dx(ix)+dx(Ix+1)))  &
                        -(fpx(ix,iy)-fpx(ix-1,iy))/(0.5d0*(dx(ix)+dx(Ix-1))))/dx(ix)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix-1,iy)*rho(ix-1,iy)*dif(ix-1,iy))  &
                        *(fpx(ix,iy)-fpx(ix-1,iy))/(0.5d0*(dx(ix)+dx(Ix-1)))/dx(ix)  &
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*((fpx(ix,iy+1)-fpx(ix,iy))/(0.5d0*(dy(iy)+dy(iy+1)))  &
                        -(fpx(ix,iy)-fsw)/dy(iy))/dy(iy)  & 
                    - rhob(ix,iy)*ms*mw*kex(ix,iy)*rxn(ix,iy) 
                tmprl(1) = &
                    ! + poro(ix,iy)*rho(ix,iy)*mw*(1d0)/dt   &
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(1d0)/dx(ix)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*(-1d0)/dx(ix)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(1d0)/dy(iy)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*(-1d0)/dy(iy)  & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(-1d0/(0.5d0*(dx(ix)+dx(Ix+1)))-1d0/(0.5d0*(dx(Ix)+dx(ix-1))))/dx(ix)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix-1,iy)*rho(ix-1,iy)*dif(ix-1,iy))*(1d0)  &
                        /(0.5d0*(dx(ix)+dx(Ix-1)))/dx(ix)  &
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(-1d0/(0.5d0*(dy(iy)+dy(iy+1)))-1d0/dy(iy))/dy(iy)  & 
                    - rhob(ix,iy)*ms*mw*kex(ix,iy)*drxn_dfp(ix,iy) 
                tmprl(2) =  & 
                    + rhob(ix,iy)*ms*mw*kex(ix,iy)*drxn_dfp(ix,iy)
                tmprl(3) =    &
                    + poro(ix-1,iy)*rho(ix-1,iy)*mw*(vmx(ix-1,iy)-abs(vmx(ix-1,iy)))*0.5d0*(1d0)/dx(ix-1)  & 
                    - poro(ix-1,iy)*rho(ix-1,iy)*mw*dif(ix-1,iy)*(1d0)/(0.5d0*(dx(ix)+dx(Ix-1)))/dx(ix-1)  
                tmprl(4) =    &
                    + poro(ix+1,iy)*rho(ix+1,iy)*mw*(vmx(ix+1,iy)+abs(vmx(ix+1,iy)))*0.5d0*(-1d0)/dx(ix+1)  & 
                    - poro(ix+1,iy)*rho(ix+1,iy)*mw*dif(ix+1,iy)*(1d0)/(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix+1)  &
                    - mw*(poro(ix+1,iy)*rho(ix+1,iy)*dif(ix+1,iy)-poro(ix,iy)*rho(ix,iy)*dif(ix,iy))*(-1d0)  &
                        /(0.5d0*(dx(ix)+dx(Ix+1)))/dx(ix+1)  
                tmprl(5) = &
                    + poro(ix,iy+1)*rho(ix,iy+1)*mw*(vmy(ix,iy+1)+abs(vmy(ix,iy+1)))*0.5d0*(-1d0)/dy(iy+1)  &
                    - poro(ix,iy+1)*rho(ix,iy+1)*mw*dif(ix,iy+1)*(1d0)/(0.5d0*(dy(iy)+dy(iy+1)))/dy(iy+1)  &
                    - mw*(poro(ix,iy+1)*rho(ix,iy+1)*dif(ix,iy+1)-poro(ix,iy)*rho(ix,iy)*dif(ix,iy))*(-1d0)  &
                        /(0.5d0*(dy(iy)+dy(iy+1)))/dy(iy+1)  
                
                tmpint(2) = row-1
                tmpint(3) = row-2
                tmpint(4) = row+2
                tmpint(5) = row+2*nx
                cnt2 = cnt2 + 4   
            
            endif
        
        else if (iy == ny) then 
            if (ix==1) then 
                bx(row)= &
                    ! + poro(ix,iy)*rho(ix,iy)*mw*(fpx(ix,iy)-fp(ix,iy))/dt   &
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*(fpx(ix+1,iy)-fpx(ix,iy))/dx(ix)  &
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(fpx(ix,iy)-fpx(ix,iy-1))/dy(iy)  & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(fpx(ix+1,iy)-1d0*fpx(ix,iy))/(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix)  & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(fpx(ix,iy-1)-1d0*fpx(ix,iy))/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix,iy-1)*rho(ix,iy-1)*dif(ix,iy-1))  &
                        *(fpx(ix,iy)-fpx(ix,iy-1))/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)  &
                    - rhob(ix,iy)*ms*mw*kex(ix,iy)*rxn(ix,iy) 
                tmprl(1) = &
                    ! + poro(ix,iy)*rho(ix,iy)*mw*(1d0)/dt   &
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*(-1d0)/dx(ix)  &
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(1d0)/dy(iy)  & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(-1d0)/(0.5d0*(dx(ix)+dx(Ix+1)))/dx(ix)  & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(-1d0)/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix,iy-1)*rho(ix,iy-1)*dif(ix,iy-1))*(1d0)  &
                        /(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)  &
                    - rhob(ix,iy)*ms*mw*kex(ix,iy)*drxn_dfp(ix,iy) 
                tmprl(2) =  & 
                    + rhob(ix,iy)*ms*mw*kex(ix,iy)*drxn_dfp(ix,iy)
                tmprl(3) =    &
                    + poro(ix+1,iy)*rho(ix+1,iy)*mw*(vmx(ix+1,iy)+abs(vmx(ix+1,iy)))*0.5d0*(-1d0)/dx(ix+1)  & 
                    - poro(ix+1,iy)*rho(ix+1,iy)*mw*dif(ix+1,iy)*(1d0)/(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix+1)  &
                    - mw*(poro(ix+1,iy)*rho(ix+1,iy)*dif(ix+1,iy)-poro(ix,iy)*rho(ix,iy)*dif(ix,iy))*(-1d0)   &
                        /(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix+1)  
                tmprl(4) = &
                    + poro(ix,iy-1)*rho(ix,iy-1)*mw*(vmy(ix,iy-1)-abs(vmy(ix,iy-1)))*0.5d0*(1d0)/dy(iy-1)  &
                    - poro(ix,iy-1)*rho(ix,iy-1)*mw*dif(ix,iy-1)*(1d0)/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy-1)  
                
                tmpint(2) = row-1
                tmpint(3) = row+2
                tmpint(4) = row-2*nx
                cnt2 = cnt2 + 3  
                
            else if (ix==nx) then 
                bx(row)= &
                    ! + poro(ix,iy)*rho(ix,iy)*mw*(fpx(ix,iy)-fp(ix,iy))/dt   &
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(fpx(ix,iy)-fpx(ix-1,iy))/dx(ix)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(fpx(ix,iy)-fpx(ix,iy-1))/dy(iy)  & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(fpx(ix-1,iy)-1d0*fpx(ix,iy))/(0.5d0*(dx(ix)+dx(Ix-1)))/dx(ix)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix-1,iy)*rho(ix-1,iy)*dif(ix-1,iy)) &
                        *(fpx(ix,iy)-fpx(ix-1,iy))/(0.5d0*(dx(ix)+dx(Ix-1)))/dx(ix)  &
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(fpx(ix,iy-1)-1d0*fpx(ix,iy))/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix,iy-1)*rho(ix,iy-1)*dif(ix,iy-1))  &
                        *(fpx(ix,iy)-fpx(ix,iy-1))/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)  &
                    - rhob(ix,iy)*ms*mw*kex(ix,iy)*rxn(ix,iy) 
                tmprl(1) = &
                    ! + poro(ix,iy)*rho(ix,iy)*mw*(1d0)/dt   &
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(1d0)/dx(ix)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(1d0)/dy(iy)  & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(-1d0)/(0.5d0*(dx(ix)+dx(Ix-1)))/dx(ix)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix-1,iy)*rho(ix-1,iy)*dif(ix-1,iy))*(1d0)  &
                        /(0.5d0*(dx(ix)+dx(Ix-1)))/dx(ix)  &
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(-1d0)/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix,iy-1)*rho(ix,iy-1)*dif(ix,iy-1))*(1d0)  &
                        /(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)  &
                    - rhob(ix,iy)*ms*mw*kex(ix,iy)*drxn_dfp(ix,iy) 
                tmprl(2) =  & 
                    + rhob(ix,iy)*ms*mw*kex(ix,iy)*drxn_dfp(ix,iy)
                tmprl(3) =    &
                    + poro(ix-1,iy)*rho(ix-1,iy)*mw*(vmx(ix-1,iy)-abs(vmx(ix-1,iy)))*0.5d0*(1d0)/dx(ix-1)  & 
                    - poro(ix-1,iy)*rho(ix-1,iy)*mw*dif(ix-1,iy)*(1d0)/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix-1)  
                tmprl(4) = &
                    + poro(ix,iy-1)*rho(ix,iy-1)*mw*(vmy(ix,iy-1)-abs(vmy(ix,iy-1)))*0.5d0*(1d0)/dy(iy-1)  &
                    - poro(ix,iy-1)*rho(ix,iy-1)*mw*dif(ix,iy-1)*(1d0)/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy-1)  
                
                tmpint(2) = row-1
                tmpint(3) = row-2
                tmpint(4) = row-2*nx
                cnt2 = cnt2 + 3  
                
            else 
                bx(row)= &
                    ! + poro(ix,iy)*rho(ix,iy)*mw*(fpx(ix,iy)-fp(ix,iy))/dt   &
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(fpx(ix,iy)-fpx(ix-1,iy))/dx(ix)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*(fpx(ix+1,iy)-fpx(ix,iy))/dx(ix)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(fpx(ix,iy)-fpx(ix,iy-1))/dy(iy)  & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*((fpx(ix+1,iy)-fpx(ix,iy))/(0.5d0*(dx(ix)+dx(ix+1)))  &
                        -(fpx(ix,iy)-fpx(ix-1,iy))/(0.5d0*(dx(ix)+dx(Ix-1))))/dx(ix)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix-1,iy)*rho(ix-1,iy)*dif(ix-1,iy))  &
                        *(fpx(ix,iy)-fpx(ix-1,iy))/(0.5d0*(dx(ix)+dx(Ix-1)))/dx(ix)  &
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(fpx(ix,iy-1)-1d0*fpx(ix,iy))/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix,iy-1)*rho(ix,iy-1)*dif(ix,iy-1))  &
                        *(fpx(ix,iy)-fpx(ix,iy-1))/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)  &
                    - rhob(ix,iy)*ms*mw*kex(ix,iy)*rxn(ix,iy) 
                tmprl(1) = &
                    ! + poro(ix,iy)*rho(ix,iy)*mw*(1d0)/dt   &
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(1d0)/dx(ix)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*(-1d0)/dx(ix)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(1d0)/dy(iy)  & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(-1d0/(0.5d0*(dx(ix)+dx(Ix+1)))-1d0/(0.5d0*(dx(ix)+dx(ix-1))))/dx(ix)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix-1,iy)*rho(ix-1,iy)*dif(ix-1,iy))*(1d0)  &
                        /(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix)  &
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(-1d0)/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix,iy-1)*rho(ix,iy-1)*dif(ix,iy-1))*(1d0)  &
                        /(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)  &
                    - rhob(ix,iy)*ms*mw*kex(ix,iy)*drxn_dfp(ix,iy) 
                tmprl(2) =  & 
                    + rhob(ix,iy)*ms*mw*kex(ix,iy)*drxn_dfp(ix,iy)
                tmprl(3) =    &
                    + poro(ix-1,iy)*rho(ix-1,iy)*mw*(vmx(ix-1,iy)-abs(vmx(ix-1,iy)))*0.5d0*(1d0)/dx(ix-1)  & 
                    - poro(ix-1,iy)*rho(ix-1,iy)*mw*dif(ix-1,iy)*(1d0)/(0.5d0*(dx(ix)+dx(Ix-1)))/dx(ix-1)  
                tmprl(4) =    &
                    + poro(ix+1,iy)*rho(ix+1,iy)*mw*(vmx(ix+1,iy)+abs(vmx(ix+1,iy)))*0.5d0*(-1d0)/dx(ix+1)  & 
                    - poro(ix+1,iy)*rho(ix+1,iy)*mw*dif(ix+1,iy)*(1d0)/(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix+1)  &
                    - mw*(poro(ix+1,iy)*rho(ix+1,iy)*dif(ix+1,iy)-poro(ix,iy)*rho(ix,iy)*dif(ix,iy))*(-1d0)  &
                        /(0.5d0*(dx(ix)+dx(Ix+1)))/dx(ix+1) 
                tmprl(5) = &
                    + poro(ix,iy-1)*rho(ix,iy-1)*mw*(vmy(ix,iy-1)-abs(vmy(ix,iy-1)))*0.5d0*(1d0)/dy(iy-1)  &
                    - poro(ix,iy-1)*rho(ix,iy-1)*mw*dif(ix,iy-1)*(1d0)/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy-1)  
                
                tmpint(2) = row-1
                tmpint(3) = row-2
                tmpint(4) = row+2
                tmpint(5) = row-2*nx
                cnt2 = cnt2 + 4 
            
            endif
        
        else 
            if (ix==1) then 
                bx(row)= &
                    ! + poro(ix,iy)*rho(ix,iy)*mw*(fpx(ix,iy)-fp(ix,iy))/dt   & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*(fpx(ix+1,iy)-fpx(ix,iy))/dx(ix)  &
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(fpx(ix,iy)-fpx(ix,iy-1))/dy(iy)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*(fpx(ix,iy+1)-fpx(ix,iy))/dy(iy)  & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(fpx(ix+1,iy)-1d0*fpx(ix,iy))/(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix)  & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*((fpx(ix,iy+1)-fpx(ix,iy))/(0.5d0*(dy(iy)+dy(iy+1)))  &
                        -(fpx(ix,iy)-fpx(ix,iy-1))/(0.5d0*(dy(iy)+dy(iy-1))))/dy(iy)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix,iy-1)*rho(ix,iy-1)*dif(ix,iy-1))  &
                        *(fpx(ix,iy)-fpx(ix,iy-1))/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)  &
                    - rhob(ix,iy)*ms*mw*kex(ix,iy)*rxn(ix,iy) 
                tmprl(1) = &
                    ! + poro(ix,iy)*rho(ix,iy)*mw*(1d0)/dt   &
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*(-1d0)/dx(ix)  &
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(1d0)/dy(iy)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*(-1d0)/dy(iy)  & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(-1d0)/(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix)  & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(-1d0/(0.5d0*(dy(iy)+dy(iy+1)))-1d0/(0.5d0*(dy(iy)+dy(iy-1))))/dy(iy)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix,iy-1)*rho(ix,iy-1)*dif(ix,iy-1))*(1d0)  &
                        /(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)  &
                    - rhob(ix,iy)*ms*mw*kex(ix,iy)*drxn_dfp(ix,iy) 
                tmprl(2) =  & 
                    + rhob(ix,iy)*ms*mw*kex(ix,iy)*drxn_dfp(ix,iy)
                tmprl(3) =    &
                    + poro(ix+1,iy)*rho(ix+1,iy)*mw*(vmx(ix+1,iy)+abs(vmx(ix+1,iy)))*0.5d0*(-1d0)/dx(ix+1)  & 
                    - poro(ix+1,iy)*rho(ix+1,iy)*mw*dif(ix+1,iy)*(1d0)/(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix+1)  &
                    - mw*(poro(ix+1,iy)*rho(ix+1,iy)*dif(ix+1,iy)-poro(ix,iy)*rho(ix,iy)*dif(ix,iy))*(-1d0)  &
                        /(0.5d0*(dx(ix)+dx(Ix+1)))/dx(ix+1) 
                tmprl(4) = &
                    + poro(ix,iy-1)*rho(ix,iy-1)*mw*(vmy(ix,iy-1)-abs(vmy(ix,iy-1)))*0.5d0*(1d0)/dy(iy-1)  &
                    - poro(ix,iy-1)*rho(ix,iy-1)*mw*dif(ix,iy-1)*(1d0)/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy-1)  
                tmprl(5) = &
                    + poro(ix,iy+1)*rho(ix,iy+1)*mw*(vmy(ix,iy+1)+abs(vmy(ix,iy+1)))*0.5d0*(-1d0)/dy(iy+1)  &
                    - poro(ix,iy+1)*rho(ix,iy+1)*mw*dif(ix,iy+1)*(1d0)/(0.5d0*(dy(iy)+dy(iy+1)))/dy(iy+1)  &
                    - mw*(poro(ix,iy+1)*rho(ix,iy+1)*dif(ix,iy+1)-poro(ix,iy)*rho(ix,iy)*dif(ix,iy))*(-1d0)  &
                        /(0.5d0*(dy(iy)+dy(iy+1)))/dy(iy+1)  
                
                tmpint(2) = row-1
                tmpint(3) = row+2
                tmpint(4) = row-2*nx
                tmpint(5) = row+2*nx
                cnt2 = cnt2 + 4
                
            else if (ix==nx) then 
                bx(row)= &
                    ! + poro(ix,iy)*rho(ix,iy)*mw*(fpx(ix,iy)-fp(ix,iy))/dt   &
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(fpx(ix,iy)-fpx(ix-1,iy))/dx(ix)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(fpx(ix,iy)-fpx(ix,iy-1))/dy(iy)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*(fpx(ix,iy+1)-fpx(ix,iy))/dy(iy)  & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(fpx(ix-1,iy)-1d0*fpx(ix,iy))/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix-1,iy)*rho(ix-1,iy)*dif(ix-1,iy))  &
                        *(fpx(ix,iy)-fpx(ix-1,iy))/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix)  &
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*((fpx(ix,iy+1)-fpx(ix,iy))/(0.5d0*(dy(iy)+dy(iy+1)))  &
                        -(fpx(ix,iy)-fpx(ix,iy-1))/(0.5d0*(dy(iy)+dy(iy-1))))/dy(iy)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix,iy-1)*rho(ix,iy-1)*dif(ix,iy-1))  &
                        *(fpx(ix,iy)-fpx(ix,iy-1))/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)  &
                    - rhob(ix,iy)*ms*mw*kex(ix,iy)*rxn(ix,iy) 
                tmprl(1) = &
                    ! + poro(ix,iy)*rho(ix,iy)*mw*(1d0)/dt   &
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(1d0)/dx(ix)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(1d0)/dy(iy)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*(-1d0)/dy(iy)  & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(-1d0)/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix-1,iy)*rho(ix-1,iy)*dif(ix-1,iy))*(1d0)   &
                        /(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix)  &
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(-1d0/(0.5d0*(dy(iy)+dy(iy+1)))-1d0/(0.5d0*(dy(iy)+dy(iy-1))))/dy(iy)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix,iy-1)*rho(ix,iy-1)*dif(ix,iy-1))*(1d0)  &
                        /(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)  &
                    - rhob(ix,iy)*ms*mw*kex(ix,iy)*drxn_dfp(ix,iy) 
                tmprl(2) =  & 
                    + rhob(ix,iy)*ms*mw*kex(ix,iy)*drxn_dfp(ix,iy)
                tmprl(3) =    &
                    + poro(ix-1,iy)*rho(ix-1,iy)*mw*(vmx(ix-1,iy)-abs(vmx(ix-1,iy)))*0.5d0*(1d0)/dx(ix-1)  & 
                    - poro(ix-1,iy)*rho(ix-1,iy)*mw*dif(ix-1,iy)*(1d0)/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix-1)  
                tmprl(4) = &
                    + poro(ix,iy-1)*rho(ix,iy-1)*mw*(vmy(ix,iy-1)-abs(vmy(ix,iy-1)))*0.5d0*(1d0)/dy(iy-1)  &
                    - poro(ix,iy-1)*rho(ix,iy-1)*mw*dif(ix,iy-1)*(1d0)/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy-1)  
                tmprl(5) = &
                    + poro(ix,iy+1)*rho(ix,iy+1)*mw*(vmy(ix,iy+1)+abs(vmy(ix,iy+1)))*0.5d0*(-1d0)/dy(iy+1)  &
                    - poro(ix,iy+1)*rho(ix,iy+1)*mw*dif(ix,iy+1)*(1d0)/(0.5d0*(dy(iy)+dy(iy+1)))/dy(iy+1)  &
                    - mw*(poro(ix,iy+1)*rho(ix,iy+1)*dif(ix,iy+1)-poro(ix,iy)*rho(ix,iy)*dif(ix,iy))*(-1d0)  &
                        /(0.5d0*(dy(iy)+dy(iy+1)))/dy(iy+1)  
                
                tmpint(2) = row-1
                tmpint(3) = row-2
                tmpint(4) = row-2*nx
                tmpint(5) = row+2*nx
                cnt2 = cnt2 + 4  
                
            else 
                bx(row)= &
                    ! + poro(ix,iy)*rho(ix,iy)*mw*(fpx(ix,iy)-fp(ix,iy))/dt   &
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(fpx(ix,iy)-fpx(ix-1,iy))/dx(ix)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*(fpx(ix+1,iy)-fpx(ix,iy))/dx(ix)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(fpx(ix,iy)-fpx(ix,iy-1))/dy(iy)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*(fpx(ix,iy+1)-fpx(ix,iy))/dy(iy)  & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*((fpx(ix+1,iy)-fpx(ix,iy))/(0.5d0*(dx(ix)+dx(ix+1)))  &
                        -(fpx(ix,iy)-fpx(ix-1,iy))/(0.5d0*(dx(ix)+dx(ix-1))))/dx(ix)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix-1,iy)*rho(ix-1,iy)*dif(ix-1,iy))  &
                        *(fpx(ix,iy)-fpx(ix-1,iy))/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix)  &
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*((fpx(ix,iy+1)-fpx(ix,iy))/(0.5d0*(dy(iy)+dy(iy+1)))  &
                        -(fpx(ix,iy)-fpx(ix,iy-1))/(0.5d0*(dy(iy)+dy(iy-1))))/dy(iy)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix,iy-1)*rho(ix,iy-1)*dif(ix,iy-1))  &
                        *(fpx(ix,iy)-fpx(ix,iy-1))/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)  &
                    - rhob(ix,iy)*ms*mw*kex(ix,iy)*rxn(ix,iy) 
                tmprl(1) = &
                    ! + poro(ix,iy)*rho(ix,iy)*mw*(1d0)/dt   &
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(1d0)/dx(ix)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*(-1d0)/dx(ix)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(1d0)/dy(iy)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*(-1d0)/dy(iy)  & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(-1d0/(0.5d0*(dx(ix)+dx(ix+1)))-1d0/(0.5d0*(dx(ix)+dx(ix-1))))/dx(ix)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix-1,iy)*rho(ix-1,iy)*dif(ix-1,iy))*(1d0)   &
                        /(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix)  &
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(-1d0/(0.5d0*(dy(iy)+dy(iy+1)))-1d0/(0.5d0*(dy(iy)+dy(iy-1))))/dy(iy)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix,iy-1)*rho(ix,iy-1)*dif(ix,iy-1))*(1d0)  &
                        /(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)  &
                    - rhob(ix,iy)*ms*mw*kex(ix,iy)*drxn_dfp(ix,iy) 
                tmprl(2) =  & 
                    + rhob(ix,iy)*ms*mw*kex(ix,iy)*drxn_dfp(ix,iy)
                tmprl(3) =    &
                    + poro(ix-1,iy)*rho(ix-1,iy)*mw*(vmx(ix-1,iy)-abs(vmx(ix-1,iy)))*0.5d0*(1d0)/dx(ix-1)  & 
                    - poro(ix-1,iy)*rho(ix-1,iy)*mw*dif(ix-1,iy)*(1d0)/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix-1)  
                tmprl(4) =    &
                    + poro(ix+1,iy)*rho(ix+1,iy)*mw*(vmx(ix+1,iy)+abs(vmx(ix+1,iy)))*0.5d0*(-1d0)/dx(ix+1)  & 
                    - poro(ix+1,iy)*rho(ix+1,iy)*mw*dif(ix+1,iy)*(1d0)/(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix+1)  &
                    - mw*(poro(ix+1,iy)*rho(ix+1,iy)*dif(ix+1,iy)-poro(ix,iy)*rho(ix,iy)*dif(ix,iy))*(-1d0)  &
                        /(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix+1) 
                tmprl(5) = &
                    + poro(ix,iy-1)*rho(ix,iy-1)*mw*(vmy(ix,iy-1)-abs(vmy(ix,iy-1)))*0.5d0*(1d0)/dy(iy-1)  &
                    - poro(ix,iy-1)*rho(ix,iy-1)*mw*dif(ix,iy-1)*(1d0)/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy-1)  
                tmprl(6) = &
                    + poro(ix,iy+1)*rho(ix,iy+1)*mw*(vmy(ix,iy+1)+abs(vmy(ix,iy+1)))*0.5d0*(-1d0)/dy(iy+1)  &
                    + poro(ix,iy+1)*rho(ix,iy+1)*mw*dif(ix,iy+1)*(1d0)/(0.5d0*(dy(iy)+dy(iy+1)))/dy(iy+1)  &
                    - mw*(poro(ix,iy+1)*rho(ix,iy+1)*dif(ix,iy+1)-poro(ix,iy)*rho(ix,iy)*dif(ix,iy))*(-1d0)  &
                        /(0.5d0*(dy(iy)+dy(iy+1)))/dy(iy+1)  
                
                tmpint(2) = row-1
                tmpint(3) = row-2
                tmpint(4) = row+2
                tmpint(5) = row-2*nx
                tmpint(6) = row+2*nx
                cnt2 = cnt2 + 5  
            
            endif 
            
        endif
        
        ! print*,cnt2
        if (cnt2 < 1 .or. cnt2 > 6) then 
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
    do ix=1,nx
        row = 2*(ix-1) + 1 + (iy-1)*nx*2
        ! if (kai(row)>frx(ix,iy))kai(row)=0.5d0*frx(ix,iy)
        ! if (kai(row)<-frx(ix,iy))kai(row)=-0.5d0*frx(ix,iy)
        ! if (kai(row+1)>fpx(ix,iy))kai(row+1)=0.5d0*fpx(ix,iy)
        ! if (kai(row+1)<-fpx(ix,iy))kai(row+1)=-0.5d0*fpx(ix,iy)
        frx(ix,iy) = frx(ix,iy) + kai(row) 
        fpx(ix,iy) = fpx(ix,iy) + kai(row+1) 
        if (frx(ix,iy)/=0d0) emx(row) = abs(kai(row)/frx(ix,iy))
        if (fpx(ix,iy)/=0d0) emx(row+1) = abs(kai(row+1)/fpx(ix,iy))
    enddo
enddo


error = maxval(emx)

print*,'iteration, error',itr, error
itr = itr+1
enddo

! it = it + 1
time = time + dt

enddo

do iy=1,ny
    do ix=1,nx
        o18r(ix,iy) = f2d(frx(ix,iy))
        o18p(ix,iy) = f2d(fpx(ix,iy))
    enddo
enddo

omega = frx*(1d0-fpx)/(1d0-frx)/fpx/alfa

open(unit=100,file=trim(adjustl(workdir))//'o18r-'//trim(adjustl(intsw))//'.txt',action='write',status='replace')
open(unit=200,file=trim(adjustl(workdir))//'o18p-'//trim(adjustl(intsw))//'.txt',action='write',status='replace')
open(unit=300,file=trim(adjustl(workdir))//'omage-'//trim(adjustl(intsw))//'.txt',action='write',status='replace')
do iy=1,ny
    write(100,*) (o18r(ix,iy),ix=1,nx) 
    write(200,*) (o18p(ix,iy),ix=1,nx) 
    write(300,*) (log10(omega(ix,iy)),ix=1,nx) 
enddo
close(100)
close(200)
close(300)

rflxadv = 0d0
rflxrxn = 0d0 
rflxt = 0d0
pflxadv = 0d0 
pflxrxn = 0d0 
pflxdif = 0d0
pflxt = 0d0

rflxrxnpls =0d0 
rflxrxnmns =0d0

do iy=1,ny
    do ix = 1,nx
        ! rock phase 
        
        if (frx(ix,iy)>=fri) rflxrxnpls(ix,iy) = rflxrxnpls(ix,iy) + rhob(ix,iy)*ms*mw*kex(ix,iy)*rxn(ix,iy)
        if (frx(ix,iy)< fri) rflxrxnmns(ix,iy) = rflxrxnmns(ix,iy) + rhob(ix,iy)*ms*mw*kex(ix,iy)*rxn(ix,iy)
        
        if (ix==1) then 
            ! rflxt = rflxt + &
                ! + (1d0-poro(ix,iy))*rhom*ms*(frx(ix,iy)-fr(ix,iy))/dt &
            
            rflxadv(ix,iy) = rflxadv(ix,iy) &
                + (1d0-poro(ix,iy))*rhom*ms*w*(frx(ix,iy)-fri)/dx(ix) 
            
            rflxrxn(ix,iy) = rflxrxn(ix,iy) &
                + rhob(ix,iy)*ms*mw*kex(ix,iy)*rxn(ix,iy)
            
        else if (ix==nx) then 
            ! rflxt = rflxt &
                ! + (1d0-poro(ix,iy))*rhom*ms*(frx(ix,iy)-fr(ix,iy))/dt &
            
            rflxadv(ix,iy) = rflxadv(ix,iy)  &
                + (1d0-poro(ix,iy))*rhom*ms*w*(frx(ix,iy)-frx(ix-1,iy))/dx(ix)   
                
            rflxrxn(ix,iy) = rflxrxn(ix,iy) &
                + rhob(ix,iy)*ms*mw*kex(ix,iy)*rxn(ix,iy) 
            
        else 
            ! rflxt = rflxt &
                ! + (1d0-poro(ix,iy))*rhom*ms*(frx(ix,iy)-fr(ix,iy))/dt &
                
            rflxadv(ix,iy) = rflxadv(ix,iy) &
                + (1d0-poro(ix,iy))*rhom*ms*w*(frx(ix,iy)-frx(ix-1,iy))/dx(ix)   
            
            rflxrxn(ix,iy) = rflxrxn(ix,iy) & 
                + rhob(ix,iy)*ms*mw*kex(ix,iy)*rxn(ix,iy) 
            
        endif
        
        !  water phase ---
        
        if (iy == 1) then 
            if (ix==1) then 
                ! pflxt = pflxt &
                    ! + poro(ix,iy)*rho(ix,iy)*mw*(fpx(ix,iy)-fp(ix,iy))/dt   &
                
                pflxadv(ix,iy) = pflxadv(ix,iy) &
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*(fpx(ix+1,iy)-fpx(ix,iy))/dx(ix)  &
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(fpx(ix,iy)-fsw)/dy(iy)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*(fpx(ix,iy+1)-fpx(ix,iy))/dy(iy)  
                
                pflxdif(ix,iy) = pflxdif(ix,iy) & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(fpx(ix+1,iy)-1d0*fpx(ix,iy))/(0.5d0*(dx(ix)+dx(Ix+1)))/dx(ix)  & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*((fpx(ix,iy+1)-fpx(ix,iy))/(0.5d0*(dy(iy)+dy(iy+1)))  &
                        -(fpx(ix,iy)-fsw)/dy(iy))/dy(iy)  
                
                pflxrxn(ix,iy) = pflxrxn(ix,iy) & 
                    - rhob(ix,iy)*ms*mw*kex(ix,iy)*rxn(ix,iy) 
                
            else if (ix==nx) then 
                ! pflxt = pflxt &
                    ! + poro(ix,iy)*rho(ix,iy)*mw*(fpx(ix,iy)-fp(ix,iy))/dt   &
                
                pflxadv(ix,iy) = pflxadv(ix,iy) &
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(fpx(ix,iy)-fpx(ix-1,iy))/dx(ix)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(fpx(ix,iy)-fsw)/dy(iy)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*(fpx(ix,iy+1)-fpx(ix,iy))/dy(iy)  
                
                pflxdif(ix,iy) = pflxdif(ix,iy) & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(fpx(ix-1,iy)-1d0*fpx(ix,iy))/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix-1,iy)*rho(ix-1,iy)*dif(ix-1,iy)) &
                        *(fpx(ix,iy)-fpx(ix-1,iy))/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix)  &
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*((fpx(ix,iy+1)-fpx(ix,iy))/(0.5d0*(dy(iy)+dy(iy+1)))  &
                        -(fpx(ix,iy)-fsw)/dy(iy))/dy(iy)  
                
                pflxrxn(ix,iy) = pflxrxn(ix,iy) & 
                    - rhob(ix,iy)*ms*mw*kex(ix,iy)*rxn(ix,iy) 
                
            else 
                ! pflxt = pflxt &
                    ! + poro(ix,iy)*rho(ix,iy)*mw*(fpx(ix,iy)-fp(ix,iy))/dt   &
                
                pflxadv(ix,iy) = pflxadv(ix,iy) &
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(fpx(ix,iy)-fpx(ix-1,iy))/dx(ix)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*(fpx(ix+1,iy)-fpx(ix,iy))/dx(ix)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(fpx(ix,iy)-fsw)/dy(iy)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*(fpx(ix,iy+1)-fpx(ix,iy))/dy(iy)  
                
                pflxdif(ix,iy) = pflxdif(ix,iy) & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*((fpx(ix+1,iy)-fpx(ix,iy))/(0.5d0*(dx(ix)+dx(Ix+1)))  &
                        -(fpx(ix,iy)-fpx(ix-1,iy))/(0.5d0*(dx(ix)+dx(ix-1))))/dx(ix)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix-1,iy)*rho(ix-1,iy)*dif(ix-1,iy))  &
                        *(fpx(ix,iy)-fpx(ix-1,iy))/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix)  &
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*((fpx(ix,iy+1)-fpx(ix,iy))/(0.5d0*(dy(iy)+dy(iy+1)))  &
                        -(fpx(ix,iy)-fsw)/dy(iy))/dy(iy)  
                
                pflxrxn(ix,iy) = pflxrxn(ix,iy) & 
                    - rhob(ix,iy)*ms*mw*kex(ix,iy)*rxn(ix,iy) 
            
            endif
        
        else if (iy == ny) then 
            if (ix==1) then 
                ! pflxt = pflxt &
                    ! + poro(ix,iy)*rho(ix,iy)*mw*(fpx(ix,iy)-fp(ix,iy))/dt   &
                
                pflxadv(ix,iy) = pflxadv(ix,iy) & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*(fpx(ix+1,iy)-fpx(ix,iy))/dx(ix)  &
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(fpx(ix,iy)-fpx(ix,iy-1))/dy(iy)  
                
                pflxdif(ix,iy) = pflxdif(ix,iy) &
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(fpx(ix+1,iy)-1d0*fpx(ix,iy))/(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix)  & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(fpx(ix,iy-1)-1d0*fpx(ix,iy))/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix,iy-1)*rho(ix,iy-1)*dif(ix,iy-1))  &
                        *(fpx(ix,iy)-fpx(ix,iy-1))/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)  
                
                pflxrxn(ix,iy) = pflxrxn(ix,iy) & 
                    - rhob(ix,iy)*ms*mw*kex(ix,iy)*rxn(ix,iy) 
                
            else if (ix==nx) then 
                ! pflxt = pflxt &
                    ! + poro(ix,iy)*rho(ix,iy)*mw*(fpx(ix,iy)-fp(ix,iy))/dt   &
                
                pflxadv(ix,iy) = pflxadv(ix,iy) & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(fpx(ix,iy)-fpx(ix-1,iy))/dx(ix)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(fpx(ix,iy)-fpx(ix,iy-1))/dy(iy)  
                
                pflxdif(ix,iy) = pflxdif(ix,iy) &
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(fpx(ix-1,iy)-1d0*fpx(ix,iy))/(0.5d0*(dx(ix)+dx(Ix-1)))/dx(ix)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix-1,iy)*rho(ix-1,iy)*dif(ix-1,iy)) &
                        *(fpx(ix,iy)-fpx(ix-1,iy))/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix)  &
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(fpx(ix,iy-1)-1d0*fpx(ix,iy))/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix,iy-1)*rho(ix,iy-1)*dif(ix,iy-1))  &
                        *(fpx(ix,iy)-fpx(ix,iy-1))/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)  
                
                pflxrxn(ix,iy) = pflxrxn(ix,iy) &
                    - rhob(ix,iy)*ms*mw*kex(ix,iy)*rxn(ix,iy) 
                
            else 
                ! pflxt = pflxt &
                    ! + poro(ix,iy)*rho(ix,iy)*mw*(fpx(ix,iy)-fp(ix,iy))/dt   
                
                pflxadv(ix,iy) = pflxadv(ix,iy) & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(fpx(ix,iy)-fpx(ix-1,iy))/dx(ix)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*(fpx(ix+1,iy)-fpx(ix,iy))/dx(ix)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(fpx(ix,iy)-fpx(ix,iy-1))/dy(iy)  
                
                pflxdif(ix,iy) = pflxdif(ix,iy) & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*((fpx(ix+1,iy)-fpx(ix,iy))/(0.5d0*(dx(ix)+dx(ix+1)))  &
                        -(fpx(ix,iy)-fpx(ix-1,iy))/(0.5d0*(dx(ix)+dx(Ix-1))))/dx(ix)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix-1,iy)*rho(ix-1,iy)*dif(ix-1,iy))  &
                        *(fpx(ix,iy)-fpx(ix-1,iy))/(0.5d0*(dx(ix)+dx(Ix-1)))/dx(ix)  &
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(fpx(ix,iy-1)-1d0*fpx(ix,iy))/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix,iy-1)*rho(ix,iy-1)*dif(ix,iy-1))  &
                        *(fpx(ix,iy)-fpx(ix,iy-1))/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)  
                
                pflxrxn(ix,iy) = pflxrxn(ix,iy) & 
                    - rhob(ix,iy)*ms*mw*kex(ix,iy)*rxn(ix,iy) 
            
            endif
        
        else 
            if (ix==1) then 
                ! pflxt = pflxt &
                    ! + poro(ix,iy)*rho(ix,iy)*mw*(fpx(ix,iy)-fp(ix,iy))/dt   
                
                pflxadv(ix,iy) = pflxadv(ix,iy) & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*(fpx(ix+1,iy)-fpx(ix,iy))/dx(ix)  &
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(fpx(ix,iy)-fpx(ix,iy-1))/dy(iy)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*(fpx(ix,iy+1)-fpx(ix,iy))/dy(iy)  
                
                pflxdif(ix,iy) = pflxdif(ix,iy) & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(fpx(ix+1,iy)-1d0*fpx(ix,iy))/(0.5d0*(dx(ix)+dx(ix+1)))/dx(ix)  & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*((fpx(ix,iy+1)-fpx(ix,iy))/(0.5d0*(dy(iy)+dy(iy+1)))  &
                        -(fpx(ix,iy)-fpx(ix,iy-1))/(0.5d0*(dy(iy)+dy(iy-1))))/dy(iy)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix,iy-1)*rho(ix,iy-1)*dif(ix,iy-1))  &
                        *(fpx(ix,iy)-fpx(ix,iy-1))/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)  
                
                pflxrxn(ix,iy) = pflxrxn(ix,iy) & 
                    - rhob(ix,iy)*ms*mw*kex(ix,iy)*rxn(ix,iy) 
                
            else if (ix==nx) then 
                ! pflxt = pflxt &
                    ! + poro(ix,iy)*rho(ix,iy)*mw*(fpx(ix,iy)-fp(ix,iy))/dt   
                
                pflxadv(ix,iy) = pflxadv(ix,iy) & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(fpx(ix,iy)-fpx(ix-1,iy))/dx(ix)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(fpx(ix,iy)-fpx(ix,iy-1))/dy(iy)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*(fpx(ix,iy+1)-fpx(ix,iy))/dy(iy)  
                pflxdif(ix,iy) = pflxdif(ix,iy) & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*(fpx(ix-1,iy)-1d0*fpx(ix,iy))/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix-1,iy)*rho(ix-1,iy)*dif(ix-1,iy))  &
                        *(fpx(ix,iy)-fpx(ix-1,iy))/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix)  &
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*((fpx(ix,iy+1)-fpx(ix,iy))/(0.5d0*(dy(iy)+dy(iy+1)))  &
                        -(fpx(ix,iy)-fpx(ix,iy-1))/(0.5d0*(dy(iy)+dy(iy-1))))/dy(iy)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix,iy-1)*rho(ix,iy-1)*dif(ix,iy-1))  &
                        *(fpx(ix,iy)-fpx(ix,iy-1))/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)  
                pflxrxn(ix,iy) = pflxrxn(ix,iy) &
                    - rhob(ix,iy)*ms*mw*kex(ix,iy)*rxn(ix,iy) 
                
            else 
                ! pflxt = pflxt &
                    ! + poro(ix,iy)*rho(ix,iy)*mw*(fpx(ix,iy)-fp(ix,iy))/dt   
                    
                pflxadv(ix,iy) = pflxadv(ix,iy) & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)+abs(vmx(ix,iy)))*0.5d0*(fpx(ix,iy)-fpx(ix-1,iy))/dx(ix)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmx(ix,iy)-abs(vmx(ix,iy)))*0.5d0*(fpx(ix+1,iy)-fpx(ix,iy))/dx(ix)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)+abs(vmy(ix,iy)))*0.5d0*(fpx(ix,iy)-fpx(ix,iy-1))/dy(iy)  & 
                    + poro(ix,iy)*rho(ix,iy)*mw*(vmy(ix,iy)-abs(vmy(ix,iy)))*0.5d0*(fpx(ix,iy+1)-fpx(ix,iy))/dy(iy)  
                
                pflxdif(ix,iy) = pflxdif(ix,iy) & 
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*((fpx(ix+1,iy)-fpx(ix,iy))/(0.5d0*(dx(ix)+dx(ix+1)))  &
                        -(fpx(ix,iy)-fpx(ix-1,iy))/(0.5d0*(dx(ix)+dx(ix-1))))/dx(ix)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix-1,iy)*rho(ix-1,iy)*dif(ix-1,iy))  &
                        *(fpx(ix,iy)-fpx(ix-1,iy))/(0.5d0*(dx(ix)+dx(ix-1)))/dx(ix)  &
                    - poro(ix,iy)*rho(ix,iy)*mw*dif(ix,iy)*((fpx(ix,iy+1)-fpx(ix,iy))/(0.5d0*(dy(iy)+dy(iy+1)))  &
                        -(fpx(ix,iy)-fpx(ix,iy-1))/(0.5d0*(dy(iy)+dy(iy-1))))/dy(iy)  & 
                    - mw*(poro(ix,iy)*rho(ix,iy)*dif(ix,iy)-poro(ix,iy-1)*rho(ix,iy-1)*dif(ix,iy-1))  &
                        *(fpx(ix,iy)-fpx(ix,iy-1))/(0.5d0*(dy(iy)+dy(iy-1)))/dy(iy)  
                
                pflxrxn(ix,iy) = pflxrxn(ix,iy) & 
                    - rhob(ix,iy)*ms*mw*kex(ix,iy)*rxn(ix,iy) 
            
            endif 
            
        endif
        
    enddo
enddo

! rflxadv=rflxadv*dx*dy*rl
! rflxrxn=rflxrxn*dx*dy*rl
! rflxt=rflxt*dx*dy*rl
! pflxadv=pflxadv*dx*dy*rl
! pflxrxn=pflxrxn*dx*dy*rl
! pflxdif=pflxdif*dx*dy*rl
! pflxt=pflxt*dx*dy*rl

! rflxrxnpls = rflxrxnpls*dx*dy*rl
! rflxrxnmns = rflxrxnmns*dx*dy*rl
rflxadv=rflxadv*ds*rl
rflxrxn=rflxrxn*ds*rl
rflxt=rflxt*ds*rl
pflxadv=pflxadv*ds*rl
pflxrxn=pflxrxn*ds*rl
pflxdif=pflxdif*ds*rl
pflxt=pflxt*ds*rl

rflxrxnpls = rflxrxnpls*ds*rl
rflxrxnmns = rflxrxnmns*ds*rl

open(unit=100,file=trim(adjustl(workdir))//'sflxsum.txt',action='write',status='unknown',position='append')
open(unit=200,file=trim(adjustl(workdir))//'pflxsum.txt',action='write',status='unknown',position='append')
print*,'rflxadv, rflxrxn, rflxt, res, rflxrxnpls, rflxrxnmns'
print'(6E9.2)',sum(rflxadv), sum(rflxrxn), sum(rflxt),  sum(rflxadv)+sum(rflxrxn)+sum(rflxt), sum(rflxrxnpls), sum(rflxrxnmns)
write(100,*) o18swi, sum(rflxadv), sum(rflxrxn), sum(rflxt),  sum(rflxadv)+sum(rflxrxn)+sum(rflxt), sum(rflxrxnpls), sum(rflxrxnmns)
print*,'pflxadv, pflxrxn, pflxdif, pflxt, res '
print'(5E9.2)', sum(pflxadv), sum(pflxrxn), sum(pflxdif), sum(pflxt), sum(pflxadv)+sum(pflxrxn)+sum(pflxdif)+sum(pflxt)
write(200,*) o18swi, sum(pflxadv), sum(pflxrxn), sum(pflxdif), sum(pflxt), sum(pflxadv)+sum(pflxrxn)+sum(pflxdif)+sum(pflxt)
close(100)
close(200)

enddo

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



!**************************************************************************************************************************************
function d2r(d)
implicit none
real(kind=8), parameter :: rsmow = 0.0020052d0
real(kind=8) d, d2r

d2r = (1d0+d*1d-3)*rsmow 

end function d2r
!**************************************************************************************************************************************
!**************************************************************************************************************************************
function d2f(d)
implicit none
real(kind=8), parameter :: rsmow = 0.0020052d0
real(kind=8) d, d2r, d2f

d2f = d2r(d)/(1d0+d2r(d))

end function d2f
!**************************************************************************************************************************************
!**************************************************************************************************************************************
function r2d(r)
implicit none
real(kind=8), parameter :: rsmow = 0.0020052d0
real(kind=8) r, r2d

r2d = 1d3*(r-rsmow)/rsmow

end function r2d
!**************************************************************************************************************************************
!**************************************************************************************************************************************
function f2d(f)
implicit none
real(kind=8), parameter :: rsmow = 0.0020052d0
real(kind=8) f, r2d, r, f2d

r = f/(1d0-f)
f2d = r2d(r)

end function f2d
!**************************************************************************************************************************************