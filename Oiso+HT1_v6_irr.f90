program Oiso_HT_main
implicit none

integer(kind=4),parameter :: nx = 200, ny = 320
real(kind=8),parameter :: xmax = 30d3   !  m
real(kind=8),parameter :: ymax = 12d3    !  m 
real(kind=8),parameter :: rg = 8.3144d-3 !  kJ mol^-1 K^-1
real(kind=8),parameter :: E = 50d0 ! kJ
real(kind=8),parameter :: kref_0 = 10d0**(-8.5d0) ! mol-1 kg yr-1
real(kind=8),parameter :: tempref = 5d0 ! C
real(kind=8),parameter :: tempk_0 = 273.15d0  !  K
real(kind=8),parameter :: yr2sec = 60d0*60d0*24d0*365.25d0 ! sec/yr
real(kind=8),parameter :: ms = 0.5d0/16.0d0*1.0d3 ! mol kg-1
real(kind=8),parameter :: mw = 16d0/18d0/16.0d0*1.0d3  !  mol kg-1
real(kind=8) :: w = 3d-2 ! m/yr spreading rate 
real(kind=8) :: lambda = 0.5305d0
real(kind=8) :: gamma = 0d0
real(kind=8) :: poroi = 0.05d0 !  
real(kind=8) :: disp = 10d0 !  dispersivity [m]
real(kind=8) :: beta = 0.876002291d0 ! anddesite (Zhao and Zheng, 2003)
real(kind=8) :: age_lim = 1d3 ! ages (yr) until which rate is the same as that in the laboratory; inferred from Maher et al. (2014)
real(kind=8) :: rhom = 3d3 ! rock density [kg m-3]
real(kind=8),dimension(2) :: fsw, fri,dri, dswi, rsmow
real(kind=8),dimension(2,nx,ny) :: rflxadv, rflxrxn, rflxt, pflxadv, pflxrxn, pflxdif, pflxt  &
    & ,rflxrxnpls, rflxrxnmns, frx, fpx, dr, dp, omega,kex,alfa,dif
real(kind=8),dimension(nx,ny) :: capDp17,capDr17 &
    & ,ds,temp,rho,vmx,vmy,visc,vabs,poro,rhob,kref,theta_eq,theta_kin,wr
real(kind=8),dimension(nx) :: dx, xm
real(kind=8),dimension(ny) :: dy, ym
real(kind=8),dimension(nx+1) :: xe 
real(kind=8),dimension(ny+1) :: ye 
logical:: kref_variable = .false.  ! default
! logical:: kref_variable = .true.  ! rate reciprocal of age is imposed after Maher et al. (2004)
character*500 workdir,basefield
integer isw,iiso,ix,iy
character*2 intsw, isochr
character*1 signsw
real(kind=8) beta_grid 

real(kind=8) dp2d,mwl_Luz,d2r,r2f,f2r,r2d,r2dp
!-------------------------------------------------------------------------------------

rsmow(1) = 2.0052d-3
rsmow(2) = 3.799d-4 
dri(1) = 5.7d0
dri(2) = 2.86d0
fri(1) = r2f(d2r(dri(1),rsmow(1)),d2r(dri(2),rsmow(2)))
fri(2) = r2f(d2r(dri(2),rsmow(2)),d2r(dri(1),rsmow(1)))

basefield = 'perm_expexp_-16_8-11_8_zh300_spx1_200x320_irr-20200830'

workdir = '../ht-oiso_output/'//trim(adjustl(basefield))//'/'
 
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

! modified Stokes-Einstein (Krynicki et al., 1978)
dif(1,:,:) = 6.9d-15*(tempk_0+temp(:,:))/visc(:,:)  ! water diffusion [m2 s-1]
dif(2,:,:) = 6.9d-15*(tempk_0+temp(:,:))/visc(:,:)  ! water diffusion [m2 s-1]
dif(1,:,:) = dif(1,:,:)/sqrt(20d0/18d0)  ! correction for H218O (Harris et al., 1979)
dif(2,:,:) = dif(2,:,:)/sqrt(19d0/18d0)  ! correction for H217O (just guessing)
dif = dif*yr2sec   ! converting sec to yr 
dif(1,:,:) = dif(1,:,:)*poro(:,:)**1.4d0 + disp*vabs(:,:)  !  assuming homogeneous dispersivity 
dif(2,:,:) = dif(2,:,:)*poro(:,:)**1.4d0 + disp*vabs(:,:)  !  assuming homogeneous dispersivity 

rhob = (1d0-poro)*rhom + poro*rho

wr = poro*rho*mw*vabs/((1d0-poro)*rhom*ms*w)

alfa(1,:,:) = (6.673d6/(tempk_0+temp(:,:))**2.0d0+10.398d3/(tempk_0+temp(:,:))-4.78d0) &
    & *exp((1.0d0-beta)/(rg*(tempk_0+temp(:,:))))*beta  &
    & -(2.194d6/(tempk_0+temp(:,:))**2.0d0+15.163d3/(tempk_0+temp(:,:))-4.72d0)+1.767d0*(2.0d0*beta-1.0d0)
alfa(1,:,:) = exp(alfa(1,:,:)/1d3) ! changing from o/oo to fractionation factor

theta_eq = -1.85d0/(tempk_0+temp) + 0.5305d0 ! sharp16
! theta_eq = -740d0/(tempk_0+temp)/(tempk_0+temp) + 0.5305d0 ! pack14

alfa(2,:,:) = alfa(1,:,:)**theta_eq(:,:)

kref = kref_0  ! constant kref
!  asuming change rate with age 
if (kref_variable)then 
    do ix=1,nx
        if (xm(ix)/w <= age_lim) then 
            kref(ix,:) = 10d0**(-6.9d0)
        elseif (xm(ix)/w > age_lim) then 
            kref(ix,:) = 10d0**(-6.9d0)/(xm(ix)/w/age_lim)
        endif 
    enddo
endif 
kex(1,:,:) = kref(:,:)*exp(-E*(1.0d0/(tempk_0+temp(:,:))-1.0d0/(tempk_0+tempref))/rg)

theta_kin = 1d0
kex(2,:,:) = kex(1,:,:)**theta_kin(:,:)

open(unit=100,file=trim(adjustl(workdir))//'k1.txt',action='write',status='replace')
open(unit=200,file=trim(adjustl(workdir))//'k2.txt',action='write',status='replace')
open(unit=300,file=trim(adjustl(workdir))//'alfa1.txt',action='write',status='replace')
open(unit=400,file=trim(adjustl(workdir))//'alfa2.txt',action='write',status='replace')
open(unit=500,file=trim(adjustl(workdir))//'dif1.txt',action='write',status='replace')
open(unit=600,file=trim(adjustl(workdir))//'dif2.txt',action='write',status='replace')
open(unit=700,file=trim(adjustl(workdir))//'vabs.txt',action='write',status='replace')
open(unit=800,file=trim(adjustl(workdir))//'wr.txt',action='write',status='replace')
do iy=1,ny
    write(100,*) (log10(kex(1,ix,iy)),ix=1,nx) 
    write(200,*) (log10(kex(2,ix,iy)),ix=1,nx) 
    write(300,*) (alfa(1,ix,iy),ix=1,nx) 
    write(400,*) (alfa(2,ix,iy),ix=1,nx) 
    write(500,*) (dif(1,ix,iy),ix=1,nx) 
    write(600,*) (dif(2,ix,iy),ix=1,nx) 
    write(700,*) (vabs(ix,iy),ix=1,nx) 
    write(800,*) (wr(ix,iy),ix=1,nx)
enddo
close(100)
close(200)
close(300)
close(400)
close(500)
close(600)
close(700)
close(800)

do isw = -5,13, 2
! do isw = 1,13, 2
! do isw = 1,1

    dswi(1) = -1d0*(isw-1)
    write(intsw,'(I2.2)') int(abs(dswi(1)))
    if (dswi(1) > 0d0) then 
        signsw = 'p'
    elseif (dswi(1) == 0d0) then
        signsw = '0'
    elseif (dswi(1) < 0d0) then 
        signsw = 'n'
    endif 

    dswi(2) = dp2d(mwl_Luz(dswi(1),rsmow(1)),rsmow(2))

    fsw(1) = r2f(d2r(dswi(1),rsmow(1)),d2r(dswi(2),rsmow(2)))
    fsw(2) = r2f(d2r(dswi(2),rsmow(2)),d2r(dswi(1),rsmow(1)))

    do iiso = 1,2
        
        call Oiso_HT1( &
            & nx,ny,xmax,ymax,w,fsw(iiso),fri(iiso),basefield &! input 
            & ,dx,dy,ds,rho,vmx,vmy,poro,rhob,kex(iiso,:,:),alfa(iiso,:,:),dif(iiso,:,:)  &! input
            & ,frx(iiso,:,:),fpx(iiso,:,:),omega(iiso,:,:) &! output
            & ,rflxadv(iiso,:,:),rflxrxn(iiso,:,:),rflxt(iiso,:,:),rflxrxnpls(iiso,:,:),rflxrxnmns(iiso,:,:) &! output
            & ,pflxadv(iiso,:,:),pflxrxn(iiso,:,:),pflxdif(iiso,:,:),pflxt(iiso,:,:) &! output
            )
    enddo 

    do iy=1,ny
        do ix=1,nx
            dr(1,ix,iy) = r2d(f2r(frx(1,ix,iy),frx(2,ix,iy)),rsmow(1))
            dp(1,ix,iy) = r2d(f2r(fpx(1,ix,iy),fpx(2,ix,iy)),rsmow(1))
            dr(2,ix,iy) = r2d(f2r(frx(2,ix,iy),frx(1,ix,iy)),rsmow(2))
            dp(2,ix,iy) = r2d(f2r(fpx(2,ix,iy),fpx(1,ix,iy)),rsmow(2))
            capDr17(ix,iy) = r2dp(f2r(frx(2,ix,iy),frx(1,ix,iy)),rsmow(2))  &
                & - (lambda*r2dp(f2r(frx(1,ix,iy),frx(2,ix,iy)),rsmow(1)) + gamma)
            capDp17(ix,iy) = r2dp(f2r(fpx(2,ix,iy),fpx(1,ix,iy)),rsmow(2))  &
                & - (lambda*r2dp(f2r(fpx(1,ix,iy),fpx(2,ix,iy)),rsmow(1)) + gamma) 
        enddo
    enddo
    
    do iiso = 1,2
        isochr = '18'
        if (iiso ==2) isochr = '17'
        
        open(unit=100,file=trim(adjustl(workdir))//'d'//isochr//'r-'//trim(adjustl(signsw))//trim(adjustl(intsw))//'.txt'  &
            ,action='write',status='replace')
        open(unit=200,file=trim(adjustl(workdir))//'d'//isochr//'p-'//trim(adjustl(signsw))//trim(adjustl(intsw))//'.txt'  &
            ,action='write',status='replace')
        open(unit=300,file=trim(adjustl(workdir))//'omage'//isochr//'-'//trim(adjustl(signsw))//trim(adjustl(intsw))//'.txt'  &
            ,action='write',status='replace')
        do iy=1,ny
            write(100,*) (dr(iiso,ix,iy),ix=1,nx) 
            write(200,*) (dp(iiso,ix,iy),ix=1,nx) 
            write(300,*) (log10(omega(iiso,ix,iy)),ix=1,nx) 
        enddo
        close(100)
        close(200)
        close(300)

        open(unit=100,file=trim(adjustl(workdir))//'sflxsum'//isochr//'.txt',action='write',status='unknown',position='append')
        open(unit=200,file=trim(adjustl(workdir))//'pflxsum'//isochr//'.txt',action='write',status='unknown',position='append')
        print*,'rflxadv, rflxrxn, rflxt, res, rflxrxnpls, rflxrxnmns'
        print'(6E9.2)',sum(rflxadv(iiso,:,:)), sum(rflxrxn(iiso,:,:)), sum(rflxt(iiso,:,:)) &
            & ,sum(rflxadv(iiso,:,:))+sum(rflxrxn(iiso,:,:))+sum(rflxt(iiso,:,:)) &
            & ,sum(rflxrxnpls(iiso,:,:)), sum(rflxrxnmns(iiso,:,:))
        write(100,*) dswi(iiso), sum(rflxadv(iiso,:,:)), sum(rflxrxn(iiso,:,:)), sum(rflxt(iiso,:,:))  &
            & ,sum(rflxadv(iiso,:,:))+sum(rflxrxn(iiso,:,:))+sum(rflxt(iiso,:,:)) &
            & ,sum(rflxrxnpls(iiso,:,:)), sum(rflxrxnmns(iiso,:,:))
        print*,'pflxadv, pflxrxn, pflxdif, pflxt, res '
        print'(5E9.2)', sum(pflxadv(iiso,:,:)), sum(pflxrxn(iiso,:,:)), sum(pflxdif(iiso,:,:)), sum(pflxt(iiso,:,:))  &
            & ,sum(pflxadv(iiso,:,:))+sum(pflxrxn(iiso,:,:))+sum(pflxdif(iiso,:,:))+sum(pflxt(iiso,:,:))
        write(200,*) dswi(iiso), sum(pflxadv(iiso,:,:)), sum(pflxrxn(iiso,:,:)), sum(pflxdif(iiso,:,:)), sum(pflxt(iiso,:,:)) &
            & ,sum(pflxadv(iiso,:,:))+sum(pflxrxn(iiso,:,:))+sum(pflxdif(iiso,:,:))+sum(pflxt(iiso,:,:))
        close(100)
        close(200)
    
    enddo 

        
    open(unit=100,file=trim(adjustl(workdir))//'capd17r-'//trim(adjustl(signsw))//trim(adjustl(intsw))//'.txt'  &
        ,action='write',status='replace')
    open(unit=200,file=trim(adjustl(workdir))//'capd17p-'//trim(adjustl(signsw))//trim(adjustl(intsw))//'.txt'  &
        ,action='write',status='replace')
    do iy=1,ny
        write(100,*) (capDr17(ix,iy),ix=1,nx) 
        write(200,*) (capDp17(ix,iy),ix=1,nx) 
    enddo
    close(100)
    close(200)
    
enddo

EndProgram

!**************************************************************************************************************************************
!**************************************************************************************************************************************
!**************************************************************************************************************************************

!**************************************************************************************************************************************
subroutine Oiso_HT1( &
    & nx,ny,xmax,ymax,w,fsw,fri,basefield  &! input 
    & ,dx,dy,ds,rho,vmx,vmy,poro,rhob,kex,alfa,dif  &! input
    & ,frx,fpx,omega &! output
    & ,rflxadv,rflxrxn,rflxt,rflxrxnpls,rflxrxnmns &! output
    & ,pflxadv,pflxrxn,pflxdif,pflxt &! output
    )
! oxygen isotope calculation using 2D data for flow, temperature and density
! looking working 
! modified from v3 to include the change in rho for diffusion term
! modified from v4 to include the temp dependence of diffusion coefficient
! remove non-sparse matrix solver to save the memory
! from v5 to remove porosity before q
! from v6 trying irregular grid 
! enabling application to 17O
implicit none

integer(kind=4),intent(in) :: nx,ny
real(kind=8),intent(in) :: xmax,ymax,w,fsw,fri
real(kind=8),dimension(nx),intent(in) :: dx
real(kind=8),dimension(ny),intent(in) :: dy
real(kind=8),dimension(nx,ny),intent(in) :: ds,rho,vmx,vmy,poro,rhob,kex,alfa,dif
character(500),intent(in) :: basefield

real(kind=8),dimension(nx,ny),intent(out)::frx,fpx,omega  &
    & ,rflxadv,rflxrxn,rflxt,rflxrxnpls,rflxrxnmns &
    & ,pflxadv,pflxrxn,pflxdif,pflxt
    
real(kind=8) :: rl = 1d8 ! m ridge length
real(kind=8) dt, time  ! m,m, yr, yr
real(kind=8) rxn(nx,ny) 
real(kind=8) drxn_dfr(nx,ny),drxn_dfp(nx,ny) 
real(kind=8) :: rhom = 3d3 ! rock density [kg m-3]
integer(kind = 4) ix,iy
! o isotope 
real(kind=8) fr(nx,ny), fp(nx,ny)
real(kind=8), parameter :: ms = 0.5d0/16.0d0*1.0d3 ! mol kg-1
real(kind=8), parameter :: mw = 16d0/18d0/16.0d0*1.0d3  !  mol kg-1
! matrix solver
integer(kind = 4) nmx 
real(kind=8), allocatable :: emx(:)
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
real(kind=8),parameter :: tol = 1d-6
integer(kind=4) cnt, cnt2, ixp, ixn, iyp, iyn, tmpint(6),tmpint2(6)
real(kind=8) tmprl(6)  
!-------------------------------------------------------------------------------------

fr=fri
fp=fsw

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


omega = frx*(1d0-fpx)/(1d0-frx)/fpx/alfa

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

Endsubroutine Oiso_HT1
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



! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function r2d(ratio,rstd)
implicit none
real(kind=8) r2d,ratio,rstd
r2d=(ratio/rstd-1d0)*1d3
endfunction r2d

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function d2r(delta,rstd)
implicit none
real(kind=8) d2r,delta,rstd
d2r=(delta*1d-3+1d0)*rstd
endfunction d2r

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function r2dp(ratio,rstd)
implicit none
real(kind=8) r2dp,ratio,rstd
r2dp=1d3*log(ratio/rstd)
endfunction r2dp

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function dp2r(dp,rstd)
implicit none
real(kind=8) dp2r,dp,rstd
dp2r=rstd*exp(dp*1d-3)
endfunction dp2r

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function dp2d(dp,rstd)
implicit none
real(kind=8) dp2d,dp,rstd,dp2r, r2d
dp2d = r2d(dp2r(dp,rstd),rstd)  
endfunction dp2d

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function d2dp(d,rstd)
implicit none
real(kind=8) d2dp,d,rstd,d2r, r2dp
d2dp = r2dp(d2r(d,rstd),rstd)  
endfunction d2dp

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function r2f(r1,r2)
implicit none
real(kind=8) r2f,r1,r2
r2f = r1/(1d0+r1+r2)
endfunction r2f

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function f2r(f1,f2)
implicit none
real(kind=8) f2r,f1,f2
f2r = f1/(1d0-f1-f2)
endfunction f2r

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function mwl_Luz(d18o,rstd)  ! returning d'17 based on d18O and meteoric water line by Luz and Barkan (2010) 
implicit none
real(kind=8) mwl_Luz,d18o,d2r,r2dp,rstd
mwl_Luz = 0.528d0*r2dp(d2r(d18o,rstd),rstd) + 0.033d0
endfunction mwl_Luz

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
