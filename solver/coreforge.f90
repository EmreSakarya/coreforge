!=======================================================================
!  COREFORGE v8.6
!  2-D / 3-D multigroup neutron diffusion eigenvalue engine (k_eff)
!  for arbitrary reactor cores.  NZ = 1 (default) reproduces the 2-D
!  engine bit-for-bit; NZ > 1 solves full x-y-z cores with per-layer
!  material maps (partially inserted control rods, axial reflectors).
!
!  Physics
!    -div( D_g grad phi_g ) + Sr_g phi_g
!        = chi_g/k * SUM_g' nuSf_g' phi_g'  +  SUM_g'/=g  Ss(g'->g) phi_g'
!    Sr_g = Sa_g + SUM_g'/=g Ss(g->g')          (removal)
!
!  Numerics
!    * mesh-centred finite differences, harmonic-mean interface D
!    * per-face BCs on all 6 faces: reflective, or Robin vacuum
!         J_out = gamma * phi_face   (0.4692 transport-corrected default)
!    * power iteration; per-group red-black SOR inner sweeps,
!      OpenMP-parallel over (z,y) rows; full NG x NG scattering matrix
!
!  Input : keyword text file. 3-D additions (all optional):
!            NZ <layers>   DZ <cm>
!            BC W E S N [Bottom Top]      (0 vacuum, 1 reflective)
!            MAP [l1 l2]                  (layer range, bottom-based;
!                                          bare MAP fills all layers)
!  Output: stdout summary; flux.csv / power.csv with x,y,z columns.
!
!  Build : ifx      -O3 -qopenmp coreforge.f90 -o coreforge
!          gfortran -O3 -fopenmp coreforge.f90 -o coreforge
!=======================================================================
module cf
  use iso_fortran_env, only: error_unit
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  ! ---- problem definition -------------------------------------------
  integer  :: ng = 0, nx = 0, ny = 0, nz = 1, nmat = 0
  real(dp) :: dx = 0.0_dp, dy = 0.0_dp, dz = 0.0_dp
  integer  :: bc(6) = (/0,0,0,0,1,1/)   ! W,E,S,N,Bottom,Top
  real(dp) :: gam    = 0.4692_dp        ! vacuum J/phi ratio
  real(dp) :: omega  = 1.6_dp           ! SOR over-relaxation factor
  real(dp) :: tolk   = 1.0e-7_dp
  real(dp) :: tols   = 1.0e-5_dp
  real(dp) :: intol  = 1.0e-4_dp        ! inner within-group convergence
  integer  :: maxout = 8000, ninner = 2 ! ninner = minimum inner sweeps
  integer  :: ninmax = 400              ! cap on adaptive inner sweeps

  ! ---- cross sections (g,mat) / (from,to,mat) ------------------------
  real(dp), allocatable :: xD(:,:), xSa(:,:), xNf(:,:), xChi(:,:), xSs(:,:,:)
  integer,  allocatable :: matmap(:,:,:)        ! (i,j,k) k=1 bottom
  logical,  allocatable :: matfis(:)
  logical,  allocatable :: fuelmask(:,:,:)

  ! ---- fields ---------------------------------------------------------
  real(dp), allocatable :: phi(:,:,:,:)         ! (i,j,k,g)
  real(dp), allocatable :: fis(:,:,:), fisold(:,:,:), src(:,:,:)

  ! ---- precomputed 7-point coefficients (i,j,k,g) ---------------------
  real(dp), allocatable :: dg(:,:,:,:)
  real(dp), allocatable :: cw(:,:,:,:), ce(:,:,:,:)
  real(dp), allocatable :: cs(:,:,:,:), cn(:,:,:,:)
  real(dp), allocatable :: cb(:,:,:,:), ct(:,:,:,:)

contains

  subroutine fatal(msg)
    character(len=*), intent(in) :: msg
    write(error_unit,'(a)') 'ERROR: '//trim(msg)
    error stop 1
  end subroutine fatal

  pure function upcase(s) result(t)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: t
    integer :: k, c
    t = s
    do k = 1, len_trim(s)
       c = iachar(s(k:k))
       if (c >= iachar('a') .and. c <= iachar('z')) t(k:k) = achar(c-32)
    end do
  end function upcase

  subroutine next_line(u, line, ok, lineno)
    integer, intent(in)    :: u
    character(len=*), intent(out) :: line
    logical, intent(out)   :: ok
    integer, intent(inout) :: lineno
    integer :: ios, p
    ok = .false.
    do
       read(u,'(A)',iostat=ios) line
       if (ios /= 0) return
       lineno = lineno + 1
       p = index(line, char(13)); if (p > 0) line(p:) = ' '
       p = index(line, '#');      if (p > 0) line(p:) = ' '
       if (len_trim(line) > 0) then
          ok = .true.
          return
       end if
    end do
  end subroutine next_line

  !--------------------------------------------------------------------
  subroutine read_input(fname)
    character(len=*), intent(in) :: fname
    character(len=65536) :: line
    character(len=32)    :: key, key2
    integer :: u, ios, lineno, mid, gf, gt, i, r, m, l1, l2, kk
    logical :: ok
    real(dp) :: csum
    integer, allocatable :: rowbuf(:)

    open(newunit=u, file=fname, status='old', action='read', iostat=ios)
    if (ios /= 0) call fatal('cannot open input file: '//trim(fname))
    lineno = 0

    do
       call next_line(u, line, ok, lineno); if (.not. ok) exit
       read(line,*,iostat=ios) key
       if (ios /= 0) call fatal('unreadable line in input')
       key = upcase(key)

       select case (trim(key))
       case ('NG');     read(line,*,iostat=ios) key, ng
       case ('NX');     read(line,*,iostat=ios) key, nx
       case ('NY');     read(line,*,iostat=ios) key, ny
       case ('NZ');     read(line,*,iostat=ios) key, nz
       case ('DX');     read(line,*,iostat=ios) key, dx
       case ('DY');     read(line,*,iostat=ios) key, dy
       case ('DZ');     read(line,*,iostat=ios) key, dz
       case ('BC')
          read(line,*,iostat=ios) key, bc(1:6)
          if (ios /= 0) then                 ! legacy 4-value form
             read(line,*,iostat=ios) key, bc(1:4)
             bc(5:6) = 1
          end if
       case ('GAMMA');  read(line,*,iostat=ios) key, gam
       case ('OMEGA');  read(line,*,iostat=ios) key, omega
       case ('TOLK');   read(line,*,iostat=ios) key, tolk
       case ('TOLS');   read(line,*,iostat=ios) key, tols
       case ('NINNER'); read(line,*,iostat=ios) key, ninner
       case ('NINMAX'); read(line,*,iostat=ios) key, ninmax
       case ('INTOL');  read(line,*,iostat=ios) key, intol
       case ('MAXOUT'); read(line,*,iostat=ios) key, maxout

       case ('NMAT')
          read(line,*,iostat=ios) key, nmat
          if (ios /= 0 .or. nmat < 1) call fatal('bad NMAT')
          if (ng < 1) call fatal('NG must be given before NMAT')
          allocate(xD(ng,nmat), xSa(ng,nmat), xNf(ng,nmat), xChi(ng,nmat), &
                   xSs(ng,ng,nmat), matfis(nmat))
          xD=0; xSa=0; xNf=0; xChi=0; xSs=0; matfis=.false.

       case ('MAT')
          if (nmat < 1) call fatal('NMAT must come before MAT blocks')
          read(line,*,iostat=ios) key, mid
          if (ios /= 0 .or. mid < 1 .or. mid > nmat) call fatal('bad MAT id')
          call next_line(u,line,ok,lineno); if (.not.ok) call fatal('EOF in MAT block')
          read(line,*,iostat=ios) key2, (xD(gf,mid),  gf=1,ng)
          if (ios/=0 .or. upcase(trim(key2))/='D')    call fatal('expected D line in MAT block')
          call next_line(u,line,ok,lineno); if (.not.ok) call fatal('EOF in MAT block')
          read(line,*,iostat=ios) key2, (xSa(gf,mid), gf=1,ng)
          if (ios/=0 .or. upcase(trim(key2))/='SA')   call fatal('expected SA line in MAT block')
          call next_line(u,line,ok,lineno); if (.not.ok) call fatal('EOF in MAT block')
          read(line,*,iostat=ios) key2, (xNf(gf,mid), gf=1,ng)
          if (ios/=0 .or. upcase(trim(key2))/='NUSF') call fatal('expected NUSF line in MAT block')
          call next_line(u,line,ok,lineno); if (.not.ok) call fatal('EOF in MAT block')
          read(line,*,iostat=ios) key2, (xChi(gf,mid), gf=1,ng)
          if (ios/=0 .or. upcase(trim(key2))/='CHI')  call fatal('expected CHI line in MAT block')
          call next_line(u,line,ok,lineno); if (.not.ok) call fatal('EOF in MAT block')
          read(line,*,iostat=ios) key2
          if (ios/=0 .or. upcase(trim(key2))/='SCAT') call fatal('expected SCAT line in MAT block')
          do gf = 1, ng
             call next_line(u,line,ok,lineno); if (.not.ok) call fatal('EOF in SCAT matrix')
             read(line,*,iostat=ios) (xSs(gf,gt,mid), gt=1,ng)
             if (ios /= 0) call fatal('bad SCAT row')
          end do
          matfis(mid) = any(xNf(:,mid) > 0.0_dp)

       case ('MAP')
          if (nx<1 .or. ny<1) call fatal('NX/NY must come before MAP')
          if (nmat < 1)       call fatal('NMAT must come before MAP')
          if (.not. allocated(matmap)) then
             allocate(matmap(nx,ny,max(nz,1)))
             matmap = 0
          end if
          l1 = 1; l2 = nz
          read(line,*,iostat=ios) key, l1, l2
          if (ios /= 0) then
             l1 = 1; l2 = nz
             ios = 0
          end if
          if (l1 < 1 .or. l2 > nz .or. l1 > l2) call fatal('bad MAP layer range')
          if (.not. allocated(rowbuf)) allocate(rowbuf(nx))
          do r = 1, ny                      ! file row 1 = TOP (max y)
             call next_line(u,line,ok,lineno); if (.not.ok) call fatal('EOF inside MAP')
             read(line,*,iostat=ios) (rowbuf(i), i=1,nx)
             if (ios /= 0) call fatal('bad MAP row')
             do kk = l1, l2
                matmap(:, ny-r+1, kk) = rowbuf
             end do
          end do

       case default
          call fatal('unknown keyword: '//trim(key))
       end select
       if (ios /= 0) call fatal('bad value after keyword '//trim(key))
    end do
    close(u)

    ! ---- validation ---------------------------------------------------
    if (ng<1 .or. ng>20)          call fatal('NG out of range 1..20')
    if (nx<3 .or. ny<3)           call fatal('mesh must be at least 3x3')
    if (nz<1)                     call fatal('NZ must be >= 1')
    if (dx<=0 .or. dy<=0)         call fatal('DX/DY must be positive')
    if (nz > 1 .and. dz <= 0)     call fatal('DZ must be positive when NZ > 1')
    if (nz == 1 .and. dz <= 0)    dz = 1.0_dp
    if (any(bc<0) .or. any(bc>1)) call fatal('BC values must be 0 or 1')
    if (.not. allocated(matmap))  call fatal('MAP block missing')
    if (minval(matmap)<1 .or. maxval(matmap)>nmat) &
       call fatal('MAP does not cover every layer with valid material ids')
    do m = 1, nmat
       if (any(xD(:,m) <= 0.0_dp)) call fatal('all D values must be positive')
       if (any(xSa(:,m) < 0.0_dp)) call fatal('SA values must be non-negative')
       csum = sum(xChi(:,m))
       if (matfis(m) .and. abs(csum-1.0_dp) > 1.0e-6_dp) &
          write(error_unit,'(a,i0,a,f10.6)') 'WARNING: chi of material ', m, &
               ' does not sum to 1: ', csum
    end do
    if (omega <= 0.0_dp .or. omega >= 2.0_dp) call fatal('OMEGA must be in (0,2)')
  end subroutine read_input

  !--------------------------------------------------------------------
  subroutine setup_coefficients()
    integer :: i, j, k, g, m
    real(dp) :: ax, ay, az, rem, dd, dn, c
    allocate(dg(nx,ny,nz,ng))
    allocate(cw(nx,ny,nz,ng), ce(nx,ny,nz,ng))
    allocate(cs(nx,ny,nz,ng), cn(nx,ny,nz,ng))
    allocate(cb(nx,ny,nz,ng), ct(nx,ny,nz,ng))
    ax = 1.0_dp/(dx*dx);  ay = 1.0_dp/(dy*dy);  az = 1.0_dp/(dz*dz)

    do g = 1, ng
       do k = 1, nz
          do j = 1, ny
             do i = 1, nx
                m  = matmap(i,j,k)
                dd = xD(g,m)
                rem = xSa(g,m) + sum(xSs(g,:,m)) - xSs(g,g,m)
                cw(i,j,k,g)=0; ce(i,j,k,g)=0; cs(i,j,k,g)=0
                cn(i,j,k,g)=0; cb(i,j,k,g)=0; ct(i,j,k,g)=0
                dg(i,j,k,g) = rem
                if (i > 1) then
                   dn = xD(g,matmap(i-1,j,k)); c = 2.0_dp*dd*dn/(dd+dn)*ax
                   cw(i,j,k,g) = c;  dg(i,j,k,g) = dg(i,j,k,g) + c
                else if (bc(1) == 0) then
                   dg(i,j,k,g) = dg(i,j,k,g) + vacterm(dd,dx)
                end if
                if (i < nx) then
                   dn = xD(g,matmap(i+1,j,k)); c = 2.0_dp*dd*dn/(dd+dn)*ax
                   ce(i,j,k,g) = c;  dg(i,j,k,g) = dg(i,j,k,g) + c
                else if (bc(2) == 0) then
                   dg(i,j,k,g) = dg(i,j,k,g) + vacterm(dd,dx)
                end if
                if (j > 1) then
                   dn = xD(g,matmap(i,j-1,k)); c = 2.0_dp*dd*dn/(dd+dn)*ay
                   cs(i,j,k,g) = c;  dg(i,j,k,g) = dg(i,j,k,g) + c
                else if (bc(3) == 0) then
                   dg(i,j,k,g) = dg(i,j,k,g) + vacterm(dd,dy)
                end if
                if (j < ny) then
                   dn = xD(g,matmap(i,j+1,k)); c = 2.0_dp*dd*dn/(dd+dn)*ay
                   cn(i,j,k,g) = c;  dg(i,j,k,g) = dg(i,j,k,g) + c
                else if (bc(4) == 0) then
                   dg(i,j,k,g) = dg(i,j,k,g) + vacterm(dd,dy)
                end if
                if (k > 1) then
                   dn = xD(g,matmap(i,j,k-1)); c = 2.0_dp*dd*dn/(dd+dn)*az
                   cb(i,j,k,g) = c;  dg(i,j,k,g) = dg(i,j,k,g) + c
                else if (bc(5) == 0) then
                   dg(i,j,k,g) = dg(i,j,k,g) + vacterm(dd,dz)
                end if
                if (k < nz) then
                   dn = xD(g,matmap(i,j,k+1)); c = 2.0_dp*dd*dn/(dd+dn)*az
                   ct(i,j,k,g) = c;  dg(i,j,k,g) = dg(i,j,k,g) + c
                else if (bc(6) == 0) then
                   dg(i,j,k,g) = dg(i,j,k,g) + vacterm(dd,dz)
                end if
                if (dg(i,j,k,g) <= 0.0_dp) &
                   call fatal('non-positive diagonal (check cross sections)')
             end do
          end do
       end do
    end do
  end subroutine setup_coefficients

  pure function vacterm(Dv, del) result(t)
    real(dp), intent(in) :: Dv, del
    real(dp) :: t
    t = 2.0_dp*gam*Dv / (del*(gam*del + 2.0_dp*Dv))
  end function vacterm

  !--------------------------------------------------------------------
  subroutine build_fission(fsum)
    real(dp), intent(out) :: fsum
    integer :: i, j, k, g, m
    real(dp) :: f
    fsum = 0.0_dp
    !$omp parallel do collapse(2) private(i,j,k,g,m,f) reduction(+:fsum) schedule(static)
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             m = matmap(i,j,k)
             f = 0.0_dp
             do g = 1, ng
                f = f + xNf(g,m)*phi(i,j,k,g)
             end do
             fis(i,j,k) = f
             fsum = fsum + f
          end do
       end do
    end do
    !$omp end parallel do
  end subroutine build_fission

  subroutine build_source(g, keff)
    integer,  intent(in) :: g
    real(dp), intent(in) :: keff
    integer :: i, j, k, gp, m
    real(dp) :: s, invk
    invk = 1.0_dp/keff
    !$omp parallel do collapse(2) private(i,j,k,gp,m,s) schedule(static)
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             m = matmap(i,j,k)
             s = xChi(g,m)*fis(i,j,k)*invk
             do gp = 1, ng
                if (gp /= g) s = s + xSs(gp,g,m)*phi(i,j,k,gp)
             end do
             src(i,j,k) = s
          end do
       end do
    end do
    !$omp end parallel do
  end subroutine build_source

  !  one red-black SOR sweep for group g (colour = parity of i+j+k);
  !  returns dmax = the largest absolute change in phi over this sweep,
  !  so the caller can iterate the inner solve to convergence.
  subroutine rbsor(g, dmax)
    integer, intent(in)   :: g
    real(dp), intent(out) :: dmax
    integer :: i, j, k, color, istart
    real(dp) :: gs, dloc
    dmax = 0.0_dp
    do color = 0, 1
       !$omp parallel do collapse(2) private(i,j,k,istart,gs,dloc) reduction(max:dmax) schedule(static)
       do k = 1, nz
          do j = 1, ny
             istart = 1 + mod(j+k+color, 2)
             do i = istart, nx, 2
                gs = ( src(i,j,k)                                        &
                     + cw(i,j,k,g)*phi(max(i-1,1 ),j,k,g)                &
                     + ce(i,j,k,g)*phi(min(i+1,nx),j,k,g)                &
                     + cs(i,j,k,g)*phi(i,max(j-1,1 ),k,g)                &
                     + cn(i,j,k,g)*phi(i,min(j+1,ny),k,g)                &
                     + cb(i,j,k,g)*phi(i,j,max(k-1,1 ),g)                &
                     + ct(i,j,k,g)*phi(i,j,min(k+1,nz),g) ) / dg(i,j,k,g)
                gs = phi(i,j,k,g) + omega*(gs - phi(i,j,k,g))
                if (gs < 0.0_dp) gs = 0.0_dp
                dloc = abs(gs - phi(i,j,k,g))
                if (dloc > dmax) dmax = dloc
                phi(i,j,k,g) = gs
             end do
          end do
       end do
       !$omp end parallel do
    end do
  end subroutine rbsor

  !  max relative residual of the within-group linear system for group g:
  !  ||src + (neighbour couplings) - dg*phi||_inf / ||src||_inf.
  !  Unlike the per-sweep change, this is a true convergence measure — it
  !  goes to zero only when the flux actually solves the group equation,
  !  so the inner solver cannot stop while still crawling.  NOTE the norm
  !  is the SOURCE scale (physical); normalising by dg*phi would divide by
  !  ~1/h^2 at fine mesh and make the residual look spuriously small.
  function grp_resid(g) result(rr)
    integer, intent(in) :: g
    real(dp) :: rr, resmax, smax, res
    integer :: i, j, k
    resmax = 0.0_dp;  smax = 0.0_dp
    !$omp parallel do collapse(2) private(i,j,k,res) reduction(max:resmax,smax) schedule(static)
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             res = src(i,j,k)                                        &
                  + cw(i,j,k,g)*phi(max(i-1,1 ),j,k,g)                &
                  + ce(i,j,k,g)*phi(min(i+1,nx),j,k,g)                &
                  + cs(i,j,k,g)*phi(i,max(j-1,1 ),k,g)                &
                  + cn(i,j,k,g)*phi(i,min(j+1,ny),k,g)                &
                  + cb(i,j,k,g)*phi(i,j,max(k-1,1 ),g)                &
                  + ct(i,j,k,g)*phi(i,j,min(k+1,nz),g)                &
                  - dg(i,j,k,g)*phi(i,j,k,g)
             if (abs(res)        > resmax) resmax = abs(res)
             if (abs(src(i,j,k)) > smax  ) smax   = abs(src(i,j,k))
          end do
       end do
    end do
    rr = resmax / max(smax, 1.0e-30_dp)
  end function grp_resid

  !--------------------------------------------------------------------
  subroutine balance_report(keff)
    real(dp), intent(in) :: keff
    real(dp) :: prod(ng), ains(ng), aabs(ng), aout(ng), aleak(ng)
    integer :: i, j, k, g, gp, m
    prod=0; ains=0; aabs=0; aout=0
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             m = matmap(i,j,k)
             do g = 1, ng
                aabs(g) = aabs(g) + xSa(g,m)*phi(i,j,k,g)
                prod(g) = prod(g) + xChi(g,m)*fis(i,j,k)/keff
                do gp = 1, ng
                   if (gp /= g) then
                      aout(g) = aout(g) + xSs(g,gp,m)*phi(i,j,k,g)
                      ains(g) = ains(g) + xSs(gp,g,m)*phi(i,j,k,gp)
                   end if
                end do
             end do
          end do
       end do
    end do
    aleak = prod + ains - aabs - aout
    write(*,'(a)') 'BALANCE_BEGIN'
    write(*,'(a)') 'group,production,inscatter,absorption,outscatter,net_leakage'
    do g = 1, ng
       write(*,'(i0,5(",",es14.6))') g, prod(g), ains(g), aabs(g), aout(g), aleak(g)
    end do
    write(*,'(a)') 'BALANCE_END'
    write(*,'(a,f8.4,a)') ' leakage fraction    = ', 100.0_dp*sum(aleak)/max(sum(prod),1.0e-30_dp), ' %'
    write(*,'(a,f8.4,a)') ' absorption fraction = ', 100.0_dp*sum(aabs)/max(sum(prod),1.0e-30_dp), ' %'
  end subroutine balance_report

  !--------------------------------------------------------------------
  subroutine write_fields(scale)
    real(dp), intent(in) :: scale
    integer :: i, j, k, g, u1, u2
    real(dp) :: xc, yc, zc
    open(newunit=u1, file='flux.csv',  status='replace')
    open(newunit=u2, file='power.csv', status='replace')
    write(u1,'(a)',advance='no') 'x_cm,y_cm,z_cm,mat'
    do g = 1, ng
       write(u1,'(a,i0)',advance='no') ',phi_g', g
    end do
    write(u1,*)
    write(u2,'(a)') 'x_cm,y_cm,z_cm,mat,power_rel'
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             xc = (real(i,dp)-0.5_dp)*dx
             yc = (real(j,dp)-0.5_dp)*dy
             zc = (real(k,dp)-0.5_dp)*dz
             write(u1,'(f12.4,",",f12.4,",",f12.4,",",i0)',advance='no') &
                  xc, yc, zc, matmap(i,j,k)
             do g = 1, ng
                write(u1,'(",",es14.6)',advance='no') phi(i,j,k,g)*scale
             end do
             write(u1,*)
             write(u2,'(f12.4,",",f12.4,",",f12.4,",",i0,",",es14.6)') &
                  xc, yc, zc, matmap(i,j,k), fis(i,j,k)*scale
          end do
       end do
    end do
    close(u1); close(u2)
  end subroutine write_fields

end module cf

!=======================================================================
program coreforge
  use cf
  use ieee_arithmetic, only: ieee_is_nan
  implicit none
  character(len=256) :: infile, envs
  integer  :: outit, g, it, nfuel, iloc(3), outers, i, j, k
  integer(8) :: c0, c1, crate
  real(dp) :: keff, knew, rk, rs, fsum, dre, rkm1, tsec, fxy, dphi
  logical  :: conv

  if (command_argument_count() < 1) then
     write(error_unit,'(a)') 'usage: coreforge <input.txt>  |  coreforge --version'
     error stop 1
  end if
  call get_command_argument(1, infile)
  if (trim(infile) == '--version') then
     write(*,'(a)') 'COREFORGE 8.6 (2D/3D multigroup diffusion engine)'
     stop
  end if

  write(*,'(a)') '======================================================'
  write(*,'(a)') ' COREFORGE 5.0  --  2D/3D multigroup diffusion (k_eff)'
  write(*,'(a)') '======================================================'

  call read_input(infile)
  call setup_coefficients()

  allocate(phi(nx,ny,nz,ng), fis(nx,ny,nz), fisold(nx,ny,nz), src(nx,ny,nz), &
           fuelmask(nx,ny,nz))
  do k = 1, nz
     do j = 1, ny
        do i = 1, nx
           fuelmask(i,j,k) = matfis(matmap(i,j,k))
        end do
     end do
  end do
  phi = 1.0_dp

  call build_fission(fsum)
  if (fsum <= 0.0_dp) call fatal('no fissile material in the core (all nuSf = 0)')
  fis = fis/fsum;  phi = phi/fsum

  call get_environment_variable('OMP_NUM_THREADS', envs)
  if (nz > 1) then
     write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0)') ' mesh ', nx, ' x ', ny, ' x ', nz, &
          '   groups ', ng, '   materials ', nmat
  else
     write(*,'(a,i0,a,i0,a,i0,a,i0)') ' mesh ', nx, ' x ', ny, &
          '   groups ', ng, '   materials ', nmat
  end if
  write(*,'(a,6i2,a,f7.4,a,f5.2)') ' bc(W E S N B T) ', bc, '   gamma ', gam, &
       '   omega ', omega
  if (len_trim(envs) > 0) write(*,'(a,a)') ' OMP_NUM_THREADS = ', trim(envs)

  keff = 1.0_dp;  conv = .false.;  dre = 0.0_dp;  rkm1 = -1.0_dp;  outers = 0
  call system_clock(c0, crate)

  do outit = 1, maxout
     fisold = fis
     do g = 1, ng
        call build_source(g, keff)
        ! adaptive inner solve: sweep until the within-group linear
        ! RESIDUAL is below intol (a true convergence measure) or the cap
        ! is hit, with a floor of `ninner` sweeps.  A fixed, too-small
        ! sweep count lets the OUTER residual fall below tol while the flux
        ! is still un-converged -> a FALSE "converged" whose error GROWS
        ! with mesh refinement.  Converging the inner solve makes the outer
        ! rk/rs honest, restoring clean 2nd-order mesh convergence.
        do it = 1, ninmax
           call rbsor(g, dphi)
           if (it >= ninner .and. grp_resid(g) <= intol) exit
        end do
     end do
     call build_fission(fsum)
     knew = keff * fsum
     rk   = abs(knew - keff)
     fis  = fis/fsum;  phi = phi/fsum
     rs   = maxval(abs(fis - fisold)) / max(maxval(fis), 1.0e-30_dp)
     keff = knew
     outers = outit
     if (ieee_is_nan(keff)) call fatal('k_eff diverged (NaN) -- try a lower OMEGA')
     if (outit > 20 .and. rkm1 > 1.0e-30_dp) dre = 0.9_dp*dre + 0.1_dp*min(rk/rkm1, 1.0_dp)
     rkm1 = rk
     if (mod(outit,25) == 0) &
        write(*,'(a,i6,a,f12.7,a,es9.2,a,es9.2)') ' outer', outit, '  k =', keff, &
              '  dk =', rk, '  dS =', rs
     if (outit >= 10 .and. rk < tolk .and. rs < tols) then
        conv = .true.
        exit
     end if
  end do

  call system_clock(c1)
  tsec = real(c1-c0,dp)/real(crate,dp)

  nfuel = count(fuelmask)
  fxy   = maxval(fis, mask=fuelmask) * real(nfuel,dp)
  iloc  = maxloc(fis, mask=fuelmask)

  write(*,'(a)') '------------------------------------------------------'
  write(*,'(a,f12.7)')      ' KEFF = ', keff
  write(*,'(a,a)')          ' CONVERGED = ', merge('YES','NO ', conv)
  write(*,'(a,i0)')         ' OUTERS = ', outers
  write(*,'(a,f10.4)')      ' FXY = ', fxy
  write(*,'(a,i0,1x,i0,1x,i0)') ' FXYLOC = ', iloc(1), iloc(2), iloc(3)
  write(*,'(a,f8.4)')       ' DR = ', dre
  write(*,'(a,f10.2)')      ' TIME_S = ', tsec
  write(*,'(a)') '------------------------------------------------------'
  if (.not. conv) write(error_unit,'(a)') 'WARNING: not converged within MAXOUT outers'

  call balance_report(keff)
  call write_fields(real(nfuel,dp))
  write(*,'(a)') ' wrote flux.csv, power.csv'
end program coreforge
