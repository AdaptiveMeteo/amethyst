!
! Copyright (c) 2013 Paolo Antonelli, Tiziana Cherubini, Graziano Giuliani
! Original code from oss forward model (ir). aer inc. 2004
! No copyright on source file known.
!
! Permission is hereby granted, free of charge, to any person obtaining a
! copy of this software and associated documentation files (the "Software"),
! to deal in the Software without restriction, including without limitation
! the rights to use, copy, modify, merge, publish, distribute, sublicense,
! and/or sell copies of the Software, and to permit persons to whom the
! Software is furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included
! in all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
! OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
! THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
! DEALINGS IN THE SOFTWARE.
!
module oss_ir

  use iso_fortran_env , only : error_unit

  implicit none

  private

  ! -----------------------------------------------
  ! Public Interfaces which can be called by python
  ! ===============================================

  ! Setup of the OSS Forward Model. To be called in the below sequence
  public :: set_imols
  public :: set_solar_irradiance
  public :: set_hitran
  public :: set_hitran_absorption_coefficients

  ! Obtain the radiances and Jacobians from the OSS Forward Model
  public :: ossdrv_ir

  ! Release memory allocated by the module and cleanup
  public :: release

  ! Internal. Should not be called.
  public :: fatal

  integer , parameter :: stderr_ftn = error_unit

  !integer , parameter :: rk8 = selected_real_kind(P=13,R=300)
  !integer , parameter :: rk4 = selected_real_kind(P= 6,R=37)

  !
  ! Dimension parameter. To allow greater values, need to recompile
  ! the code.
  !
  !P.A. changes for IASI
  !integer , parameter :: mxlev = 101
  integer , parameter :: mxlev = 201
  integer , parameter :: mxlay = mxlev - 1
  integer , parameter :: mxhmol = 25
  integer , parameter :: mxmols = 20

  ! 2010 CODATA
  !real(4) , parameter :: h = 6.62606957E-27
  !real(4) , parameter :: c = 2.99792458E+10
  !real(4) , parameter :: b = 1.380662E-16
  !real(4) , parameter :: c1 = 2.0E0*h*c*c
  !real(4) , parameter :: c2 = h*c/b
  !real(4) , parameter :: pi = 3.1415926535897932384626433832795029E0
  !real(4) , parameter :: deg2rad = pi/180.0E0
  !real(4) , parameter :: rad2deg = 180.0E0/pi
  !real(4) , parameter :: grav = 9.80665E0
  !real(4) , parameter :: rair = 10.0E0/grav

  ! Original constants
  real(4) , parameter :: h = 6.626176E-27
  real(4) , parameter :: c = 2.997925E+10
  real(4) , parameter :: b = 1.380662E-16
  real(4) , parameter :: c1 = 2.0E0*h*c*c
  real(4) , parameter :: c2 = h*c/b
  real(4) , parameter :: pi = 3.14159265
  real(4) , parameter :: deg2rad = 0.017453293E0
  real(4) , parameter :: rair = 1.020408163E0

  ! Weighting factor derivative wrt optical depth, for linear-in-tau,
  ! in low-OD limit
  real(4) , parameter :: dvint = 10.0E0
  real(4) , parameter :: epsiln = 1.0E-05

  ! Filter out
  real(4) , parameter :: odfac = 0.0E0

  integer :: ntmpod = -1
  integer :: nlayod = -1
  integer :: nlev = -1
  integer :: nmol = -1
  integer :: nfsmp = -1
  integer :: nchmax = -1

  integer :: ipsfcg = -1
  integer :: itempg = -1
  integer :: itsking = -1

  ! Internal dynamical storage

  integer , dimension(:) , allocatable :: imolind
  integer , dimension(:) , allocatable :: imolid
  real(4) , dimension(:) , allocatable :: pref
  real(4) , dimension(:) , allocatable :: pavlref
  real(4) , dimension(:,:) , allocatable :: tmptab
  real(4) , dimension(:) , allocatable :: sunrad
  real(4) , dimension(:) , allocatable :: vwvn

  real(4) , dimension(:,:) , allocatable :: coef
  real(4) , dimension(:,:,:) , allocatable :: kfix
  real(4) , dimension(:,:,:) , allocatable :: kh2o
  real(4) , dimension(:,:,:) , allocatable :: dkh2o
  real(4) , dimension(:,:,:,:) , allocatable :: kvar

  integer , dimension(:) , allocatable :: nmols
  integer , dimension(:,:) , allocatable :: imols
  integer , dimension(:) , allocatable :: nch
  integer , dimension(:,:) , allocatable :: ichmap

  logical , parameter :: lookup = .false.
  logical , parameter :: linintau = .false.

  ! 'variablegrid' specifies whether the input is on the absorption
  !  coefficient grid (i.e. normal operating mode) or on a user-specified grid.
  logical :: variablegrid = .false.
  logical :: is_planck_set = .false.

  ! Internal Static storage. See dimension parameters above
  real(4) , dimension(mxlay) :: ap1
  real(4) , dimension(mxlay) :: ap2
  real(4) , dimension(mxlay) :: at1_p1
  real(4) , dimension(mxlay) :: at2_p1
  real(4) , dimension(mxlay) :: at1_p2
  real(4) , dimension(mxlay) :: at2_p2
  real(4) , dimension(mxlay) :: adt1_p1
  real(4) , dimension(mxlay) :: adt2_p1
  real(4) , dimension(mxlay) :: adt1_p2
  real(4) , dimension(mxlay) :: adt2_p2
  integer , dimension(mxlay) :: indxt_p1
  integer , dimension(mxlay) :: indxt_p2
  integer , dimension(mxlay) :: indxp
  real(4) , dimension(mxlay) :: tavl
  real(4) , dimension(mxlay) :: pavl
  real(4) , dimension(mxlay) :: wfix
  real(4) , dimension(mxlay) :: dtu
  real(4) , dimension(mxlay) :: dtl
  real(4) , dimension(mxhmol,mxlay) :: q
  real(4) , dimension(mxhmol,mxlay) :: w
  real(4) , dimension(mxlay,mxhmol) :: dwqu
  real(4) , dimension(mxlay,mxhmol) :: dwql
  real(4) , dimension(mxhmol) :: qobs
  real(4) , dimension(mxhmol) :: wobs
  real(4) , dimension(mxhmol) :: dwquobs
  real(4) , dimension(mxhmol) :: dwqlobs
  real(4) , dimension(mxmols,mxlay) :: abso
  real(4) , dimension(mxlay) :: tautot
  real(4) , dimension(mxlay) :: dtaudtmp
  real(4) , dimension(mxmols) :: absoobs
  real(4) , dimension(mxlev) :: txdn
  real(4) , dimension(mxlev) :: txup
  real(4) , dimension(mxlev) :: bbar
  real(4) , dimension(mxlev) :: dbbar
  real(4) , dimension(mxlev) :: draddtmp
  real(4) , dimension(mxlev) :: draddtau
  real(4) , dimension(mxlev) :: drdw
  real(4) , dimension(mxlev) :: draddtmpdw
  real(4) , dimension(mxlev) :: draddtmpuw
  real(4) , dimension(mxlev) :: dp
  real(4) , dimension(mxlev) :: qcor
  real(4) , dimension(mxlev) :: qr
  real(4) , dimension(mxlay) :: ploc
  real(4) , dimension(mxlay) :: b2
  real(4) , dimension(mxlay) :: db2
  real(4) , dimension(mxlay) :: ar
  real(4) , dimension(mxlay) :: br
  real(4) , dimension(mxlay) :: ad
  real(4) , dimension(mxlay) :: bd
  real(4) :: tautotobs
  real(4) :: dtaudtmpobs
  real(4) :: tavlobs
  real(4) :: wfixobs
  real(4) :: ap1obs
  real(4) :: ap2obs
  real(4) :: dtuobs
  real(4) :: dtlobs
  real(4) :: at1_p1obs
  real(4) :: at2_p1obs
  real(4) :: adt1_p1obs
  real(4) :: adt2_p1obs
  real(4) :: v2
  real(4) :: bs2
  real(4) :: dbs2
  real(4) :: ars
  real(4) :: brs
  real(4) :: ads
  real(4) :: bds

  integer :: indxt_p1obs = 0
  integer :: indxpobs

  contains

  subroutine fatal(f,l,message)
    implicit none
    character(len=*) , intent(in) :: f , message
    integer , intent(in) :: l
    write(stderr_ftn,'(a,a,a,i8,a,a)') 'Fatal in file: ',trim(f), &
      ' at line ',l,' : ',message
    stop
  end subroutine fatal

  subroutine ossdrv_ir(nparg,nchan,nspe,nsf,nuser,initempg,initsking, &
      inipsfcg,xg,pobs,obsang,sunang,sfgrd,remrf,y,xkt,xkemrf,paxkemrf,puser)
    !---------------------------------------------------------------------------
    ! purpose: Driver for the ir radiative transfer model.
    !          Computes both radiances and their jacobians wrt geophysical
    !          parameters.
    ! Input:
    !   xg       Profile vector of geophysical parmaters (containing
    !            temperature and constituents profiles, surface
    !            pressure and skin temperature, and cloud parameters
    !   remrf    Vector of mw surface emissivities
    !   pobs     Pressure at observation point
    !   obsang   Angle of observation
    !   sunang   Sun angle
    !   puser    Vertical profile pressure levels
    ! Output:
    !   y        Vector of radiances.
    !   xkt      Array of derivatives of the radiances wrt geophysical
    !            parameters.
    !   xkemrf   Array of radiance derivatives wrt surface emissivities.
    !---------------------------------------------------------------------------
    implicit none
    integer , intent(in) :: nparg , nchan , nsf , nspe , nuser
    real(4) , intent(in) , dimension(nparg) :: xg
    real(4) , intent(in) , dimension(nsf) :: sfgrd
    real(4) , intent(in) , dimension(nsf) :: remrf
    integer , intent(in) :: inipsfcg , initempg , initsking
    real(4) , intent(in) :: obsang , sunang , pobs
    real(4) , intent(in) , dimension(nuser) , optional :: puser
    real(4) , intent(inout) , dimension(nchan) :: y
    real(4) , intent(inout) , dimension(nparg,nchan) :: xkt
    real(4) , intent(inout) , dimension(2,nspe,nchan) :: xkemrf
    real(4) , intent(inout) , dimension(2,nchan) :: paxkemrf
    logical :: sun
    real(4) , dimension(mxlev) :: xgtmp
    real(4) , dimension(2) :: vemrf
    real(4) , dimension(2) :: xkemrf_tmp
    real(4) , dimension(nparg) :: xkt_tmp
    real(4) , dimension(2,nsf) :: emrf
    real(4) :: umu , umu0 , delphi , vn , a , fbeam , rad , xx
    real(4) :: tsfc
    integer :: nsurf , nobs , n1 , n2 , ip0 , nn , ich , i , ich0

    if ( nchmax < 0 ) then
      call fatal(__FILE__,__LINE__, &
        'HITRAN COEFFICIENTS HAVE BEEN NOT PROVIDED. CANNOT COMPUTE')
    end if

    ipsfcg = inipsfcg
    itempg = initempg
    itsking = initsking
    if ( imolind(1) < 1 ) then
      ! Assume "default" ordering
      do i = 1 , size(imolid)
        imolind(i) = ipsfcg+1+(i-1)*(itsking-1)
      end do
    end if

    emrf(1,:) = remrf
    emrf(2,:) = 1.0E0-emrf(1,:)

    !======================================================================
    !     initialize radiance vector and k-matrix
    !======================================================================
    y(1:nchan) = 0.0E0
    xkt(1:nparg,1:nchan) = 0.0E0
    xkemrf(1:2,1:nsf,1:nchan) = 0.0E0
    paxkemrf(1:2,1:nchan) = 0.0E0

    if ( present(puser) )then
      if ( nuser > mxlev ) then
        call fatal(__FILE__,__LINE__, &
          'err[oss_ir_module::ossdrv]: input pressure grid is too large.')
      end if
      variablegrid = .true.
    else ! initializes to zero to avoid nan.
      variablegrid = .false.
    end if
    !======================================================================
    !     compute path variables
    !======================================================================
    !     compute path geometry
    if ( variablegrid ) then
      ! n1 is defined by pref(n1-1) < pobs < pref(n1)
      call setpath_ir_vg(pobs,obsang,sunang,nsurf,nobs,    &
                         n1,n2,sun,umu,umu0,delphi,puser)
      !---compute average temperature and molecular amounts for the layers
      call fpath_ir_vg(xg,pobs,n1,n2,tsfc,puser)
      !---compute coefficients for temperature interpolation of ods
      call settabindx_ir_vg(n2)
    else
      ! n1 is defined by pref(n1-1) < pobs < pref(n1)
      call setpath_ir(pref,xg,pobs,obsang,sunang,nsurf,nobs,    &
                      n1,n2,sun,umu,umu0,delphi)
      !---compute average temperature and molecular amounts for the layers
      call fpath_ir(pref,xg,pobs,n1,n2,tsfc)
      call settabindx_ir(n2)
    end if
    if ( n1 > 1 ) then
      call settabindx_irobs(n1)
    end if
    !======================================================================
    !     loop over spectral points
    !======================================================================
    ip0 = 2
    vn = vwvn(1)-0.01E0
    if ( linintau .and. variablegrid ) then
      call planck_set(vn,xg(1:n2+1),xg(itsking),n2+1)
    else if ( linintau .and. .not. variablegrid ) then
      xgtmp(1:n2) = xg(1:n2)
      xgtmp(1+n2) = tsfc
      call planck_set(vn,xgtmp,xg(itsking),n2+1)
    else
      call planck_set(vn,tavl,xg(itsking),n2)
    endif
    if ( variablegrid ) then
      do nn = 1 , nfsmp
        !---compute molecular optical depth for all atmospheric layers
        call osstran_vg(n2,nn)
        !---interpolate input surface emissivity to node wavenumber
        vn = vwvn(nn)
        call vinterp(emrf,vn,sfgrd,nsf,2,ip0,a,vemrf)
        fbeam = sunrad(nn)
        if ( .not. sun ) fbeam = 0.0E0
        !---perform rt calculations !clear sky model
        call ossrad_vg(nn,xg,vemrf,fbeam,vn,n1,n2,sun,umu,umu0, &
                       rad,xkt_tmp,xkemrf_tmp)
        do ich = 1 , nch(nn)
          ich0 = ichmap(ich,nn)
          y(ich0) = y(ich0) + rad*coef(ich,nn)
          xkt(1:nparg,ich0) = xkt(1:nparg,ich0) + xkt_tmp(1:nparg)*coef(ich,nn)
          do i = 1 , 2
            xx = xkemrf_tmp(i)*coef(ich,nn)
            xkemrf(i,ip0-1,ich0) = xkemrf(i,ip0-1,ich0) + xx*(1.0E0-a)
            xkemrf(i,ip0,ich0) = xkemrf(i,ip0,ich0) + xx*a
            paxkemrf(i,ich0) = paxkemrf(i,ich0) + xx
          end do
        end do
      end do
    else
      do nn = 1 , nfsmp
        !---compute molecular optical depth for all atmospheric layers
        tautotobs = 0.0E0
        call osstran(n1,n2,nn)
        !---interpolate input surface emissivity to node wavenumber
        vn = vwvn(nn)
        call vinterp(emrf,vn,sfgrd,nsf,2,ip0,a,vemrf)
        fbeam = sunrad(nn)
        if ( .not. sun ) fbeam = 0.0E0
        !---perform rt calculations !clear sky model
        call ossrad(nn,xg,vemrf,fbeam,vn,n1,n2,sun,umu,umu0,rad,   &
                    xkt_tmp,xkemrf_tmp)
        do ich = 1 , nch(nn)
          ich0 = ichmap(ich,nn)
          y(ich0) = y(ich0) + rad*coef(ich,nn)
          xkt(1:nparg,ich0) = xkt(1:nparg,ich0) + &
                              xkt_tmp(1:nparg)*coef(ich,nn)
          do i = 1 , 2
            xx = xkemrf_tmp(i)*coef(ich,nn)
            xkemrf(i,ip0-1,ich0) = xkemrf(i,ip0-1,ich0) + xx*(1.0E0-a)
            xkemrf(i,ip0,ich0) = xkemrf(i,ip0,ich0) + xx*a
            paxkemrf(i,ich0) = paxkemrf(i,ich0) + xx
          end do
        end do
      end do
    end if

  contains

  subroutine planck_set(vn,t,ts,nlay)
    integer , intent(in) :: nlay
    real(4) , intent(in) :: vn , ts
    real(4) , intent(in) , dimension(:) :: t
    integer :: l
    if ( nlay > mxlay ) then
      call fatal(__FILE__,__LINE__,'nlay exceeds mxlay')
    end if
    v2 = vn
    do l = 1 , nlay
      call fplanck(vn,t(l),b2(l),db2(l))
    end do
    call fplanck(vn,ts,bs2,dbs2)
    is_planck_set = .true.
  end subroutine planck_set

  subroutine planck_int(vn,t,ts,nlay,radpl,db,bs,dbs)
    implicit none
    integer , intent(in) :: nlay
    real(4) , intent(in) :: vn , ts
    real(4) , intent(in) , dimension(:) :: t
    real(4) , intent(out) :: bs , dbs
    real(4) , intent(out) , dimension(:) :: radpl , db
    real(4) :: bs1 , dbs1
    real(4) :: v1 , v2n
    real(4) , dimension(nlay) :: b1 , db1
    if ( nlay > mxlay ) then
      call fatal(__FILE__,__LINE__,'nlay exceeds mxlay')
    end if
    if ( .not. is_planck_set ) then
      call planck_set(vn,t,ts,nlay)
    end if
    if ( vn >= v2 ) then
      v2n = v2 + dvint
      if ( vn >= v2n ) then
        call planck_set(vn,t,ts,nlay)
        v2 = vn
        v2n = v2+dvint
      end if
      v1 = v2
      b1(1:nlay)  = b2(1:nlay)
      db1(1:nlay) = db2(1:nlay)
      bs1 = bs2
      dbs1 = dbs2
      v2 = v2n
      call planck_set(v2,t,ts,nlay)
      ar(1:nlay) = (b2(1:nlay)-b1(1:nlay))/dvint
      br(1:nlay) = (b1(1:nlay)*v2-b2(1:nlay)*v1)/dvint
      ad(1:nlay) = (db2(1:nlay)-db1(1:nlay))/dvint
      bd(1:nlay) = (db1(1:nlay)*v2-db2(1:nlay)*v1)/dvint
      ars = (bs2-bs1)/dvint
      brs = (bs1*v2-bs2*v1)/dvint
      ads = (dbs2-dbs1)/dvint
      bds = (dbs1*v2-dbs2*v1)/dvint
    end if
    radpl(1:nlay) = vn*ar(1:nlay)+br(1:nlay)
    db(1:nlay) = vn*ad(1:nlay)+bd(1:nlay)
    bs = vn*ars+brs
    dbs = vn*ads+bds
  end subroutine planck_int

  subroutine fplanck(vn,t,rad,draddt)
    !------------------------------------------------------------------------
    ! purpose: calculates planck function and its derivative  wrt temperature
    !------------------------------------------------------------------------
    implicit none
    real(4) , intent(in) :: t , vn
    real(4) , intent(out) :: rad , draddt
    real(4) :: f1 , f2 , f3
    f1 = c1*vn**3
    f2 = c2*vn/t
    f3 = exp(f2)
    rad = f1/(f3-1.0E0)
    draddt =(rad*rad)*f2*f3/(f1*t)
  end subroutine fplanck

  subroutine lint_log(xinp,pgrid,n0,p0,x0,x1)
    !---------------------------------------
    ! purpose: interpolation in log pressure
    !---------------------------------------
    integer , intent(in) :: n0
    real(4) , intent(in) :: p0
    real(4) , dimension(:) , intent(in) :: xinp
    real(4) , dimension(:) , intent(in) :: pgrid
    real(4) , intent(out) :: x0
    real(4) , intent(out) , optional :: x1
    real(4) :: xx
    xx = log(xinp(n0)/xinp(n0-1))/log(pgrid(n0)/pgrid(n0-1))
    x0 = xinp(n0-1)*(p0/pgrid(n0-1))**xx
    if ( present(x1) ) x1 = xinp(n0)*(p0/pgrid(n0))**xx
  end subroutine lint_log

  subroutine lpsum_log(pu,pl,xu,xl,scal,xint,dxu,dxl)
    !------------------------------------------------------------------
    ! purpose: computes average layer quantities (or integrated amount)
    !          using a log-x dependence on log-p.
    !------------------------------------------------------------------
    implicit none
    real(4) , intent(in) :: pu , pl , scal
    real(4) , intent(in) :: xu , xl
    real(4) , intent(out) :: xint , dxu , dxl
    real(4) :: hp , x0 , zeta , alza , alpha
    hp = log(pl/pu)
    x0 = pl*xl*hp
    zeta = pu*xu/(pl*xl)
    if ( abs(zeta-1.0E0) > epsiln ) then
      alza  = log(zeta)
      xint  = x0*(zeta-1.0E0)/alza
      alpha = zeta/(zeta-1.0E0)-1.0E0/alza
    else
      xint  = x0*2.0E0/(3.0E0-zeta)
      alpha = zeta/(3.0E0-zeta)
    end if
    xint = xint*scal
    dxu = xint*alpha/xu
    dxl = xint*(1.0E0-alpha)/xl
  end subroutine lpsum_log

  subroutine ossrad_vg(nn,xg,emrf,bs_sun,vn,n1,n2,sun,umu,umu0,rad,xkt,xkemrf)
    !------------------------------------------------------------------
    ! purpose: compute radiances (in mw/m2/str/cm-1) and derivatives of
    !   radiances with respect to atmospheric and surface parameters
    !------------------------------------------------------------------
    implicit none
    logical , intent(in) :: sun
    integer , intent(in) :: n1 , n2 , nn
    real(4) , intent(in) :: umu
    real(4) , intent(in) :: umu0
    real(4) , intent(in) :: bs_sun
    real(4) , intent(in) :: vn
    real(4) , dimension(:) , intent(in) :: xg
    real(4) , dimension(:) , intent(in) :: emrf
    real(4) , intent(out) :: rad
    real(4) , intent(out) , dimension(:) :: xkt , xkemrf
    integer :: l , ks , k , ixoff
    real(4) :: sec , secdif , sec0 , em , sun_rsfc , bs , dbs
    real(4) :: tausfc , radsun , draddrsfc , draddtau_sun , draddemis
    real(4) :: draddtskn , dtran , txsun , rsfc
    real(4) :: sumtau0_up = 0.0E0
    real(4) :: sumtau_dwn = 0.0E0
    sec = 1.0E0/umu
    sec0 = 1.0E0/umu0
    secdif = sec
    em = emrf(1)
    sun_rsfc = emrf(2)
    !-----------------------------------------------------------------------
    !     compute planck function and its derivative wrt temperature
    !-----------------------------------------------------------------------
    ! linintau is skipped in this version
    call planck_int(vn,tavl,xg(itsking),n2,bbar,dbbar,bs,dbs)
    !-----------------------------------------------------------------------
    !     compute transmittance profile along viewing path down to surface
    !-----------------------------------------------------------------------
    txdn(1) = 1.0E0
    sumtau_dwn = 0.0E0
    if ( n1 > 1 ) then
      txdn(1:n1-1) = 1.0E0
      sumtau_dwn = tautotobs*sec
      txdn(n1) = exp(-sumtau_dwn)
    end if
    do l = n1 , n2
      sumtau_dwn = sumtau_dwn+tautot(l)*sec
      txdn(l+1) = exp(-sumtau_dwn)
    end do
    tausfc = txdn(n2+1)
    !---initialize radiance and derivative arrays
    rad = 0.0E0
    radsun = 0.0E0
    draddrsfc = 0.0E0
    draddtau_sun = 0.0E0
    draddtmp(1:n2) = 0.0E0
    draddtau(1:n2) = 0.0E0
    draddtmpdw(1:n2+1)= 0.0E0
    draddtmpuw(1:n2+1)= 0.0E0
    draddemis = 0.0E0
    draddtskn = 0.0E0
    !-----------------------------------------------------------------------
    !     1- downwelling thermal radiance calculation:
    !-----------------------------------------------------------------------
    if ( tausfc > 1.0E-06 ) then
      txup(n2+1) = tausfc
      sumtau0_up = 0.0E0
      do l = n2 , 1 , -1
        sumtau0_up = sumtau0_up+tautot(l)
        txup(l) = exp(-(sumtau_dwn+sumtau0_up*secdif))
      end do
      do l = 1 , n2
        dtran = txup(l+1)-txup(l)
        !--- derivative of downwelling emission wrt to
        !--- temperature and constituents
        draddtmp(l) = dtran*dbbar(l)
        draddtau(l) = (txup(l)*bbar(l)-rad)*secdif
        rad = rad+dtran*bbar(l)
      end do
    end if
    !-----------------------------------------------------------------------
    !     2- add solar component
    !-----------------------------------------------------------------------
    if ( sun ) then
      txsun = exp(-(sumtau_dwn+sumtau0_up*sec0))
      !---derivatives wrt sfc solar reflectance
      draddrsfc = txsun*bs_sun*umu0/pi
      radsun = sun_rsfc*draddrsfc
      draddtau_sun = -radsun*(sec+sec0)
    end if
    !-----------------------------------------------------------------------
    !     3- add surface terms
    !-----------------------------------------------------------------------
    rsfc = (1.0E0-em)
    !---derivatives wrt emissivity and sfc skin temperature:
    draddemis = tausfc*bs-rad
    draddtskn = em*tausfc*dbs
    rad = rad*rsfc+em*tausfc*bs+radsun
    !-----------------------------------------------------------------------
    !     4- upwelling thermal radiance calculation
    !-----------------------------------------------------------------------
    do l = n2 , n1 , -1
      dtran = txdn(l)-txdn(l+1)
      draddtau(l) = draddtau(l)*rsfc+(txdn(l+1)*bbar(l)-rad)*sec+draddtau_sun
      draddtmp(l) = draddtmp(l)*rsfc+dtran*dbbar(l)+draddtau(l)*dtaudtmp(l)
      rad = rad+dtran*bbar(l)
    end do
    do l = 1 , n1-1
      draddtau(l) = draddtau(l)*rsfc
      draddtmp(l) = draddtmp(l)*rsfc
    end do
    !-----------------------------------------------------------------------
    !     compute level derivatives and and map to array xkt:
    !-----------------------------------------------------------------------
    xkt = 0.0E0
    !  air temperature
    do l = n2 , 1 , -1
      xkt(itempg+l) = xkt(itempg+l)+draddtmp(l)*dtl(l)
      xkt(itempg+l-1) = xkt(itempg+l-1)+draddtmp(l)*dtu(l)
    end do
    !---molecular concentrations
    do ks = 1 , nmols(nn)
      k = imols(ks,nn)
      ixoff = imolind(k)-1
      do l = n2 , 1 , -1
        drdw(l) = draddtau(l)*abso(ks,l)
        xkt(ixoff+l+1) = xkt(ixoff+l+1)+drdw(l)*dwql(l,k)
        xkt(ixoff+l) = xkt(ixoff+l)+drdw(l)*dwqu(l,k)
      end do
    end do
    !---surface terms
    xkt(itsking) = draddtskn !--tskin
    xkemrf(1) = draddemis
    xkemrf(2) = draddrsfc
  end subroutine ossrad_vg

  subroutine ossrad(nn,xg,emrf,bs_sun,vn,n1,n2,sun,umu,umu0,   &
       rad,xkt,xkemrf)
    !------------------------------------------------------------------
    ! purpose: compute radiances (in mw/m2/str/cm-1) and derivatives of
    !   radiances with respect to atmospheric and surface parameters
    !------------------------------------------------------------------
    implicit none
    logical , intent(in) :: sun
    integer , intent(in) :: n1 , n2 , nn
    real(4) , intent(in) :: umu
    real(4) , intent(in) :: umu0
    real(4) , intent(in) :: bs_sun
    real(4) , intent(in) :: vn
    real(4) , dimension(:) , intent(in) :: xg
    real(4) , dimension(:) , intent(in) :: emrf
    real(4) , intent(out) :: rad
    real(4) , intent(out) , dimension(:) :: xkt , xkemrf
    integer :: l , ks , k , ixoff
    real(4) :: sec , secdif , sec0 , em , sun_rsfc , bs , dbs
    real(4) :: tausfc , radsun , draddrsfc , draddtau_sun , draddemis
    real(4) :: draddtskn , dtran , txsun , rsfc
    real(4) :: bobs , dbobs
    real(4) :: drdwobs2
    real(4) :: draddtmpobs1 = 0.0E0
    real(4) :: draddtmpobs2 = 0.0E0
    real(4) :: draddtauobs1 = 0.0E0
    real(4) :: draddtauobs2 = 0.0E0
    real(4) :: sumtau0_up = 0.0E0
    real(4) :: sumtau_dwn = 0.0E0
    sec = 1.0E0/umu
    sec0 = 1.0E0/umu0
    secdif = sec
    em = emrf(1)
    sun_rsfc = emrf(2)
    !-----------------------------------------------------------------------
    !     compute planck function and its derivative wrt temperature
    !-----------------------------------------------------------------------
    ! linintau is skipped in this version
    call planck_int(vn,tavl,xg(itsking),n2,bbar,dbbar,bs,dbs)
    !-----------------------------------------------------------------------
    !     compute transmittance profile along viewing path down to surface
    !-----------------------------------------------------------------------
    txdn(1) = 1.0E0
    sumtau_dwn = 0.0E0
    if ( n1 > 1 ) then
      txdn(1:n1-1) = 1.0E0
      sumtau_dwn = tautotobs*sec
      txdn(n1) = exp(-sumtau_dwn)
    end if
    do l = n1 , n2
      sumtau_dwn = sumtau_dwn+tautot(l)*sec
      txdn(l+1) = exp(-sumtau_dwn)
    end do
    tausfc = txdn(n2+1)
    !---initialize radiance and derivative arrays
    rad = 0.0E0
    radsun = 0.0E0
    draddrsfc = 0.0E0
    draddtau_sun = 0.0E0
    draddtmp(1:n2) = 0.0E0
    draddtau(1:n2) = 0.0E0
    draddtmpdw(1:n2+1)= 0.0E0
    draddtmpuw(1:n2+1)= 0.0E0
    draddemis = 0.0E0
    draddtskn = 0.0E0
    !-----------------------------------------------------------------------
    !     1- downwelling thermal radiance calculation:
    !-----------------------------------------------------------------------
    if ( tausfc > 1.0E-06 ) then
      txup(n2+1) = tausfc
      sumtau0_up = 0.0E0
      do l = n2 , 1 , -1
        sumtau0_up = sumtau0_up+tautot(l)
        txup(l) = exp(-(sumtau_dwn+sumtau0_up*secdif))
      end do
      do l = 1 , n2
        dtran = txup(l+1)-txup(l)
        !--- derivative of downwelling emission wrt to
        !--- temperature and constituents
        draddtmp(l) = dtran*dbbar(l)
        draddtau(l) = (txup(l)*bbar(l)-rad)*secdif
        rad = rad+dtran*bbar(l)
      end do
    end if
    !-----------------------------------------------------------------------
    !     2- add solar component
    !-----------------------------------------------------------------------
    if ( sun ) then
      txsun = exp(-(sumtau_dwn+sumtau0_up*sec0))
      !---derivatives wrt sfc solar reflectance
      draddrsfc = txsun*bs_sun*umu0/pi
      radsun = sun_rsfc*draddrsfc
      draddtau_sun = -radsun*(sec+sec0)
    end if
    !-----------------------------------------------------------------------
    !     3- add surface terms
    !-----------------------------------------------------------------------
    rsfc = (1.0E0-em)
    !---derivatives wrt emissivity and sfc skin temperature:
    draddemis = tausfc*bs-rad
    draddtskn = em*tausfc*dbs
    rad = rad*rsfc+em*tausfc*bs+radsun
    !-----------------------------------------------------------------------
    !     4- upwelling thermal radiance calculation
    !-----------------------------------------------------------------------
    do l = n2 , n1 , -1
      dtran = txdn(l)-txdn(l+1)
      draddtau(l) = draddtau(l)*rsfc+(txdn(l+1)*bbar(l)-rad)*sec+draddtau_sun
      draddtmp(l) = draddtmp(l)*rsfc+dtran*dbbar(l)+draddtau(l)*dtaudtmp(l)
      rad = rad+dtran*bbar(l)
    end do
    do l = 1 , n1-1
      draddtau(l) = draddtau(l)*rsfc
      draddtmp(l) = draddtmp(l)*rsfc
    end do
    if ( n1 > 1 ) then
      dtran          = txdn(n1-1)-txdn(n1)
      call fplanck(vn,tavlobs,bobs,dbobs)
      draddtauobs1 = draddtau(l)*rsfc
      draddtauobs2 = (txdn(n1)*bobs-rad)*sec+draddtau_sun
      draddtmpobs1 = draddtmp(l)*rsfc+draddtauobs1*dtaudtmp(n1-1)
      draddtmpobs2 = dtran*dbobs+draddtauobs2*dtaudtmpobs
      rad = rad+dtran*bobs
    end if
    !-----------------------------------------------------------------------
    !     compute level derivatives and and map to array xkt:
    !-----------------------------------------------------------------------
    xkt = 0.0E0
    !  air temperature
    do l = n2 , 1 , -1
      xkt(itempg+l) = xkt(itempg+l)+draddtmp(l)*dtl(l)
      xkt(itempg+l-1) = xkt(itempg+l-1)+draddtmp(l)*dtu(l)
    end do
    if ( n1 > 1 ) then
      xkt(itempg+n1-1) = xkt(itempg+n1-1)+draddtmpobs2*dtlobs
      xkt(itempg+n1-2) = xkt(itempg+n1-2)+draddtmpobs2*dtuobs
    end if
    !---molecular concentrations
    do ks = 1 , nmols(nn)
      k = imols(ks,nn)
      ixoff = imolind(k)-1
      do l = n2 , 1 , -1
        drdw(l) = draddtau(l)*abso(ks,l)
        xkt(ixoff+l+1) = xkt(ixoff+l+1)+drdw(l)*dwql(l,k)
        xkt(ixoff+l) = xkt(ixoff+l)+drdw(l)*dwqu(l,k)
      end do
      if ( n1 > 1 ) then
        ! drdwobs1 = draddtauobs1*abso(ks,n1-1)
        drdwobs2 = draddtauobs2*absoobs(ks)
        xkt(ixoff+n1) = xkt(ixoff+n1)+drdwobs2*dwqlobs(ks)
        xkt(ixoff+n1-1) = xkt(ixoff+n1-1)+drdwobs2*dwquobs(k)
      end if
    end do
    !---surface terms
    xkt(itsking) = draddtskn !--tskin
    xkemrf(1) = draddemis
    xkemrf(2) = draddrsfc
  end subroutine ossrad

  subroutine vinterp(datin,vn,sfgrd,nsf,nxdim,ip0,coefint,datout)
    !-------------------------------------
    ! purpose: interpolation in wavenumber
    !-------------------------------------
    implicit none
    integer , intent(in) :: nxdim , nsf
    real(4) , intent(in) , dimension(:,:) :: datin
    real(4) , intent(in) :: vn
    real(4) , intent(in) , dimension(:) :: sfgrd
    integer , intent(inout) :: ip0
    real(4) , intent(out) :: coefint
    real(4) , intent(out) , dimension(:) :: datout
    integer :: ip , i
    iploop: &
    do ip = ip0 , nsf-1
      if ( sfgrd(ip) > vn ) exit iploop
    end do iploop
    ip0 = ip
    coefint = (vn-sfgrd(ip0-1))/(sfgrd(ip0)-sfgrd(ip0-1))
    do i = 1 , nxdim
      datout(i) = coefint*(datin(i,ip0)-datin(i,ip0-1))+datin(i,ip0-1)
    end do
  end subroutine vinterp

  subroutine settabindx_ir_vg(n2)
    !--------------------------------------------------------------------
    ! purpose: computes temperature/water vapor indexes and interpolation
    !          coefficients.
    !--------------------------------------------------------------------
    implicit none
    integer , intent(in) :: n2
    integer :: l , i , j1 , j2 , lp1 , lp2
    real(4) :: dent1_p1 , dent2_p1 , dent1_p2 , dent2_p2

    !---compute coefficients for temperature interpolation of ods
    !!! indxp: is related to the index of od pressure layer located just above
    !!!  (by altitude) the current profile layer
    do l = 1 , n2
      lp1 = 1
      do while ( pavlref(lp1) < pavl(l) )
        lp1 = lp1+1
      end do
      lp2 = lp1
      lp1 = lp1-1
      if ( lp1 == 0 ) then
        lp1 = 1
        lp2 = 2
      end if
      if ( lp2 >= nlev ) then
        lp1 = nlev-2
        lp2 = nlev-1
      end if
      !!! if pavl(l) becomes located between lp1 and lp2:
      ap2(l) = (pavl(l)-pavlref(lp1))/(pavlref(lp2)-pavlref(lp1))
      ap1(l) = (pavlref(lp2)-pavl(l))/(pavlref(lp2)-pavlref(lp1))
      indxp(l) = lp1
      !!! indxt_p1: is related to the midst temperature index within the
      !!!           upper (by altitude) od pressure layer indxp
      !!! indxt_p2: is related to the midst temperature index within the
      !!!           lower (by altitude) od pressure layer indxp+1
      temploop1: &
      do i = 3 , ntmpod-2
        if ( tmptab(i,lp1) >= tavl(l) ) exit temploop1
      end do temploop1
      if ( abs(tavl(l)-tmptab(i-1,lp1)) <= abs(tavl(l)-tmptab(i,lp1)) ) then
        j1 = i-1
      else
        j1 = i
      end if
      indxt_p1(l) = j1
      dent1_p1 = (tmptab(j1-1,lp1)-tmptab(j1,lp1)) * &
                 (tmptab(j1-1,lp1)-tmptab(j1+1,lp1))
      dent2_p1 = (tmptab(j1,lp1)-tmptab(j1-1,lp1)) * &
                 (tmptab(j1,lp1)-tmptab(j1+1,lp1))
      at1_p1(l) = (tavl(l)-tmptab(j1,lp1)) * &
                  (tavl(l)-tmptab(j1+1,lp1))/dent1_p1
      at2_p1(l) = (tavl(l)-tmptab(j1-1,lp1)) * &
                  (tavl(l)-tmptab(j1+1,lp1))/dent2_p1
      adt1_p1(l) = (2.0E0*tavl(l)-tmptab(j1,lp1)-tmptab(j1+1,lp1))/dent1_p1
      adt2_p1(l) = (2.0E0*tavl(l)-tmptab(j1-1,lp1)-tmptab(j1+1,lp1))/dent2_p1

      temploop2: &
      do i = 3 , ntmpod-2
        if ( tmptab(i,lp2) >= tavl(l) ) exit temploop2
      end do temploop2
      if ( abs(tavl(l)-tmptab(i-1,lp2)) <= abs(tavl(l)-tmptab(i,lp2)) ) then
        j2 = i-1
      else
        j2 = i
      end if
      indxt_p2(l) = j2
      !!! aty_px: interpolation coefficients, which are related to pressure
      !!!         layer (x=1,or,2, upper, or, lower) and temperature
      !!!         interpolation points (y=1,or,2, lower,or, midst)
      dent1_p2 = (tmptab(j2-1,lp2)-tmptab(j2,lp2)) * &
                 (tmptab(j2-1,lp2)-tmptab(j2+1,lp2))
      dent2_p2 = (tmptab(j2,lp2)-tmptab(j2-1,lp2)) * &
                 (tmptab(j2,lp2)-tmptab(j2+1,lp2))
      at1_p2(l) = (tavl(l)-tmptab(j2,lp2)) * &
                  (tavl(l)-tmptab(j2+1,lp2))/dent1_p2
      at2_p2(l) = (tavl(l)-tmptab(j2-1,lp2)) * &
                  (tavl(l)-tmptab(j2+1,lp2))/dent2_p2
      adt1_p2(l) = (2.0E0*tavl(l)-tmptab(j2,lp2)-tmptab(j2+1,lp2))/dent1_p2
      adt2_p2(l) = (2.0E0*tavl(l)-tmptab(j2-1,lp2)-tmptab(j2+1,lp2))/dent2_p2
    end do
  end subroutine settabindx_ir_vg

  subroutine settabindx_ir(n2)
    !--------------------------------------------------------------------
    ! purpose: computes temperature/water vapor indexes and interpolation
    !          coefficients.
    !--------------------------------------------------------------------
    implicit none
    integer , intent(in) :: n2
    integer :: l , i , j1 , lp1
    real(4) :: dent1_p1 , dent2_p1

    !---compute coefficients for temperature interpolation of ods
    !!! indxp: is related to the index of od pressure layer located just above
    !!!  (by altitude) the current profile layer
    do l = 1 , n2
      lp1 = l
      ap2(l) = 0.0E0
      ap1(l) = 1.0E0
      indxp(l) = lp1
      !!! indxt_p1: is related to the midst temperature index within the
      !!!           upper (by altitude) od pressure layer indxp
      !!! indxt_p2: is related to the midst temperature index within the
      !!!           lower (by altitude) od pressure layer indxp+1
      temploop1: &
      do i = 3 , ntmpod-2
        if ( tmptab(i,lp1) >= tavl(l) ) exit temploop1
      end do temploop1
      if ( abs(tavl(l)-tmptab(i-1,lp1)) <= abs(tavl(l)-tmptab(i,lp1)) ) then
        j1 = i-1
      else
        j1 = i
      end if
      indxt_p1(l) = j1
      dent1_p1 = (tmptab(j1-1,lp1)-tmptab(j1,lp1)) * &
                 (tmptab(j1-1,lp1)-tmptab(j1+1,lp1))
      dent2_p1 = (tmptab(j1,lp1)-tmptab(j1-1,lp1)) * &
                 (tmptab(j1,lp1)-tmptab(j1+1,lp1))
      at1_p1(l) = (tavl(l)-tmptab(j1,lp1)) * &
                  (tavl(l)-tmptab(j1+1,lp1))/dent1_p1
      at2_p1(l) = (tavl(l)-tmptab(j1-1,lp1)) * &
                  (tavl(l)-tmptab(j1+1,lp1))/dent2_p1
      adt1_p1(l) = (2.0E0*tavl(l)-tmptab(j1,lp1)-tmptab(j1+1,lp1))/dent1_p1
      adt2_p1(l) = (2.0E0*tavl(l)-tmptab(j1-1,lp1)-tmptab(j1+1,lp1))/dent2_p1
    end do
  end subroutine settabindx_ir

  subroutine settabindx_irobs(n1)
    implicit none
    integer , intent(in) :: n1
    integer :: i , j1 , lp1
    real(4) :: dent1_p1 , dent2_p1

    lp1 = n1-1
    indxpobs = lp1
    ap2obs = 0.0E0
    ap1obs = 1.0E0
    temploop1: &
    do i = 3 , ntmpod-2
      if ( tmptab(i,lp1) >= tavlobs ) exit temploop1
    end do temploop1
    if ( abs(tavlobs-tmptab(i-1,lp1)) <= abs(tavlobs-tmptab(i,lp1)) ) then
      j1 = i-1
    else
      j1 = i
    end if
    indxt_p1obs = j1
    dent1_p1 = (tmptab(j1-1,lp1)-tmptab(j1,lp1)) * &
               (tmptab(j1-1,lp1)-tmptab(j1+1,lp1))
    dent2_p1 = (tmptab(j1,lp1)-tmptab(j1-1,lp1)) * &
               (tmptab(j1,lp1)-tmptab(j1+1,lp1))
    at1_p1obs = (tavlobs-tmptab(j1,lp1)) * &
                (tavlobs-tmptab(j1+1,lp1))/dent1_p1
    at2_p1obs = (tavlobs-tmptab(j1-1,lp1)) * &
                (tavlobs-tmptab(j1+1,lp1))/dent2_p1
    adt1_p1obs = (2.0E0*tavlobs-tmptab(j1,lp1)-tmptab(j1+1,lp1))/dent1_p1
    adt2_p1obs = (2.0E0*tavlobs-tmptab(j1-1,lp1)-tmptab(j1+1,lp1))/dent2_p1
  end subroutine settabindx_irobs

  subroutine setpath_ir_vg(pobs,obsang,sunang,nsurf,nobs, &
      n1,n2,sun,umu,umu0,delphi,puser)
    !---------------------------------------------------------
    ! purpose: computes items related to the viewing geometry.
    !---------------------------------------------------------
    implicit none
    real(4) , intent(in) :: obsang , sunang , pobs
    real(4) , intent(in) , dimension(:) :: puser
    integer , intent(out) :: nobs , nsurf , n1 , n2
    logical , intent(out) :: sun
    real(4) , intent(out) :: umu , umu0 , delphi
    real(4) :: viewang
    integer :: icnt
    !---compute the surface level
    nsurf = size(puser)
    icnt = 1
    do while ( abs(pobs - puser(icnt)) > 1.0E-4 )
      icnt = icnt+1
    end do
    nobs = icnt
    !---set viewing geometry parameters
    delphi = 0.0E0
    viewang = obsang
    if ( lookup ) then
      n1 = nsurf
      n2 = nobs-1
    else
      if ( viewang > 90.0E0 ) then
        call fatal(__FILE__,__LINE__,'eia must be less than 90')
      end if
      n1 = nobs
      n2 = nsurf-1
    end if
    umu = cos(viewang*deg2rad)
    if ( sunang < 80.0E0 ) then
      sun = .true.
      umu0 = cos(sunang*deg2rad)
    else
      sun = .false.
      umu0 = 1.0E0
    end if
  end subroutine setpath_ir_vg

  subroutine setpath_ir(pref,xg,pobs,obsang,sunang,nsurf,nobs, &
      n1,n2,sun,umu,umu0,delphi)
    !---------------------------------------------------------
    ! purpose: computes items related to the viewing geometry.
    !---------------------------------------------------------
    implicit none
    real(4) , intent(in) , dimension(:) :: pref , xg
    real(4) , intent(in) :: obsang , sunang , pobs
    integer , intent(out) :: nobs , nsurf , n1 , n2
    logical , intent(out) :: sun
    real(4) , intent(out) :: umu , umu0 , delphi
    real(4) :: psurf , viewang
    integer :: i , icnt
    !---compute the surface level
    !nobs=1
    psurf = xg(ipsfcg)
    nlevloop: &
    do i = 2 , nlev-1
      if ( pref(i) >= psurf ) exit nlevloop
    end do nlevloop
    nsurf = i
    icnt = 1
    do while ( pobs > pref(icnt) )
      icnt = icnt+1
    end do
    nobs = icnt
    !---set viewing geometry parameters
    delphi = 0.0E0
    viewang = obsang
    if ( lookup ) then
      n1 = nsurf
      n2 = nobs-1
    else
      if ( viewang > 90.0E0 ) then
        call fatal(__FILE__,__LINE__,'eia must be less than 90')
      end if
      n1 = nobs
      n2 = nsurf-1
    end if
    umu = cos(viewang*deg2rad)
    if ( sunang < 80.0E0 ) then
      sun = .true.
      umu0 = cos(sunang*deg2rad)
    else
      sun = .false.
      umu0 = 1.0E0
    end if
  end subroutine setpath_ir

  subroutine fpath_ir_vg(xg,pobs,nobs,nlast,tsfc,puser)
    !---------------------------------------------------------------------
    ! purpose:  this subroutine calculates the average temperature and
    !   molecular amounts for all layers for given profiles of temperature
    !   and molecular concentrations. it also calculates the derivatives
    !   of tavl with respect to a change in the lower and upper boundary
    !   temperatures and the derivatives of the molecular amounts with
    !   respect to a change in the mixing ratios at the layer
    !   boundaries. molecular amounts are in molec./cm**2.
    !   integration assumes that t is linear in z (lnt linear in lnp)
    !   and lnq linear in lnp.
    !---------------------------------------------------------------------
    implicit none
    integer , intent(in) :: nobs , nlast
    real(4) , intent(in) :: pobs
    real(4) , intent(in) , dimension(:) :: xg
    real(4) , intent(in) , dimension(:) :: puser
    real(4) , intent(out) :: tsfc
    integer :: nsurf , l , ih2o , k , ixoff , n , nend
    real(4) :: scal , wtot , wdry , tobs
    real(4) :: wdryobs , wtotobs , dpobs

    wfix = 0.0E0
    q = 0.0E0
    w = 0.0E0
    dtu = 0.0E0
    dtl = 0.0E0
    dwqu = 0.0E0
    dwql = 0.0E0

    nsurf = nlast+1
    !------------------------------------------------------------------------
    !     compute average pressure for the layers
    !------------------------------------------------------------------------
    pavl(1:nlast) = 0.5E0*(puser(1:nlast)+puser(2:nsurf))
    ploc(1:nsurf) = puser(1:nsurf)
    nend = nlast
    !------------------------------------------------------------------------
    !     compute average temperature for the layers
    !------------------------------------------------------------------------
    do l = 1 , nend
      dp(l) = (ploc(l+1)-ploc(l))
      scal = 1.0E0/dp(l)
      call lpsum_log(ploc(l),ploc(l+1),xg(itempg+l-1),xg(itempg+l), &
                     scal,tavl(l),dtu(l),dtl(l))
    end do
    tsfc = xg(itempg+nend)
    tobs = xg(itempg+nobs-1)
    !--------------------------------------------------------------------
    ! calculate amounts for individual species and derivatives wrt mixing
    ! ratios for retrieved constituents.
    !--------------------------------------------------------------------
    scal = rair
    ih2o = imolind(1)-1
    qcor(1:nsurf) = 1.0E0/(1.0E0+xg(ih2o+1:ih2o+nsurf))
    do k = 1 , nmol
      ixoff = imolind(k)-1
      !---transform mix. ratios into mass fractions
      !--- (relative to total air mass)
      qr(1:nsurf) = xg(ixoff+1:ixoff+nsurf)*qcor(1:nsurf)
      do l = 1 , nend
        call lpsum_log(ploc(l),ploc(l+1),qr(l),qr(l+1),&
                       scal,w(k,l),dwqu(l,k),dwql(l,k))
      end do
    end do
    !---derivatives of amount wrt dry mixing ratios
    do n = 1 , nlast
      dwqu(n,1) = dwqu(n,1)*qcor(n)**2
      dwql(n,1) = dwql(n,1)*qcor(n+1)**2
      dwqu(n,2:nmol) = dwqu(n,2:nmol)*qcor(n)
      dwql(n,2:nmol) = dwql(n,2:nmol)*qcor(n+1)
    end do
    if ( nobs > 1 ) then
      dwquobs(1) = dwquobs(1)*qcor(nobs-1)**2
      dwqlobs(1) = dwqlobs(1)*qcor(nobs)**2
      dwquobs(2:nmol) = dwquobs(2:nmol)*qcor(nobs-1)
      dwqlobs(2:nmol) = dwqlobs(2:nmol)*qcor(nobs)
    end if
    !---compute fix gas amount and average layer mixing ratios
    do l = 1 , nlast
      wtot = dp(l)*rair
      wdry = wtot-w(1,l)
      wfix(l) = wdry
      q(1,l) = w(1,l)/wtot
      q(2:nmol,l) = w(2:nmol,l)/wdry
    end do
    if ( nobs > 1 ) then
      dpobs = ploc(nobs)-pobs
      wtotobs = dpobs*rair
      wdryobs = wtotobs-wobs(1)
      wfixobs = wdryobs
      qobs(1) = wobs(1)/wtotobs
      qobs(2:nmol) = wobs(2:nmol)/wdryobs
    endif
  end subroutine fpath_ir_vg

  subroutine fpath_ir(pref,xg,pobs,nobs,nlast,tsfc)
    !---------------------------------------------------------------------
    ! purpose:  this subroutine calculates the average temperature and
    !   molecular amounts for all layers for given profiles of temperature
    !   and molecular concentrations. it also calculates the derivatives
    !   of tavl with respect to a change in the lower and upper boundary
    !   temperatures and the derivatives of the molecular amounts with
    !   respect to a change in the mixing ratios at the layer
    !   boundaries. molecular amounts are in molec./cm**2.
    !   integration assumes that t is linear in z (lnt linear in lnp)
    !   and lnq linear in lnp.
    !---------------------------------------------------------------------
    implicit none
    integer , intent(in) :: nobs , nlast
    real(4) , intent(in) :: pobs
    real(4) , intent(in) , dimension(:) :: pref
    real(4) , intent(in) , dimension(:) :: xg
    real(4) , intent(out) :: tsfc
    integer :: nsurf , l , ih2o , k , ixoff , n , nend
    real(4) :: psfc , tobs
    real(4) :: scal , scal2 , qsfc , wtot , wdry
    real(4) :: qtmp
    real(4) :: dpobs = 0.0E0
    real(4) :: scalobs2 = 0.0E0
    real(4) :: scalobs
    real(4) :: wdryobs , wtotobs
    real(4) :: tobsdum , qtmpobs

    nsurf = nlast+1

    !------------------------------------------------------------------------
    !     compute average pressure for the layers
    !------------------------------------------------------------------------
    ploc(1:nsurf) = pref(1:nsurf)
    nend = nlast-1
    psfc = xg(ipsfcg)
    !------------------------------------------------------------------------
    !     compute average temperature for the layers
    !------------------------------------------------------------------------
    do l = 1 , nend
      dp(l) = (ploc(l+1)-ploc(l))
      scal = 1.0E0/dp(l)
      call lpsum_log(ploc(l),ploc(l+1),xg(itempg+l-1),xg(itempg+l), &
                     scal,tavl(l),dtu(l),dtl(l))
    end do
    tsfc = xg(itempg+nend)
    tobs = xg(itempg+nobs-1)
    !---surface layer
    dp(nlast) = (psfc-ploc(nsurf-1))
    scal = 1.0E0/dp(nlast)
    scal2 = log(psfc/ploc(nsurf))/log(ploc(nsurf-1)/ploc(nsurf))
    call lint_log(xg,ploc,nsurf,psfc,tsfc)
    call lpsum_log(ploc(nsurf-1),psfc,xg(nsurf-1),tsfc, &
                   scal,tavl(nlast),dtu(nlast),dtl(nlast))
    dtu(nlast) = dtu(nlast)+dtl(nlast)*scal2*tsfc/xg(nsurf-1)
    dtl(nlast) = dtl(nlast)*(1.0E0-scal2)*tsfc/xg(nsurf)
    if ( nobs > 1 ) then
      dpobs = ploc(nobs)-pobs
      scalobs = 1.0E0/dpobs
      scalobs2 = log(pobs/ploc(nobs-1))/log(ploc(nobs)/ploc(nobs-1))
      call lint_log(xg,ploc,nobs,pobs,tobsdum,tobs)
      call lpsum_log(pobs,ploc(nobs),tobs,xg(nobs),&
                     scalobs,tavlobs,dtuobs,dtlobs)
      dtuobs = dtuobs*tobs/xg(nobs-1)*(1.0E0-scalobs2)
      dtlobs = dtlobs+dtuobs*scalobs2*tobs/xg(nobs)
    endif
    !--------------------------------------------------------------------
    ! calculate amounts for individual species and derivatives wrt mixing
    ! ratios for retrieved constituents.
    !--------------------------------------------------------------------
    scal = rair
    ih2o = imolind(1)-1
    qcor(1:nsurf) = 1.0E0/(1.0E0+xg(ih2o+1:ih2o+nsurf))
    do k = 1 , nmol
      ixoff = imolind(k)-1
      !---transform mix. ratios into mass fractions
      !--- (relative to total air mass)
      qr(1:nsurf) = xg(ixoff+1:ixoff+nsurf)*qcor(1:nsurf)
      do l = 1 , nend
        call lpsum_log(ploc(l),ploc(l+1),qr(l),qr(l+1),&
                       scal,w(k,l),dwqu(l,k),dwql(l,k))
      end do
      !---surface layer
      call lint_log(qr,pref,nsurf,psfc,qsfc)
      call lpsum_log(pref(nsurf-1),psfc,qr(nsurf-1),qsfc,&
                     scal,w(k,nlast),dwqu(nlast,k),dwql(nlast,k))
      dwqu(nlast,k) = dwqu(nlast,k) + dwql(nlast,k)*scal2*qsfc/qr(nsurf-1)
      dwql(nlast,k) = dwql(nlast,k)*(1.0E0-scal2)*qsfc/qr(nsurf)
      if ( nobs > 1 ) then
        call lint_log(qr,ploc,nobs,pobs,qtmpobs,qtmp)
        call lpsum_log(pobs,ploc(nobs),qtmp,qr(nobs),&
                       scal,wobs(k),dwquobs(k),dwqlobs(k))
        dwquobs(k) = dwquobs(k)*qtmp/qr(nobs-1)*(1.0E0-scalobs2)
        dwqlobs(k) = dwqlobs(k)+dwquobs(k)*scalobs2*qtmp/qr(nobs)
      end if
    end do
    !---derivatives of amount wrt dry mixing ratios
    do n = 1 , nlast
      dwqu(n,1) = dwqu(n,1)*qcor(n)**2
      dwql(n,1) = dwql(n,1)*qcor(n+1)**2
      dwqu(n,2:nmol) = dwqu(n,2:nmol)*qcor(n)
      dwql(n,2:nmol) = dwql(n,2:nmol)*qcor(n+1)
    end do
    dwquobs(1) = dwquobs(1)*qcor(nobs-1)**2
    dwqlobs(1) = dwqlobs(1)*qcor(nobs)**2
    dwquobs(2:nmol) = dwquobs(2:nmol)*qcor(nobs-1)
    dwqlobs(2:nmol) = dwqlobs(2:nmol)*qcor(nobs)
    !---compute fix gas amount and average layer mixing ratios
    do l = 1 , nlast
      wtot = dp(l)*rair
      wdry = wtot-w(1,l)
      wfix(l) = wdry
      q(1,l) = w(1,l)/wtot
      q(2:nmol,l) = w(2:nmol,l)/wdry
    end do
    if ( nobs > 1 ) then
      wtotobs = dpobs*rair
      wdryobs = wtotobs-wobs(1)
      wfixobs = wdryobs
      qobs(1) = wobs(1)/wtotobs
      qobs(2:nmol) = wobs(2:nmol)/wdryobs
    endif
  end subroutine fpath_ir

  subroutine osstran_vg(n2,nn)
    !---------------------------------------------------------------
    ! purpose:  given the luts, the function computes the layer
    !   optical depths. refer to ossrad_mw for the transmittance and
    !   brightness temperature calculation
    !---------------------------------------------------------------
    implicit none
    integer , intent(in) :: nn , n2
    integer :: indxt1_p1 , indxt2_p1 , indxt3_p1
    integer :: indxt1_p2 , indxt2_p2 , indxt3_p2
    integer :: l , ks , imol , indxp1 , indxp2
    real(4) :: dft1_p1 , dft2_p1 , dft1_p2 , dft2_p2
    real(4) :: abs0_p1 , abs0_p2 , dabs0_p1 , dabs0_p2
    real(4) :: absh2ot1_p1 , absh2ot2_p1 , absh2ot3_p1
    real(4) :: absh2o_p1 , absh2o_p2
    real(4) :: absh2ot1_p2 , absh2ot2_p2 , absh2ot3_p2
    real(4) :: dabsh2o_p1 , dabsh2o_p2
    real(4) :: dabsh2odq_p1 , dabsh2odq_p2
    real(4) :: at1_p1l , at2_p1l , adt1_p1l , adt2_p1l
    real(4) :: adt1_p2l , adt2_p2l
    real(4) :: ap1l , ap2l
 
    do l = 1 , n2
      indxt1_p1 = indxt_p1(l)-1
      indxt2_p1 = indxt1_p1+1
      indxt3_p1 = indxt1_p1+2
      indxp1 = indxp(l)
      at1_p1l = at1_p1(l)
      at2_p1l = at2_p1(l)
      adt1_p1l = adt1_p1(l)
      adt2_p1l = adt2_p1(l)
      adt1_p2l = adt1_p2(l)
      adt2_p2l = adt2_p2(l)
      ap1l = ap1(l)
      ap2l = ap2(l)

      !---fixed gases
      dft1_p1 = kfix(indxp1,indxt1_p1,nn) - kfix(indxp1,indxt3_p1,nn)
      dft2_p1 = kfix(indxp1,indxt2_p1,nn) - kfix(indxp1,indxt3_p1,nn)
      abs0_p1 = at1_p1l*dft1_p1 + at2_p1l*dft2_p1 + kfix(indxp1,indxt3_p1,nn)
      dabs0_p1 = adt1_p1l*dft1_p1 + adt2_p1l*dft2_p1
      !---water vapor
      absh2ot1_p1 = q(1,l) * dkh2o(indxp1,indxt1_p1,nn) + &
                    kh2o(indxp1,indxt1_p1,nn)
      absh2ot2_p1 = q(1,l) * dkh2o(indxp1,indxt2_p1,nn) + &
                    kh2o(indxp1,indxt2_p1,nn)
      absh2ot3_p1 = q(1,l) * dkh2o(indxp1,indxt3_p1,nn) + &
                    kh2o(indxp1,indxt3_p1,nn)
      dft1_p1 = absh2ot1_p1 - absh2ot3_p1
      dft2_p1 = absh2ot2_p1 - absh2ot3_p1
      absh2o_p1 = at1_p1l*dft1_p1 + at2_p1l*dft2_p1 + absh2ot3_p1
      dabsh2odq_p1 = at1_p1l*(dkh2o(indxp1,indxt1_p1,nn) -  &
                              dkh2o(indxp1,indxt3_p1,nn)) + &
                     at2_p1l*(dkh2o(indxp1,indxt2_p1,nn) -  &
                              dkh2o(indxp1,indxt3_p1,nn)) + &
                     dkh2o(indxp1,indxt3_p1,nn)
      dabsh2o_p1 = adt1_p1l*dft1_p1 + adt2_p1l*dft2_p1
      indxt1_p2 = indxt_p2(l)-1
      indxt2_p2 = indxt1_p2+1
      indxt3_p2 = indxt1_p2+2
      indxp2 = indxp1+1
      !---fixed gases
      dft1_p2 = kfix(indxp2,indxt1_p2,nn) - kfix(indxp2,indxt3_p2,nn)
      dft2_p2 = kfix(indxp2,indxt2_p2,nn) - kfix(indxp2,indxt3_p2,nn)
      abs0_p2 = at1_p2(l)*dft1_p2 + at2_p2(l)*dft2_p2 + &
                kfix(indxp2,indxt3_p2,nn)
      dabs0_p2 = adt1_p2l*dft1_p2 + adt2_p2l*dft2_p2
      !---water vapor
      absh2ot1_p2 = q(1,l) * dkh2o(indxp2,indxt1_p2,nn) + &
                    kh2o(indxp2,indxt1_p2,nn)
      absh2ot2_p2 = q(1,l) * dkh2o(indxp2,indxt2_p2,nn) + &
                    kh2o(indxp2,indxt2_p2,nn)
      absh2ot3_p2 = q(1,l) * dkh2o(indxp2,indxt3_p2,nn) + &
                    kh2o(indxp2,indxt3_p2,nn)
      dft1_p2 = absh2ot1_p2 - absh2ot3_p2
      dft2_p2 = absh2ot2_p2 - absh2ot3_p2
      absh2o_p2 = at1_p2(l)*dft1_p2 + at2_p2(l)*dft2_p2 + absh2ot3_p2
      dabsh2odq_p2 = at1_p2(l)*(dkh2o(indxp2,indxt1_p2,nn) -  &
                                dkh2o(indxp2,indxt3_p2,nn)) + &
                     at2_p2(l)*(dkh2o(indxp2,indxt2_p2,nn) -  &
                                dkh2o(indxp2,indxt3_p2,nn)) + &
                     dkh2o(indxp2,indxt3_p2,nn)
      dabsh2o_p2 = adt1_p2l*dft1_p2 + adt2_p2l*dft2_p2
      abso(1,l) = -(abs0_p1*ap1l + abs0_p2*ap2l)
      tautot(l) = - abso(1,l) * wfix(l)
      dtaudtmp(l) = (dabs0_p1*ap1l+dabs0_p2*ap2l) * wfix(l)
      tautot(l) = tautot(l) + (absh2o_p1*ap1l+absh2o_p2*ap2l) * w(1,l)
      abso(1,l) = (dabsh2odq_p1*ap1l+dabsh2odq_p2*ap2l)*q(1,l) + &
                  (absh2o_p1*ap1l+absh2o_p2*ap2l) + abso(1,l)
      dtaudtmp(l) = dtaudtmp(l) + (dabsh2o_p1*ap1l+dabsh2o_p2*ap2l)* w(1,l)

      !---variable gases
      do ks = 2 , nmols(nn)
        imol = imols(ks,nn)
        dft1_p1 = kvar(ks-1,indxp1,indxt1_p1,nn) - &
                  kvar(ks-1,indxp1,indxt3_p1,nn)
        dft2_p1 = kvar(ks-1,indxp1,indxt2_p1,nn) - &
                  kvar(ks-1,indxp1,indxt3_p1,nn)
        abs0_p1 = at1_p1l*dft1_p1+at2_p1l*dft2_p1 + &
                  kvar(ks-1,indxp1,indxt3_p1,nn)
        dabs0_p1 = adt1_p1l*dft1_p1 + adt2_p1l*dft2_p1
        dft1_p2 = kvar(ks-1,indxp2,indxt1_p2,nn) - &
                  kvar(ks-1,indxp2,indxt3_p2,nn)
        dft2_p2 = kvar(ks-1,indxp2,indxt2_p2,nn) - &
                  kvar(ks-1,indxp2,indxt3_p2,nn)
        abs0_p2 = at1_p2(l)*dft1_p2+at2_p2(l)*dft2_p2 + &
                  kvar(ks-1,indxp2,indxt3_p2,nn)
        dabs0_p2= adt1_p2l*dft1_p2 + adt2_p2l*dft2_p2
        abso(ks,l) = (abs0_p1*ap1l + abs0_p2*ap2l)
        tautot(l) = tautot(l) + abso(ks,l) * w(imol,l)
        dtaudtmp(l) = dtaudtmp(l)+(dabs0_p1*ap1l+dabs0_p2*ap2l)*w(imol,l)
        abso(1,l) = abso(1,l) - abso(ks,l) * q(imol,l)
       end do
    end do
  end subroutine osstran_vg

  subroutine osstran(n1,n2,nn)
    !---------------------------------------------------------------
    ! purpose:  given the luts, the function computes the layer
    !   optical depths. refer to ossrad_mw for the transmittance and
    !   brightness temperature calculation
    !---------------------------------------------------------------
    implicit none
    integer , intent(in) :: n1 , n2 , nn
    integer :: indxt1_p1 , indxt2_p1 , indxt3_p1
    integer :: l , ks , imol , indxp1
    real(4) :: dft1_p1 , dft2_p1
    real(4) :: abs0_p1 , abs0_p2 , dabs0_p1 , dabs0_p2
    real(4) :: absh2ot1_p1 , absh2ot2_p1 , absh2ot3_p1
    real(4) :: absh2o_p1 , absh2o_p2
    real(4) :: dabsh2o_p1 , dabsh2o_p2
    real(4) :: dabsh2odq_p1 , dabsh2odq_p2
    real(4) :: at1_p1l , at2_p1l , adt1_p1l , adt2_p1l
    real(4) :: adt1_p2l , adt2_p2l
    real(4) :: ap1l , ap2l

    abs0_p2 = 0.0E0
    dabs0_p1 = 0.0E0
    dabs0_p2 = 0.0E0
    absh2o_p2 = 0.0E0
    dabsh2o_p2 = 0.0E0
    dabsh2odq_p2 = 0.0E0

    if ( n1 > 1 ) then
      indxt1_p1 = indxt_p1obs-1
      indxt2_p1 = indxt1_p1+1
      indxt3_p1 = indxt1_p1+2
      indxp1 = indxpobs
      at1_p1l = at1_p1obs
      at2_p1l = at2_p1obs
      adt1_p1l = adt1_p1obs
      adt2_p1l = adt2_p1obs
      ap1l = ap1obs
      ap2l = ap2obs
      !---fixed gases
      dft1_p1 = kfix(indxp1,indxt1_p1,nn) - kfix(indxp1,indxt3_p1,nn)
      dft2_p1 = kfix(indxp1,indxt2_p1,nn) - kfix(indxp1,indxt3_p1,nn)
      abs0_p1 = at1_p1l*dft1_p1 + at2_p1l*dft2_p1 + kfix(indxp1,indxt3_p1,nn)
      !---water vapor
      absh2ot1_p1 = qobs(1) * dkh2o(indxp1,indxt1_p1,nn) + &
                    kh2o(indxp1,indxt1_p1,nn)
      absh2ot2_p1 = qobs(1) * dkh2o(indxp1,indxt2_p1,nn) + &
                    kh2o(indxp1,indxt2_p1,nn)
      absh2ot3_p1 = qobs(1) * dkh2o(indxp1,indxt3_p1,nn) + &
                    kh2o(indxp1,indxt3_p1,nn)
      dft1_p1 = absh2ot1_p1 - absh2ot3_p1
      dft2_p1 = absh2ot2_p1 - absh2ot3_p1
      absh2o_p1 = at1_p1l*dft1_p1 + at2_p1l*dft2_p1 + absh2ot3_p1
      dabsh2odq_p1 = at1_p1l*(dkh2o(indxp1,indxt1_p1,nn) -  &
                              dkh2o(indxp1,indxt3_p1,nn)) + &
                     at2_p1l*(dkh2o(indxp1,indxt2_p1,nn) -  &
                              dkh2o(indxp1,indxt3_p1,nn)) + &
                     dkh2o(indxp1,indxt3_p1,nn)
      dabsh2o_p1 = adt1_p1l*dft1_p1 + adt2_p1l*dft2_p1
      absoobs(1) = -(abs0_p1*ap1l + abs0_p2*ap2l)
      tautotobs = -absoobs(1) * wfixobs
      dtaudtmpobs = (dabs0_p1*ap1l+dabs0_p2*ap2l) * wfixobs
      tautotobs  = tautotobs + (absh2o_p1*ap1l+absh2o_p2*ap2l) * wobs(1)
      absoobs(1) = (dabsh2odq_p1*ap1l+dabsh2odq_p2*ap2l)*qobs(1) + &
                   (absh2o_p1*ap1l+absh2o_p2*ap2l) + absoobs(1)
      dtaudtmpobs = dtaudtmpobs + (dabsh2o_p1*ap1l+dabsh2o_p2*ap2l)* wobs(1)
      !---variable gases
      do ks = 2 , nmols(nn)
        imol = imols(ks,nn)
        dft1_p1 = kvar(ks-1,indxp1,indxt1_p1,nn) - &
                  kvar(ks-1,indxp1,indxt3_p1,nn)
        dft2_p1 = kvar(ks-1,indxp1,indxt2_p1,nn) - &
                  kvar(ks-1,indxp1,indxt3_p1,nn)
        abs0_p1 = at1_p1l*dft1_p1+at2_p1l*dft2_p1 + &
                  kvar(ks-1,indxp1,indxt3_p1,nn)
        dabs0_p1 = adt1_p1l*dft1_p1 + adt2_p1l*dft2_p1
        absoobs(ks) = (abs0_p1*ap1l + abs0_p2*ap2l)
        tautotobs = tautotobs + absoobs(ks) * wobs(imol)
        dtaudtmpobs = dtaudtmpobs+(dabs0_p1*ap1l+dabs0_p2*ap2l)*wobs(imol)
        absoobs(1) = absoobs(1) - absoobs(ks) * qobs(imol)
      end do
    end if
    do l = 1 , n2
      indxt1_p1 = indxt_p1(l)-1
      indxt2_p1 = indxt1_p1+1
      indxt3_p1 = indxt1_p1+2
      indxp1 = indxp(l)
      at1_p1l = at1_p1(l)
      at2_p1l = at2_p1(l)
      adt1_p1l = adt1_p1(l)
      adt2_p1l = adt2_p1(l)
      adt1_p2l = adt1_p2(l)
      adt2_p2l = adt2_p2(l)
      ap1l = ap1(l)
      ap2l = ap2(l)
      !---fixed gases
      dft1_p1 = kfix(indxp1,indxt1_p1,nn) - kfix(indxp1,indxt3_p1,nn)
      dft2_p1 = kfix(indxp1,indxt2_p1,nn) - kfix(indxp1,indxt3_p1,nn)
      abs0_p1 = at1_p1l*dft1_p1 + at2_p1l*dft2_p1 + kfix(indxp1,indxt3_p1,nn)
      dabs0_p1 = adt1_p1l*dft1_p1 + adt2_p1l*dft2_p1
      !---water vapor
      absh2ot1_p1 = q(1,l) * dkh2o(indxp1,indxt1_p1,nn) + &
                    kh2o(indxp1,indxt1_p1,nn)
      absh2ot2_p1 = q(1,l) * dkh2o(indxp1,indxt2_p1,nn) + &
                    kh2o(indxp1,indxt2_p1,nn)
      absh2ot3_p1 = q(1,l) * dkh2o(indxp1,indxt3_p1,nn) + &
                    kh2o(indxp1,indxt3_p1,nn)
      dft1_p1 = absh2ot1_p1 - absh2ot3_p1
      dft2_p1 = absh2ot2_p1 - absh2ot3_p1
      absh2o_p1 = at1_p1l*dft1_p1 + at2_p1l*dft2_p1 + absh2ot3_p1
      dabsh2odq_p1 = at1_p1l*(dkh2o(indxp1,indxt1_p1,nn) -  &
                              dkh2o(indxp1,indxt3_p1,nn)) + &
                     at2_p1l*(dkh2o(indxp1,indxt2_p1,nn) -  &
                              dkh2o(indxp1,indxt3_p1,nn)) + &
                     dkh2o(indxp1,indxt3_p1,nn)
      dabsh2o_p1 = adt1_p1l*dft1_p1 + adt2_p1l*dft2_p1
      abso(1,l) = -(abs0_p1*ap1l + abs0_p2*ap2l)
      tautot(l) = - abso(1,l) * wfix(l)
      dtaudtmp(l) = (dabs0_p1*ap1l+dabs0_p2*ap2l) * wfix(l)
      tautot(l) = tautot(l) + (absh2o_p1*ap1l+absh2o_p2*ap2l) * w(1,l)
      abso(1,l) = (dabsh2odq_p1*ap1l+dabsh2odq_p2*ap2l)*q(1,l) + &
                  (absh2o_p1*ap1l+absh2o_p2*ap2l) + abso(1,l)
      dtaudtmp(l) = dtaudtmp(l) + (dabsh2o_p1*ap1l+dabsh2o_p2*ap2l)* w(1,l)

      !---variable gases
      do ks = 2 , nmols(nn)
        imol = imols(ks,nn)
        dft1_p1 = kvar(ks-1,indxp1,indxt1_p1,nn) - &
                  kvar(ks-1,indxp1,indxt3_p1,nn)
        dft2_p1 = kvar(ks-1,indxp1,indxt2_p1,nn) - &
                  kvar(ks-1,indxp1,indxt3_p1,nn)
        abs0_p1 = at1_p1l*dft1_p1+at2_p1l*dft2_p1 + &
                  kvar(ks-1,indxp1,indxt3_p1,nn)
        dabs0_p1 = adt1_p1l*dft1_p1 + adt2_p1l*dft2_p1
        abso(ks,l) = (abs0_p1*ap1l + abs0_p2*ap2l)
        tautot(l) = tautot(l) + abso(ks,l) * w(imol,l)
        dtaudtmp(l) = dtaudtmp(l)+(dabs0_p1*ap1l+dabs0_p2*ap2l)*w(imol,l)
        abso(1,l) = abso(1,l) - abso(ks,l) * q(imol,l)
       end do
    end do
  end subroutine osstran

  end subroutine ossdrv_ir

  subroutine set_imols(nm,iimolid,iimolind)
    implicit none
    integer , intent(in) :: nm
    integer , intent(in) , dimension(nm) :: iimolid
    integer , intent(in) , dimension(nm) , optional :: iimolind
    integer :: k , n
    if ( nm > mxhmol ) then
      call fatal(__FILE__,__LINE__,'Exceeded mxhmol size!')
    end if
    if ( iimolid(1) /= 1 ) then
      call fatal(__FILE__,__LINE__, &
        'Water vapor must be specified in molecular selection!')
    end if
    n = nm
    do k = 2 , n
      if ( iimolid(k) == 0 ) then
        n = k -1
        exit
      end if
      if ( iimolid(k) <= iimolid(k-1) ) &
        call fatal(__FILE__,__LINE__, &
          'Molecular selection of HITRAN ids must be in ascending order!')
    end do
    if ( nmol > 0 ) then
      if ( n /= nmol ) then
        write(stderr_ftn,*) 'Resetting NMOL !'
        deallocate(imolid)
        deallocate(imolind)
        nmol = n
        allocate(imolid(nmol))
        allocate(imolind(nmol))
      end if
    else
      nmol = n
      allocate(imolid(nmol))
      allocate(imolind(nmol))
    end if
    imolid(1:nmol) = iimolid(1:nmol)
    if ( present(iimolind) ) then
      imolind(1:nmol) = iimolind(1:nmol)
    else
      imolind(1:nmol) = -1
    end if
  end subroutine set_imols

  subroutine set_hitran(rpref,rtmptab)
    implicit none
    real(4) , intent(in) , dimension(:) :: rpref
    real(4) , intent(in) , dimension(:,:) :: rtmptab
    integer :: i , k
    if ( nlev > 0 ) then
      if ( size(pref) /= nlev ) then
        write(stderr_ftn,*) 'Resetting NLAYOD , NTMPOD !'
        deallocate(pref)
        deallocate(pavlref)
        deallocate(tmptab)
        nlev = size(rpref)
        nlayod = nlev - 1
        if ( nlayod > mxlay ) then
          call fatal(__FILE__,__LINE__,'Dimension exceeds fixed in mxlay')
        end if
        if ( size(rtmptab,1) /= nlayod ) then
          call fatal(__FILE__,__LINE__,'Dimension mismatch pref vs tmptab')
        end if
        ntmpod = size(rtmptab,2)
        allocate(pref(nlev))
        allocate(pavlref(nlayod))
        allocate(tmptab(ntmpod,nlayod))
      end if
    else
      nlev = size(rpref)
      nlayod = nlev - 1
      if ( nlayod > mxlay ) then
        call fatal(__FILE__,__LINE__,'Dimension exceeds fixed in mxlay')
      end if
      if ( size(rtmptab,1) /= nlayod ) then
        call fatal(__FILE__,__LINE__,'Dimension mismatch pref vs tmptab')
      end if
      ntmpod = size(rtmptab,2)
      allocate(pref(nlev))
      allocate(pavlref(nlayod))
      allocate(tmptab(ntmpod,nlayod))
    end if
    pref = rpref
    do k = 1 , nlayod
      do i = 1 , ntmpod
        tmptab(i,k) = rtmptab(k,i)
      end do
    end do
    pavlref(1:nlev-1) = 0.5E0*(pref(1:nlev-1)+pref(2:nlev))
  end subroutine set_hitran

  subroutine set_solar_irradiance(rvwvn,rsunrad)
    implicit none
    real(4) , intent(in) , dimension(:) :: rvwvn
    real(4) , intent(in) , dimension(:) :: rsunrad
    if ( nfsmp > 0 ) then
      if ( size(vwvn) /= nfsmp ) then
        write(stderr_ftn,*) 'Resetting NFSMP !'
        deallocate(vwvn)
        deallocate(sunrad)
        nfsmp = size(rvwvn)
        allocate(vwvn(nfsmp))
        allocate(sunrad(nfsmp))
      end if
    else
      nfsmp = size(rvwvn)
      allocate(vwvn(nfsmp))
      allocate(sunrad(nfsmp))
    end if
    vwvn = rvwvn
    sunrad = rsunrad
  end subroutine set_solar_irradiance

  subroutine set_hitran_absorption_coefficients(himolid,himols, &
      wvptab,rcoef,rkfix,rkh2o,rdkh2o,rkvar,innch,inichmap)
    implicit none
    integer(2) , dimension(:) , intent(in) :: himolid
    integer(2) , dimension(:,:) , intent(in) :: himols
    real(4) , dimension(:,:) , intent(in) :: wvptab
    real(4) , dimension(:,:) , intent(in) :: rcoef
    real(4) , dimension(:,:,:) , intent(in) :: rkfix
    real(4) , dimension(:,:,:) , intent(in) :: rkh2o
    real(4) , dimension(:,:,:) , intent(in) :: rdkh2o
    real(4) , dimension(:,:,:,:) , intent(in) :: rkvar
    integer(2) , dimension(:) , intent(in) :: innch
    integer , dimension(:,:) , intent(in) :: inichmap
    real(4) :: ktot
    integer , dimension(mxmols) :: map
    integer , dimension(mxmols) :: maps
    integer , dimension(mxmols) :: iflag
    integer , dimension(mxmols) :: imols_indx
    real(4) , dimension(mxmols) :: kbuf
    integer :: ismp , imol , nc , nt , k , l , kk , ks
    integer :: hinmol , iimol

    if ( nlev < 0 ) &
      call fatal(__FILE__,__LINE__, &
        'OSSTRAN: Set coefficients before setting HITRAN grid.')
    if ( nfsmp < 0 ) &
      call fatal(__FILE__,__LINE__, &
        'OSSTRAN: Set coefficients before setting solar irradiance.')
    if ( nmol < 0 ) &
      call fatal(__FILE__,__LINE__, &
      'OSSTRAN: Set coefficients before selecting molecular species.')
    if ( .not. allocated(kfix) ) then
      nchmax = size(rcoef,2)
      allocate(kfix(nlayod,ntmpod,nfsmp))
      allocate(kh2o(nlayod,ntmpod,nfsmp))
      allocate(dkh2o(nlayod,ntmpod,nfsmp))
      allocate(kvar(mxmols,nlayod,ntmpod,nfsmp))
      allocate(nmols(nfsmp))
      allocate(imols(mxmols,nfsmp))
      allocate(nch(nfsmp))
      allocate(ichmap(nchmax,nfsmp))
      allocate(coef(nchmax,nfsmp))
    else
      if ( nlayod /= size(rkfix,3) ) then
        deallocate(kfix)
        deallocate(kh2o)
        deallocate(dkh2o)
        deallocate(kvar)
        deallocate(nmols)
        deallocate(imols)
        deallocate(nch)
        deallocate(ichmap)
        deallocate(coef)
        nchmax = size(rcoef,2)
        allocate(kfix(nlayod,ntmpod,nfsmp))
        allocate(kh2o(nlayod,ntmpod,nfsmp))
        allocate(dkh2o(nlayod,ntmpod,nfsmp))
        allocate(kvar(mxmols,nlayod,ntmpod,nfsmp))
        allocate(nmols(nfsmp))
        allocate(imols(mxmols,nfsmp))
        allocate(nch(nfsmp))
        allocate(ichmap(nchmax,nfsmp))
        allocate(coef(nchmax,nfsmp))
      end if
    end if
    do ismp = 1 , nfsmp
      do nc = 1 , nchmax
        coef(nc,ismp) = rcoef(ismp,nc)
        ichmap(nc,ismp) = inichmap(ismp,nc)
      end do
    end do
    nch = innch
    do ismp = 1 , nfsmp
      do nt = 1 , ntmpod
        do l = 1 , nlayod
          kh2o(l,nt,ismp) = rkh2o(ismp,nt,l)
          dkh2o(l,nt,ismp) = rdkh2o(ismp,nt,l)
        end do
      end do
    end do
    ! Check and map the requested species
    map(:) = 0
    iflag(:) = 1
    kk = 0
    iimol = count(himolid > 0)
    do k = 1 , iimol
      if ( any(imolid(:) == himolid(k)) ) then
        kk = kk + 1
        map(k) = kk
      end if
    end do
    if ( kk /= nmol ) then
      write(stderr_ftn,*) 'Database species: ', himolid(:)
      write(stderr_ftn,*) 'Selected species: ', imolid(:)
      call fatal(__FILE__,__LINE__, &
        'Requested molecules not in the database file.')
    end if
    do ismp = 1 , nfsmp
      hinmol = count(himols(ismp,:)>0)
      do nt = 1 , ntmpod
        do l = 1 , nlayod
          kfix(l,nt,ismp) = rkfix(ismp,nt,l) * wvptab(l,1)
        end do
      end do
      iflag(2:iimol) = 0
      do nt = 1 , ntmpod
        do l = 1 , nlayod
          ktot = rkfix(ismp,nt,l)
          do ks = 2 , hinmol
            kbuf(ks-1) = rkvar(ismp,nt,l,ks-1) * wvptab(l,himols(ismp,ks)+2)
            ktot = ktot + kbuf(ks-1)
          end do
          do ks = 2 , hinmol
            if ( kbuf(ks-1) > (odfac * ktot) ) iflag(himols(ismp,ks)) = 1
          end do
        end do
      end do
      maps(1:iimol) = map(1:iimol) * iflag(1:iimol)
      kk = 0
      do ks = 1 , hinmol
        imol = himols(ismp,ks)
        if ( maps(imol) > 0 ) then
          kk = kk+1
          imols_indx(kk) = ks
        else
          do nt = 1 , ntmpod
            do l = 1 , nlayod
              kfix(l,nt,ismp) = kfix(l,nt,ismp) + &
                rkvar(ismp,nt,l,ks-1) * wvptab(l,imol+2)
            end do
          end do
        end if
      end do
      nmols(ismp) = kk
      imols(1:kk,ismp) = maps(himols(ismp,imols_indx(1:kk)))
      do nt = 1 , ntmpod
        do l = 1 , nlayod
          kvar(1:kk-1,l,nt,ismp) = rkvar(ismp,nt,l,imols_indx(2:kk)-1)
        end do
      end do
    end do
  end subroutine set_hitran_absorption_coefficients

  subroutine release
    implicit none
    ntmpod = -1
    nlayod = -1
    nlev = -1
    nmol = -1
    nfsmp = -1
    nchmax = -1
    if ( allocated(imolind) ) deallocate(imolind)
    if ( allocated(imolid) ) deallocate(imolid)
    if ( allocated(pref) ) deallocate(pref)
    if ( allocated(pavlref) ) deallocate(pavlref)
    if ( allocated(tmptab) ) deallocate(tmptab)
    if ( allocated(sunrad) ) deallocate(sunrad)
    if ( allocated(vwvn) ) deallocate(vwvn)
    if ( allocated(coef) ) deallocate(coef)
    if ( allocated(kfix) ) deallocate(kfix)
    if ( allocated(kh2o) ) deallocate(kh2o)
    if ( allocated(dkh2o) ) deallocate(dkh2o)
    if ( allocated(kvar) ) deallocate(kvar)
    if ( allocated(nmols) ) deallocate(nmols)
    if ( allocated(imols) ) deallocate(imols)
    if ( allocated(nch) ) deallocate(nch)
    if ( allocated(ichmap) ) deallocate(ichmap)
  end subroutine release

end module oss_ir
