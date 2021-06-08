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
  public :: setchanselect
  public :: planck
  public :: lzsum4tot
  public :: lpsum_log
  public :: lpsum
  public :: lzsum4w
  public :: lzsum4t
  public :: settabindx_ir
  public :: getlastselectionindex
  public :: userindex2chsetindex
  public :: setindex_selfcont
  public :: vinterp
  public :: ossrad
  public :: getcountchannel
  public :: getcountusednode
  public :: threepointinterpolation
  public :: threepointinterpolationwv
  public :: cum_fix
  public :: invertmolid
  public :: shrink_var
  public :: odthresh
  public :: findfreeunit

  integer , parameter :: stderr_ftn = error_unit

  !integer , parameter :: rk8 = selected_real_kind(P=13,R=300)
  !integer , parameter :: rk4 = selected_real_kind(P= 6,R=37)

  !
  ! Dimension parameter. To allow greater values, need to recompile
  ! the code.
  !
  !P.A. changes for IASI
  !integer , parameter :: mxlev = 101
  integer, parameter  :: mxlev=205
  integer, parameter  :: mxlay=mxlev - 1
  integer , parameter :: mxhmol = 27
  integer , parameter :: mxmols = mxhmol
  integer             :: ntself
  integer, parameter, public :: UNDEFINED_INDEX = -999

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
  real(4) , parameter :: pi = 3.1415926535
  real(4) , parameter :: deg2rad = 0.017453293E0
  real(4) , parameter :: rair = 1.020408163E0
  real,    parameter :: Rgas=8.3144621   ! Universal gas constant (J / mol / K)
  real,    parameter :: AVOGAD = 6.02214199E+23 ! Avogadro ( 1 / mol )
  integer, parameter :: MAXSPC=84
  real,    parameter :: drymwt=28.964 ! Dry air molecular mass (g / mol)
  real,    parameter :: grav_const_req1= 9.80665 ! m / s^2
  real,    parameter :: grav_const_req2=-0.02586 ! m / s^2
  real,    parameter :: rEarth = 6371.23

  real,    parameter :: pMatchLim = 1.e-3
  real,    parameter :: dbAvg2dtau=-1./12!Weighting factor derivative wrt optical depth, for linear-in-tau, in low-OD limit

  real, dimension(MAXSPC), parameter  :: molWt=(/ 18.015, 44.010, 47.998,44.010, 28.011, 16.043, 31.999, &   ! 7
						  30.010, 64.060, 46.010,17.030, 63.010, 17.000, 20.010, &   ! 14
						  36.460, 80.920,127.910,51.450, 60.080, 30.030, 52.460, &   ! 21
						  28.014, 27.030, 50.490,34.010, 26.030, 30.070, 34.000, &   ! 28
						  66.010,146.050, 34.080,46.030,  0.   ,  0.   ,  0.   , &   ! 35
						   0.   ,  0.   , 28.053,32.042,  0.   ,  0.   ,  0.   , &   ! 42
						   0.   ,  0.   , 0.0   ,0.0   ,  0.   ,  0.   ,  0.   , &   ! 49
						   0.   ,153.820, 88.000,97.460,137.37,187.380,120.910 , &   ! 56
						   0.   , 86.470,108.010,  0.  , 68.12,121.05  ,  0.   , &   ! 63
						   0.   ,  0.   , 0.0   ,0.0   ,  0.   ,  0.   ,  0.   , &   ! 70
						   0.   ,  0.   , 0.0   ,0.0   ,  0.   ,  0.   ,  0.   , &   ! 77
						   0.   ,  0.   , 0.0   ,19.02 ,  0.   ,  0.   ,  0.   /)    ! 84
                   				 !  molecular weight normalized on dry air weight
  real, dimension(MAXSPC), parameter  :: normMolWt= molWt/drymwt 
     ! calculate WV Jacobian correctly
   real, dimension(MxLay,MxHmol)     :: mydwqu, mydwql

  ! Weighting factor derivative wrt optical depth, for linear-in-tau,
  ! in low-OD limit
  real(4) , parameter :: dvint = 10.0E0
  real(4) , parameter :: epsiln = 1.0E-05

  ! Filter out
  real                      :: odfac = 0.

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
  integer , dimension(:) , allocatable :: molid
  real(4) , dimension(:) , allocatable :: pref
  real(4) , dimension(:) , allocatable :: pavlref
  real(4) , dimension(:,:) , allocatable :: tmptab
  real(4) , dimension(:) , allocatable :: sunrad
  real(4) , dimension(:) , allocatable :: vwvn
  real(4) , dimension(:,:), allocatable :: wvpTmp
  real(4) , dimension(:,:) , allocatable :: coef
  real(4) , dimension(:,:) , allocatable :: kh2o
  real(4) , dimension(:,:) , allocatable :: dkh2o
  real(4) , dimension(:,:) , allocatable :: kvar

  integer*2, dimension(:)    , allocatable :: NmolS
  integer*2, dimension(:,:)  , allocatable :: ImolS
  integer , dimension(:) , allocatable :: nch
  integer , dimension(:,:) , allocatable :: ichmap

  logical , parameter :: lookup = .false.
  logical , parameter :: linintau = .false.

  ! 'variablegrid' specifies whether the input is on the absorption
  !  coefficient grid (i.e. normal operating mode) or on a user-specified grid.
  logical       :: variablegrid = .false.
  logical       :: is_planck_set = .false.
  logical       :: sphericalGeometryFlag   = .true.
  logical, save :: zIntegrationFlag        = .false.
  logical, save :: linInTauFlag            = .true.
  logical, parameter ::    OSSLININTAU     = .true.

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
  integer,save                         :: iChSet=UNDEFINED_INDEX  ! Current channel set
  integer, allocatable                 :: nChList_arr(:)
  integer, allocatable                 :: chanList_arr(:,:)
  
  integer :: indxt_p1obs = 0
  integer :: indxpobs
  integer, parameter :: mxParG=(mxHmol+1)*mxlev
  integer,save                         :: nf_sel,nNodes
  integer, dimension(:)  , allocatable :: nNodes_arr
  REAL, ALLOCATABLE             :: chanFreq(:)
  integer,                 allocatable :: iselS_arr(:,:)
  real, dimension(:,:,:) , allocatable :: coef_arr
  integer,save                         :: lastChSet    ! Number of loaded chan sets
  integer(kind=4)   :: numberNodes, currentSelection
  integer, dimension(:),allocatable        :: userIndex2chSetIndexMap
  integer,                 allocatable :: ichMap_arr(:,:,:)
  real,      dimension (:,:) , allocatable :: kfix_ir,kh2o_ir,kvar_ir
  real,      dimension (:,:) , allocatable :: dkh2o_ir
  real,      dimension (:,:) , allocatable :: kself
  real, allocatable             :: TmpSelf(:)
  integer, allocatable          :: nch_arr(:,:)
  integer                       :: mxIndex
  INTEGER, ALLOCATABLE          :: chanIndex(:)
  INTEGER(KIND=4)               :: nChan
  integer                       :: nchanAll
  integer, parameter :: LUT_KIND = KIND(1)
  integer                              :: chSetID
  integer, allocatable                 :: varMolID(:)
  integer, allocatable                 :: varMolIndex(:)
  integer, public                      :: nSceneMol

  CONTAINS

	  subroutine ossdrv_ir(nparg_in,nsf_in,nchan_in, &
			      xG,surfEmRfGrid,surfEmRf,obsAngle, solZenith, azAngle, obsLevel, y, xkt, & 
			      xkEmRf,  paxkEmRf, lat, zSurf, pUser,dbg)
	    implicit none
	    
	    !---Input variables
	    integer                            , intent(in)      :: nparg_in, nsf_in, nchan_in
	    real(4), dimension(nparg_in)       , intent(in)      :: xG
	    real(4), dimension(nsf_in)         , intent(in)      :: surfEmRfGrid
	    real, dimension(nsf_in,2)          , intent(in)      :: surfEmRf
	    real(4)                            , intent(in)      :: obsAngle, solZenith, azAngle
	    real                               , intent(in)      :: lat 
	    integer(4)                         , intent(inout)   :: obsLevel	    	    
	    real, dimension(:), intent(in)                       :: pUser
	    real,  intent(in)                                    :: zSurf
	    logical, intent(in)                                  :: dbg

	    !---Output variables
	    real, dimension(nchan_in)          ,   intent(inout) :: y
	    real, dimension(nparg_in,nchan_in) ,   intent(inout) :: xkt
	    real, dimension(nsf_in,nchan_in,2) ,   intent(inout) :: XkEmRf
	    real, dimension(2,nchan_in),           intent(inout) :: paxkemrf

	    !---Local variables
	    logical                          :: lookup,sun,referenceGrid,interpSfc
	    integer                          :: nsf,nparG,tempIndex, tSkinIndex, pSurfIndex
	    real, dimension(10)              :: zProf=-999
	    real,    dimension(Mxlev)        :: pLoc
	    real,    dimension(MxLay)        :: tavl,pavl,wfix,	dtu,dtl,ap1,ap2
	    real,    dimension(MxHmol,MxLay) :: q,w
	    real,    dimension(MxLay,MxHmol) :: dwqu,dwql
	    integer, dimension(MxLay)        :: indxt_p1,indxt_p2,indxp,indxt
	    real,    dimension(MxLay)        :: at1_p1,at2_p1,at1_p2,at2_p2,at1,at2
	    real,    dimension(MxLay)        :: adt1_p1,adt2_p1,adt1_p2,adt2_p2
	    real,    dimension(MxLay)        :: tautot,dtaudtmp
	    real,    dimension(MxmolS,MxLay) :: abso
	    real,    dimension(MxLay)        :: umuLay,umu0Lay, alt
	    real                             :: vEmRf(2),xkt_tmp(MxParG),xkEmRf_tmp(2)
	    real, dimension(mxLay)           :: rairLay !converts p(mb) in g/cm**2
	    
	    real                             :: umu,umu0,vn,fbeam,xx,rad
	    real                             :: coefInt,tSfc,plogu,plogl
	    integer                          :: n1,n2,ip0,n,nn,ich,i,ich0,nsurf
	    integer                          :: k,ks,iXoff
	    integer                          :: nChSel, ismp
	    integer                          :: nEnd,lambertian
	    real                             :: alpha
	    real, dimension(nfsmp)           :: f1_arr,f2_arr

	    real , parameter          :: f1 = 1.19106d4 !stands in numerator in dependence of radiance on TB, units: (mW/(m^2 ster cm^-1)) / [cm^-1)^3, scaled by 1e9
    	    real , parameter          :: f2 = 1.43879d3 !stands in exponential dependence of radiance on TB, like f2*wn*1e-2/TB
	    integer                              :: fid,jj
            CHARACTER(len=120)                   :: fn
      	    integer                              :: nFix
	    integer                              ::  nSceneMol

	    nSceneMol = size(varMolID)

	    ! check if requested variable molecules are in the LUT variable molecule list
	    do k=1,nScenemol
	       if (.not. ANY(varMolID(k)==molid(1:nMol))) then
                 print*, varMolID(k)
                 print*, molid(1:nMol)
		 print*, 'Err[oss_ir_module::GetOD]: Scene variable molecule list is inconsistent with the supplied LUT'
		 call exit(1)
	       end if
	    end do

	    !----- Old Main part - moved by P. Scaccia 
	    lastChSet = 0
	    iChSet = 0
	    nlev = size(pUser)
	    lambertian = 0
	    tempIndex  = 1
	    tSkinIndex = tempIndex + nlev
	    pSurfIndex = tSkinIndex + 1

 	    do k=1,nmol
	       varMolIndex(k)   = pSurfIndex+1+(k-1)*nlev
	    end do

	    do ismp=1,nfsmp
	       f1_arr(ismp)=f1*(vwvn(ismp)*1.d-3)**3
	       f2_arr(ismp)=f2*(vwvn(ismp)*1.d-3)
	    end do

	    !---- End Main

	    nlev = nlayOD + 1
	    nsf = nsf_in
            nChSel=nChList_arr(iChSet)
	    DO currentSelection=getLastSelectionIndex(),  0, -1
		    nNodes=nNodes_arr(iChSet)
		    
		    call setChanSelect(currentSelection)

		    nChan = getCountChannel()
		    numberNodes = getCountUsedNode()

		    nChSel = nChList_arr(iChSet )
		    nNodes = nNodes_arr(iChSet )


		    if (iChSet .eq. UNDEFINED_INDEX) then
			print*, 'Err[oss_ir_module::ossdrv]: Undefined channel selection set'
			call exit(1)
		    end if
		    
		    if (size(surfEmRf,1) < nsf) then
		       print*,'Err[oss_ir_module::ossdrv]:  Hinge point dimension of surfEmRf is too small'
		       call exit(1)
		    end if
		    if (size(XkEmRf,1) < nsf) then
		       print*,'Err[oss_ir_module::ossdrv]:  Hinge point dimension of XkEmRf is too small'
		       call exit(1)
		    end if
		    if (size(XkEmRf,2) < nchSel) then
		       print*,'Err[oss_ir_module::ossdrv]:  Channel dimension of XkEmRf is too small'
		       call exit(1)
		    end if
		    if (size(xkt,1).lt.size(xG)) then
		       print*,'Err[oss_ir_module::ossdrv]:  parameter dimension of xkt is too small'
		       call exit(1)
		    end if
		    if (size(xkt,2).lt.nchSel) then
		       print*,'Err[oss_ir_module::ossdrv]:  Channel dimension of xkt is too small'
		       call exit(1)
		    end if
		    if (size(xkt,2).lt.nchSel) then
		       print*,'Err[oss_ir_module::ossdrv]:  Channel dimension of xkt is too small'
		       call exit(1)
		    end if
		    
		    if (SIZE(imolind)> mxHmol) THEN
		       print*, 'Err[oss_ir_module::ossdrv]: varMolIndex Vector too large'
		       call exit(1)
		    end if
                    
		    !if (zIntegrationFlag) then
		    !  if ( .not. (present(zSurf) .and. present(zProf))) then
		    !     print*, 'Err[oss_ir_module::ossdrv]: zIntegration requires  zSurf and zProf to be provided'
		    !     call exit(1)
		    !  end if
		    !end if
		    !if (sphericalGeometryFlag) then
		    !  if ( .not. (present(zSurf) .and. present(zProf))) then
			! print*, 'Err[oss_ir_module::ossdrv]: Spherical  geometry requires zSurf and zProf to be provided'
			 !call exit(1)
		      !end if
		    !end if
		    
		    nparG    = VarMolIndex(nmol) + tSkinIndex - tempIndex - 1

		    !======================================================================
		    !     Initialize radiance vector and k-matrix
		    !======================================================================
		    y      = 0.0
		    xkt    = 0.0
		    xkEmRf = 0.0
            paxkemrf(1:2,1:nchan_in) = 0.0E0

		    !======================================================================
		    !     Compute path variables 
		    !======================================================================
		    !     Compute path geometry


		    call setpath_ir(xG(pSurfIndex),pLoc,obsLevel,obsAngle,solZenith,Nsurf,sun, & 
			 umu,umu0,lookup,referenceGrid,interpSfc,pUser)
		    
                    !if (interpSfc) then  
			!if (present(zProf) .AND. (.NOT. present(zSurf)) ) then
			 !  print*,'Err[oss_ir_module::ossdrv]:  zProf requires zSurf to be presented as well'
			  ! call exit(1)
			!end if
		    !end if

		    n2=nsurf-1
		    n1=obsLevel
		    alt(1:N2+1)=0.
		    !end if
		    !---Compute average temperature, pressure and integrated molecular amounts for the layers
		    call layerAverage_ir(xG,N2,referenceGrid,interpSfc,&
			     tempIndex,pSurfIndex,pLoc,pavl,&
			     plogu,plogl, umu,umuLay,umu0,umu0Lay,alt, nEnd, alpha) 

		    if (zIntegrationFlag) then
		      call fpathZ_ir(xG,N2,nEnd,interpSfc,tempIndex, pSurfIndex,&
			     varMolIndex,pLoc,tSfc,tavl,dtu,dtl,alpha,wfix,q,w,dwqu,dwql,alt)
		    else
		      call fpathP_ir(xG,lat,N2,nEnd,interpSfc,&
			     pSurfIndex,tempIndex,varMolIndex,pLoc,tSfc,tavl,dtu,dtl,alpha, &
			     wfix,q,w,dwqu,dwql,alt)
		    end if
		    
                    ! DIFFERENZE
		    ! zintegrationflag set to False
	            ! wfix leggermente diverso
		    !indxt_p1 e 2 leggermente diversi


		    !---Compute coefficients for temperature interpolation of ODs  
		    call settabindx_ir(tavl,pavl,N2,referenceGrid, &
			 indxt_p1,indxt_p2,indxp,ap1,ap2,&
			 at1_p1,at2_p1,at1_p2,at2_p2,adt1_p1,adt2_p1,adt1_p2,adt2_p2)

		    call setIndex_selfCont(tavl,N2,indxt,at1,at2)

		    !======================================================================
		    !     Loop over spectral points
		    !======================================================================
		    ip0 = 2
		    
		    if (dbg) then
			    write(6,*) 'nNodes: ', nNodes
			    write(6,*) 'nChannels: ', nchan
		    end if

		    NfsmpLoop: do n=1,nNodes
		       nn=iselS_arr(n,iChSet)
		       !---Compute molecular optical depth for all atmospheric layers
		       call OSStran(kfix_ir(1,nn),kh2o_ir(1,nn),dkh2o_ir(1,nn),kvar_ir(1,nn), &
			    kself(:,nn),indxt_p1,indxt_p2,indxp,indxt, &
			    ap1,ap2,at1_p1,at2_p1,at1_p2,at2_p2,at1,at2,  &
			    adt1_p1,adt2_p1,adt1_p2,adt2_p2,N2,referenceGrid, &
			    NmolS(nn),imolS(1,nn),pavl,tavl,wfix,q,w,tautot,abso,dtaudtmp)
   	    	     
		       !---Interpolate input surface emissivity to node wavenumber
		       vn=vWvn(nn)

		       call vinterp(surfEmRf,surfEmRfGrid,vn,vEmrf,ip0,coefInt)

		       if (.NOT.sun) then 
			  fbeam=0.
		       else
			  fbeam=sunrad(nn)
		       end if
		       lambertian = 0
		       !---Perform RT calculations !clear sky model
		       call ossrad(tautot,abso,dtaudtmp,tavl,tSfc,dtu,dtl,dwqu,dwql, &
			    tempIndex,tSkinIndex,VarMolIndex,plogu,plogl,f1_arr(nn),f2_arr(nn),&
			    nmols(nn),imols(1,nn),xG,vEmRf,fbeam,N1,N2,sun,umuLay,&
			    umu0Lay,umu0,lookup,lambertian,rad,xkt_tmp,xkEmRf_tmp)
            

		       nchLoop: do ich=1,nch_arr(n,iChSet)
			  ich0              = ichMap_arr(ich,n,iChSet)
			  y(ich0)           = y(ich0)+rad*coef_arr(ich,n,iChSet)
			  xkt(tempIndex:tempIndex+N2,ich0) = xkt(tempIndex:tempIndex+N2,ich0) + &
			      				     xkt_tmp(tempIndex:tempIndex+N2)*coef_arr(ich,n,iChSet)
			  xkt(tSkinIndex,ich0) = xkt(tSkinIndex,ich0)+xkt_tmp(tSkinIndex)*coef_arr(ich,n,iChSet)
			  
			  do k=1,NmolS(nn)
			     ks=ImolS(k,nn)
			     iXoff=varMolIndex(ks)-1
			     xkt(iXoff+1:iXoff+N2+1,ich0)=xkt(iXoff+1:iXoff+N2+1,ich0)+xkt_tmp(iXoff+1:iXoff+N2+1)*coef_arr(ich,n,iChSet)
			  end do

			  do i=1,2
			     xx                   = xkEmRf_tmp(i)*coef_arr(ich,n,iChSet)
			     xkEmRf(ip0-1,ich0,i) = xkEmRf(ip0-1,ich0,i)+ xx*(1.-coefInt)
			     xkEmRf(ip0,ich0,i)   = xkEmRf(ip0,ich0,i)  + xx*coefInt
			     paxkemrf(i,ich0) = paxkemrf(i,ich0) + xx
			  end do
	       end do nchLoop
	    end do NfsmpLoop
	  END DO
	
	  if (dbg) then
		  fn = 'jac'
		  print *, '  saving Jacobian to: ', 'jac'
		  fid = findFreeUnit()
		  open(unit=fid,file=trim(fn),status='unknown', action='write')

        	  DO n=1,nChan
		 	  write(fid,'(3000(1pe28.18E4))') chanFreq(n), (xkt(jj, n), jj=1,3*nlev + 2)
		  END DO
		  close(unit=fid)

        	  fn = 'rad'
        	  print *, '  saving radiance to: ', trim(fn)
        	  fid = findFreeUnit()
          
        	  open(unit=fid,file=trim(fn),status='unknown', action='write')
        	  DO n=1,nChan
        	      write(fid,'(6(1pe28.18E4))') chanFreq(n), y(n)
        	  END DO

		  close(unit=fid)  
	  end if
	  return
	  
	  contains	!- OSSDRV_IR subroutines	  	  	

		  !----------------------------------------------------------------------------
		  ! PURPOSE: Computes items related to the viewing geometry.
		  !----------------------------------------------------------------------------
		  subroutine setpath_ir(Psurf,pLoc,obsLevel,obsAngle,solZenith,nsurf,sun,umu,umu0,lookup,referenceGrid, &
		&      interpSfc,pUser)
		    real,                    parameter     :: SUN_HORIZON = 85.0
		    !---Input variables
		    real,                       intent(in) :: Psurf,obsAngle,solZenith
		    real, dimension(:),intent(in) :: pUser
		    !---Output variables
		    real, dimension(:),      intent(inout) :: pLoc
		    integer,                 intent(inout) :: obsLevel
		    integer,                 intent(inout) :: Nsurf
		    logical,                 intent(inout) :: sun
		    logical,                 intent(inout) :: lookup
		    real,                    intent(inout) :: umu,umu0
		    logical,                 intent(inout) :: referenceGrid,interpSfc
		    !---Local variables
		    real                 :: viewang
		    integer              :: i,npLev


		    !---Set viewing geometry and pressure grid parameters
		    referenceGrid=.false. 
		    npLev=size(pUser)
		    pLoc(1:npLev)=pUser(:)
		    if (abs(Psurf-pUser(npLev)) < (pUser(npLev)-pUser(npLev-1))* &
			      pMatchLim) then
			  interpSfc =.false.
		    else
			  interpSfc = .true.
         	    end if
		    if (interpSfc) then
		       do i=2,nplev-1
			 if (pLoc(i).GE. psurf) EXIT
		       end do
		       NSurf=i
		    else
		       NSurf=npLev
		       end if

		    if (obsLevel <= -1) then
		       obsLevel = nSurf
		    elseif (obsLevel == 0) then
		       obsLevel = 1
		    elseif (obsLevel > nSurf) then
		       print*,'Err[oss_ir_module::setpath_ir]: obsLevel is beyond of its limit'
		       call exit(1)
		    end if

		    viewang = obsAngle
		    if (viewang > 90.0) then
		       lookup=.true.
		       viewang=180.-viewang
		    else
		       lookup=.false.
		    end if
		    if (lookup .and. obsLevel == 1) THEN
		       print*, 'Err[oss_ir_module::setpath_ir]: Lookup does not work for satellite'
		       call exit(1)
		    end if 
		    umu =COS(viewang*deg2rad)         
		    if (solZenith < SUN_HORIZON)THEN
		       sun=.true.
		       umu0=COS(solZenith*deg2rad)         
		    else
		       sun=.false.
		       umu0=1.
		    end if
		    return
		  end subroutine setpath_ir

		  subroutine layerAverage_ir(xG,N2,referenceGrid,interpSfc,&
			     tempIndex,pSurfIndex,pLoc,pavl, &
			     plogu,plogl, umu,umuLay,umu0,umu0Lay,alt, nEnd, alpha)
	 
		    real, dimension(:)         , intent(in) :: xG
		    real                       , intent(in) :: umu,umu0
		    real, dimension(:)         , intent(in) :: pLoc
		    integer                    , intent(in) :: N2
		    logical                    , intent(in) :: referenceGrid,interpSfc
		    integer                    , intent(in) :: tempIndex,pSurfIndex
		    real, dimension(mxlay)     , intent(in) :: alt
		    integer                    , intent(out) :: nEnd
		    !---Output variables
		    real, dimension(:)           ,intent(inout):: pavl
		    real, dimension(:)           ,intent(inout):: umuLay,umu0Lay
		    real                         ,intent(out)  :: plogl
		    real                         ,intent(out)  :: plogu
		    real                         ,intent(out)  :: alpha
		    !---Local variables
		    integer                         :: l,iH2o,k,iXoff,n
		    real, dimension(mxlev)          :: qcor,qr
		    real                            :: wtot
		    real                            :: rair
		    real                            :: rRatioSq


		   if ( .not.referenceGrid) then
		      pavl(1:n2) = &
			 (pLoc(2:n2+1)-pLoc(1:n2))/log(pLoc(2:n2+1)/pLoc(1:n2))
		   end if

		   plogu = log(xG(pSurfIndex)/pLoc(n2))
		   plogl = log(xG(pSurfIndex)/pLoc(n2+1))
		   alpha = plogl/(plogl - plogu)
		   
		   if (interpSfc) THEN
		      nEnd = N2-1
		      pavl(n2)=(xG(pSurfIndex)-pLoc(n2))/plogu
		    else
		      nEnd = N2
		   end if
		    
		    if (sphericalGeometryFlag) then
		!Here set for "secants" is done in accordance with LBL
		       do l=1,n2
			 rRatioSq = ((rEarth + alt(n2+1))/(rEarth + 0.5 * (alt(l) + alt(l + 1))))**2
			 umuLay(l) = sqrt(1.0 - rRatioSq*(1.0 - umu**2))
			 umu0Lay(l) = sqrt(1.0 - rRatioSq*(1.0 - umu0**2))
		       end do
		    else
			!plane parallel version
			umuLay(1:N2) = umu
			umu0Lay(1:N2) = umu0
		    end if    
			     
		  end subroutine layerAverage_ir             
		!----------------------------------------------------------------------------
		! PURPOSE:  This subroutine calculates the average temperature and molecular
		!           amounts for all layers for given profiles of temperature and
		!           molecular concentrations. It also calculates the derivatives of
		!           tavl with respect to a change in the lower and upper boundary
		!           temperatures and the derivatives of the molecular amounts with
		!           respect to a change in the mixing ratios at the layer
		!           boundaries. Molecular amounts are in molec./cm**2. 
		!           Integration assumes that T is linear in z (LnT linear in LnP)
		!           and LnQ linear in LnP. 
		!           it implements integration over altitude
		!----------------------------------------------------------------------------
		 subroutine fpathZ_ir(xG,N2,nEnd, interpSfc, tempIndex, pSurfIndex,&
			     varMolIndex,pLoc,tSfc,tavl,dtu,dtl, alpha, wfix,q,w,dwqu,dwql,alt)
		   !---Input variables
		   real, dimension(:)         , intent(in) :: xG
		   real, dimension(:)         , intent(in) :: pLoc
		   integer                    , intent(in) :: nEnd
		   integer                    , intent(in) :: N2
		   logical                    , intent(in) :: interpSfc
		   integer                    , intent(in) :: tempIndex
		   integer                    , intent(in) :: pSurfIndex
		   integer, dimension(:)      , intent(in) :: varMolIndex
		   real                       , intent(in) :: alpha
		   !---Output variables
		   real                         ,intent(out)  :: tSfc
		   real, dimension(:)           ,intent(inout):: tavl,dtu,dtl
		   real, dimension(:)           ,intent(inout):: wfix
		   real, dimension(:,:)         ,intent(inout):: q,w
		   real, dimension(:,:)         ,intent(inout):: dwqu,dwql
		   !---Local variables
		   integer                         :: l,iH2o,k,iXoff,n
		   real                            :: pSfc,qsfc,dzsfc,zaux
		   real, dimension(mxlev)          :: qcor,qr
		   real                            :: RHOTOT1,RHOTOT2
		   real                            :: wtot
		   real, dimension(mxlay)          :: rRatioSq,alt
		   real :: scal
		   real, parameter  :: BOLTZ = 1.3806503E-16 
		   real, parameter :: GASCON = 8.314472E+07

		    !------------------------------------------------------------------------
		    !     Compute Average Temperature for the layers
		    !------------------------------------------------------------------------
		   pSfc = xG(pSurfIndex)
		   iXoff = tempIndex-1
		    do l=1,nEnd
		      call lzsum4T(xG(iXoff+l),xG(iXoff+l+1),tavl(l),dtu(l),dtl(l))
		    end do
		    
		    if (interpSfc) then
		       !---Surface layer		      
		      tSfc=(xG(iXoff+N2)**alpha) * (xG(iXoff+N2+1)**(1.0-alpha))
		      call lzsum4T(xG(iXoff+N2),tSfc,tavl(N2),dtu(N2),dtl(N2))
		      
		      dtu(N2)=dtu(N2)+dtl(N2)*alpha*tSfc/xG(iXoff+N2)
		      dtl(N2)=dtl(N2)*(1.0-alpha)*tSfc/xG(iXoff+N2+1)
		    else
		      tSfc=xG(tempIndex + N2)
		    end if
		   
		   iH2o=varMolIndex(1)-1

		    !------------------------------------------------------------------------
		    !     Calculate amounts for individual species and derivatives wrt mixing
		    !     ratios for retrieved constituents.
		    !------------------------------------------------------------------------
		    qcor(1:N2+1)=1./(1.+xG(iH2o+1:iH2o+N2+1)/normMolWt(molid(1))) 

		    do k=1,nmol
		       iXoff=varMolIndex(k)-1
		       !---Transform mix. ratios into mass fractions (relative to total air mass)
		       qr(1:N2+1)=xG(iXoff+1:iXoff+N2+1)/normMolWt(molid(k))*qcor(1:N2+1)
		       call lzsum4tot(pLoc(1),xG(1),RHOTOT1)
		       do l=1,nEnd
			  call lzsum4tot(pLoc(l+1),xG(l+1),RHOTOT2)
			  if (k == 1) THEN
			       call lzsum4W(qr(l), qr(l + 1), alt(l) - alt(l + 1), &
				   RHOTOT1, RHOTOT2, w(k, l), dwqu(l, k), dwql(l, k), wtot)
			       wfix(l) = wtot - w(k, l)
			       q(1, l) = w(1, l)/wtot
			    else
			       call lzsum4W(qr(l), qr(l + 1), alt(l) - alt(l + 1), &
				   RHOTOT1, RHOTOT2, w(k, l), dwqu(l, k), dwql(l, k),zProf(1))
			       q(k, l) = w(k, l)/wfix(l)
			    end if
			    RHOTOT1 = RHOTOT2
		       end do
		       if (interpSfc) then
			  !---Surface Layer
			  qsfc=(qr(N2)**alpha)*(qr(N2+1)**(1.0-alpha))
			  call lzsum4tot(pSfc,tSfc,RHOTOT2)
			  if (k == 1) THEN
			       dzsfc= alt(N2) - alt(N2+1)
			       call lzsum4W(qr(N2), qsfc, dzsfc, &
				   RHOTOT1, RHOTOT2, w(k, N2), dwqu(N2, k), dwql(N2, k), wtot)
			       wfix(N2) = wtot - w(1, N2)
			       q(1, N2) = w(1, N2)/wtot
			    else
			       dzsfc= alt(N2) - alt(N2+1)
			       call lzsum4W(qr(N2), qsfc, dzsfc, &
				   RHOTOT1, RHOTOT2, w(k, N2), dwqu(N2, k), dwql(N2, k),zProf(1))
			       q(k, N2) = w(k, N2) / wfix(N2)
			    end if
			    dwqu(N2, k) = dwqu(N2, k) + dwql(N2, k) * alpha * qsfc/qr(N2)
			    dwql(N2, k) = dwql(N2, k)*(1.0 - alpha) * qsfc/qr(N2 + 1)
		       end if
		    end do
		    !---Derivatives of amount wrt dry mixing ratios
		    do n=1,N2
		       dwqu(n,1)=dwqu(n,1)*qcor(n)**2/normMolWt(molid(1))
		       dwql(n,1)=dwql(n,1)*qcor(n+1)**2/normMolWt(molid(1))
		       dwqu(n,2:nmol)=dwqu(n,2:nmol)*qcor(n)/normMolWt(molid(2:nmol))
		       dwql(n,2:nmol)=dwql(n,2:nmol)*qcor(n+1)/normMolWt(molid(2:nmol))
		    end do
		    ! additional contribution to derivatives of amount wrt water vapor dry mixing ratio
		    do k=2, nmol
			iXoff=varMolIndex(k)-1
			do n=1,N2
			    mydwqu(n,k)=dwqu(n,k)*qcor(n)  *xG(iXoff+n)/normMolWt(molid(1))
			    mydwql(n,k)=dwql(n,k)*qcor(n+1)*xG(iXoff+n+1)/normMolWt(molid(1))
			end do
		    end do
		    return
		  end subroutine fpathZ_ir

		  !----------------------------------------------------------------------------
		  ! PURPOSE:  Given the LUTs, the function computes the layer
		  !           optical depths. Refer to ossrad_mw for the transmittance and 
		  !           brightness temperature calculation
		  !----------------------------------------------------------------------------
		  subroutine OSStran(kfix,kh2o,dkh2o,kvar,kself,indxt_p1,indxt_p2,indxp,indxt,&
		    ap1,ap2,at1_p1,at2_p1,at1_p2,at2_p2,at1,at2,adt1_p1,adt2_p1,adt1_p2,adt2_p2,&
		    N2,referenceGrid,NmolS,ImolS,pavl,tavl,wfix,q,w,tautot,abso,dtaudtmp)
		    !---Input variables
		    integer*2                        ,intent(in):: NmolS,ImolS(NmolS)
		    integer                          ,intent(in):: n2
		    integer, dimension(:)            ,intent(in):: indxt_p1,indxt_p2,indxp
		    integer, dimension(:)            ,intent(in):: indxt
		    real,    dimension(:)            ,intent(in):: at1_p1,at2_p1,at1_p2,at2_p2
		    real,    dimension(:)            ,intent(in):: adt1_p1,adt2_p1
		    real,    dimension(:)            ,intent(in):: adt1_p2,adt2_p2
		    real,    dimension(:)            ,intent(in):: wfix,ap1,ap2
		    real,    dimension(:)            ,intent(in):: at1,at2
		    real,    dimension(:)            ,intent(in):: pavl,tavl
		    real,    dimension(:,:)          ,intent(in):: q,w
		    real,    dimension(NlayOD,NtmpOD),intent(in):: kfix,kh2o,dkh2o
		    real                             ,intent(in):: kvar(NmolS-1,NlayOD,NtmpOD)
		    real                             ,intent(in):: kself(ntself)
		    logical                          ,intent(in):: referenceGrid
		    !---Output variables
		    real,    dimension(:)         ,intent(inout):: tautot,dtaudtmp
		    real,    dimension(:,:)       ,intent(inout):: abso
		    !---Local variables
		    integer                            :: indxt1_p1,indxt2_p1,indxt3_p1
		    integer                            :: indxt1_p2,indxt2_p2,indxt3_p2
		    integer                            :: l,ks,imol,indxp1,indxp2
		    integer                            :: indxt1,indxt2
		    real                               :: dft1_p1,dft2_p1,dft1_p2,dft2_p2
		    real                               :: abs0_p1,abs0_p2,dabs0_p1,dabs0_p2
		    real                               :: absh2ot1_p1,absh2ot2_p1,absh2ot3_p1
		    real                               :: absh2o_p1,absh2o_p2
		    real                               :: absh2ot1_p2,absh2ot2_p2,absh2ot3_p2
		    real                               :: dabsh2o_p1,dabsh2o_p2
		    real                               :: dabsh2odq_p1,dabsh2odq_p2    
		    real                               :: at1_p1L,at2_p1L,adt1_p1L,adt2_p1L
		    real                               :: at1_p2L,at2_p2L,adt1_p2L,adt2_p2L
		    real                               :: ap1L,ap2L,at1L,at2L
		    real                               :: tauself,dtausdtmp 
		    

		    indxp2 = 0
		    indxt1_p2 = 0
		    indxt2_p2 = 0
		    indxt3_p2 = 0
		    do l=1,N2
		       
		       indxt2_p1 = indxt_p1(l)
		       indxt1_p1 = indxt2_p1-1
		       indxt3_p1 = indxt2_p1+1
		       
		       indxp1   = indxp(l)
		       
		       indxt1   =indxt(l)
		       indxt2   =indxt1+1
		       
		       at1_p1L  = at1_p1(l)
		       at2_p1L  = at2_p1(l)
		       at1_p2L  = at1_p2(l)
		       at2_p2L  = at2_p2(l)
		       adt1_p1L = adt1_p1(l)
		       adt2_p1L = adt2_p1(l)
		       adt1_p2L = adt1_p2(l)
		       adt2_p2L = adt2_p2(l)
		       ap1L     = ap1(l)
		       ap2L     = ap2(l)
		       at1L     = at1(l)
		       at2L     = at2(l)


		       !---Fixed gases
		       call threePointInterpolation(at1_p1L, at2_p1L, adt1_p1L, adt2_p1L, &
				kfix(indxp1,indxt1_p1), kfix(indxp1,indxt2_p1),kfix(indxp1,indxt3_p1), abs0_p1, dabs0_p1)
				
		       !---Water vapor
		       call threePointInterpolationWV(at1_p1L, at2_p1L, adt1_p1L, adt2_p1L, &
				kh2o(indxp1,indxt1_p1), kh2o(indxp1,indxt2_p1), kh2o(indxp1,indxt3_p1), &
				dkh2o(indxp1,indxt1_p1), dkh2o(indxp1,indxt2_p1), dkh2o(indxp1,indxt3_p1), &
				absh2o_p1, dabsh2o_p1, dabsh2odq_p1, q(1,l))
		       
		       if (referenceGrid) then
			 abs0_p2 = 0.
			 dabs0_p2 = 0.
			 absh2o_p2 = 0.
			 dabsh2odq_p2 = 0.
			 dabsh2o_p2 = 0.
		       else
			  indxt1_p2 = indxt_p2(l)-1
			  indxt2_p2 = indxt1_p2+1
			  indxt3_p2 = indxt1_p2+2
			  indxp2 = indxp1+1
			  !---Fixed gases
			  call threePointInterpolation(at1_p2L, at2_p2L, adt1_p2L, adt2_p2L, &
				kfix(indxp2,indxt1_p2), kfix(indxp2,indxt2_p2),kfix(indxp2,indxt3_p2), abs0_p2, dabs0_p2)
			  
			  !---Water vapor
			  call threePointInterpolationWV(at1_p2L, at2_p2L, adt1_p2L, adt2_p2L, &
				kh2o(indxp2,indxt1_p2), kh2o(indxp2,indxt2_p2), kh2o(indxp2,indxt3_p2), &
				dkh2o(indxp2,indxt1_p2), dkh2o(indxp2,indxt2_p2), dkh2o(indxp2,indxt3_p2), &
				absh2o_p2, dabsh2o_p2, dabsh2odq_p2, q(1,l))
		       end if	

           	       !---Fixed gases
		       abso(1,l) = -(abs0_p1*ap1L + abs0_p2*ap2L)
		       tautot(l) = - abso(1,l) * wfix(l)

		       dtaudtmp(l) = (dabs0_p1*ap1L+dabs0_p2*ap2L) * wfix(l)
		       tauself   = (kself(indxt1)*at1L+kself(indxt2)*at2L)*q(1,l)
		       dtausdtmp = (kself(indxt2)-kself(indxt1))/(tmpself(indxt2)-tmpself(indxt1))*w(1,l)*q(1,l)

		       if (referenceGrid) then
			    tauself = tauself * pavlref(l)
			    dtausdtmp = dtausdtmp * pavlref(l)
		       else
			    tauself = tauself * pavl(l)
			    dtausdtmp = dtausdtmp * pavl(l)
		       end if

		       tautot(l)  = tautot(l) + (absh2o_p1*ap1L+absh2o_p2*ap2L) * w(1,l) &
			    + tauself * w(1,l)

		       abso(1,l)   = (dabsh2odq_p1*ap1L+dabsh2odq_p2*ap2L)*q(1,l) + &
			    (absh2o_p1*ap1L+absh2o_p2*ap2L) + abso(1,l)
			    
		       abso(1,l)   = abso(1,l) + tauself *2.0
		       dtaudtmp(l) = dtaudtmp(l) + (dabsh2o_p1*ap1L+dabsh2o_p2*ap2L)* w(1,l) + dtausdtmp

		       !---Variable gases
		       do kS=2,nmolS

			  Imol=ImolS(kS)
                          call threePointInterpolation(at1_p1L, at2_p1L, adt1_p1L, adt2_p1L, &
				  kvar(kS-1,indxp1,indxt1_p1), kvar(kS-1,indxp1,indxt2_p1),kvar(kS-1,indxp1,indxt3_p1), abs0_p1, dabs0_p1)

			  if (.not.referenceGrid) then
			       call threePointInterpolation(at1_p2L, at2_p2L, adt1_p2L, adt2_p2L, &
				  kvar(kS-1,indxp2,indxt1_p2), kvar(kS-1,indxp2,indxt2_p2), kvar(kS-1,indxp2,indxt3_p2), abs0_p2, dabs0_p2)

			  end if

			  abso(kS,l)   = (abs0_p1*ap1L + abs0_p2*ap2L)
			  tautot(l)   = tautot(l) + abso(kS,l) * w(Imol,l)

			  dtaudtmp(l) = dtaudtmp(l)+(dabs0_p1*ap1L+dabs0_p2*ap2L)*w(Imol,l)

		       end do
		    end do
		    return
		  end subroutine OSStran

		!----------------------------------------------------------------------------
		  ! PURPOSE:  This subroutine calculates the average temperature and molecular
		  !           amounts for all layers for given profiles of temperature and
		  !           molecular concentrations. It also calculates the derivatives of
		  !           tavl with respect to a change in the lower and upper boundary
		  !           temperatures and the derivatives of the molecular amounts with
		  !           respect to a change in the mixing ratios at the layer
		  !           boundaries. Molecular amounts are in molec./cm**2. 
		  !           Integration assumes that T is linear in z (LnT linear in LnP)
		  !           and LnQ linear in LnP. 
		  !           it implements integration over P
		  !----------------------------------------------------------------------------
		  subroutine fpathP_ir(xG,lat,N2,nEnd,interpSfc,&
			     pSurfIndex, tempIndex, varMolIndex,pLoc, tSfc,tavl,dtu,dtl, alpha, &
			     wfix,q,w,dwqu,dwql, alt)
		    !---Input variables
		    real, dimension(:)         , intent(in) :: xG
		    real                       , intent(in) :: lat
		    real, dimension(:)         , intent(in) :: pLoc
		    integer                    , intent(in) :: nEnd
		    integer                    , intent(in) :: N2
		    integer                    , intent(in) :: pSurfIndex
		    integer                    , intent(in) :: tempIndex
		    logical                    , intent(in) :: interpSfc
		    integer, dimension(:)      , intent(in) :: varMolIndex
		    real                       , intent(in) :: alpha
		    !---Output variables
		    real                         ,intent(out)  :: tSfc
		    real, dimension(:)           ,intent(inout):: tavl,dtu,dtl
		    real, dimension(:)           ,intent(inout):: wfix
		    real, dimension(:,:)         ,intent(inout):: q,w
		    real, dimension(:,:)         ,intent(inout):: dwqu,dwql
		    !---Local variables
		    integer                         :: l,iH2o,k,iXoff,n
		    real                            :: qsfc
		    real, dimension(mxlev)          :: qcor,qr
		    real, dimension(mxlay)          :: rairLay !converts p(mb) in g/cm**2
		    real                            :: wtot
		    real, dimension(mxlay)          :: alt
		    real                            :: rair
		    real :: scal, dxu, dxl, scalu, scall
		    real, parameter  :: BOLTZ = 1.3806503E-16 
		    real, parameter :: GASCON = 8.314472E+07

		    ! - Init q and w to zero
		    do k=1,mxhmol
			do l=1,mxlev - 1
				q(k,l) = 0.0
				w(k,l) = 0.0
				dwqu(l,k) = 0.0
				dwql(l,k) = 0.0
			end do
		    end do		    

		    !------------------------------------------------------------------------
		    !     Compute Average Temperature for the layers
		    !------------------------------------------------------------------------
		    iXoff = tempIndex - 1
		    
		    do l=1,nEnd
		      scal = 1./(pLoc(l+1) - pLoc(l))
		      call lpsum_log(pLoc(l), pLoc(l+1), xG(iXoff+l),xG(iXoff+l+1), &
			     scal, tavl(l),dtu(l),dtl(l))
		    end do
		 
		    if (interpSfc) then
		      !---Surface layer
		      tSfc=(xG(iXoff+N2)**alpha) * (xG(iXoff+N2+1)**(1.0-alpha))
		      scal = 1./(xG(pSurfIndex) - pLoc(N2))
		      
		      call lpsum_log(pLoc(N2), xG(pSurfIndex), xG(iXoff+N2),tSfc, &
			     scal, tavl(N2),dtu(N2),dtl(N2))
		      
		      dtu(N2)=dtu(N2)+dtl(N2)*alpha*tSfc/xG(iXoff+N2)
		      dtl(N2)=dtl(N2)*(1.0-alpha)*tSfc/xG(iXoff+N2+1)
		    else
		      tSfc=xG(iXoff+N2+1)
		    end if

		    iH2o=varMolIndex(1)-1

		    rair=10./(grav_const_req1 + grav_const_req2*COS(2.0*deg2rad*LAT))
		    rairLay(1:N2+1) = rair/(1.0 + xg(iH2o+1:iH2o+N2+1))*(1. + alt(1:N2+1)/rEarth)**2    

		    !------------------------------------------------------------------------
		    !     Calculate amounts for individual species and derivatives wrt mixing
		    !     ratios for retrieved constituents.
		    !------------------------------------------------------------------------
		    scal =  AVOGAD/ drymwt
		    do k=1,nmol
		       iXoff=varMolIndex(k)-1

		       !---Transform mix. ratios into mass fractions (relative to total air mass)
		       qr(1:N2+1)=rairLay(1:N2+1)*xG(iXoff+1:iXoff+N2+1)/normMolWt(molid(k))
    		       do l=1,nEnd
			  call lpsum_log(pLoc(l), pLoc(l+1), qr(l), qr(l+1), scal, w(k, l), dwqu(l, k), dwql(l, k))
			  if (k == 1) THEN
			       call lpsum_log(pLoc(l), pLoc(l+1), rairLay(l), rairLay(l+1), scal, wtot, dxu, dxl)
			       wfix(l) = wtot - w(k, l)
			       q(1, l) = w(1, l)/wtot
			   else
			       q(k, l) = w(k, l)/wfix(l)
			   end if
		       end do

		       if (interpSfc) then
			  !---Surface Layer
			  qsfc=(qr(N2)**alpha)*(qr(N2+1)**(1.0-alpha))
			  call lpsum_log(pLoc(N2), xG(pSurfIndex), qr(N2), qsfc, scal, w(k, l), dwqu(N2, k), dwql(N2, k))
			  if (k == 1) THEN
			       scall=(rairLay(N2)**alpha)*(rairLay(N2+1)**(1.0-alpha))
			       call lpsum_log(pLoc(N2), xG(pSurfIndex), rairLay(N2), scall, scal, wtot, dxu, dxl)
			       wfix(N2) = wtot - w(1, N2)
			       q(1, N2) = w(1, N2)/wtot
			    else
			       q(k, N2) = w(k, N2) / wfix(N2)
			    end if
			    dwqu(N2, k) = dwqu(N2, k) + dwql(N2, k) * alpha * qsfc/qr(N2)
			    dwql(N2, k) = dwql(N2, k)*(1.0 - alpha) * qsfc/qr(N2 + 1)
		       end if
		    end do
		    !---Derivatives of amount wrt dry mixing ratios
		    qcor(1:N2+1)=1/(1.+xG(iH2o+1:iH2o+N2+1))
		    do n=1,N2
		       dwqu(n,1)=dwqu(n,1)*rairLay(n)*qcor(n)    /normMolWt(molid(1))
		       dwql(n,1)=dwql(n,1)*rairLay(n+1)*qcor(n+1)/normMolWt(molid(1))
		       dwqu(n,2:nmol)=dwqu(n,2:nmol)*rairLay(n)  /normMolWt(molid(2:nmol))
		       dwql(n,2:nmol)=dwql(n,2:nmol)*rairLay(n+1)/normMolWt(molid(2:nmol))
		    end do
		    ! calculate WV Jacobian correctly
		    do k=2, nmol
			iXoff=varMolIndex(k)-1
			do n=1,N2
			    mydwqu(n,k)=dwqu(n,k)*qcor(n)  *xG(iXoff+n)/normMolWt(molID(k))
			    mydwql(n,k)=dwql(n,k)*qcor(n+1)*xG(iXoff+n+1)/normMolWt(molID(k))
			end do
		    end do

		    return
		  end subroutine fpathP_ir  

	  end subroutine ossdrv_ir


	  subroutine fatal(f,l,message)
	    implicit none
	    character(len=*) , intent(in) :: f , message
	    integer , intent(in) :: l
	    write(stderr_ftn,'(a,a,a,i8,a,a)') 'Fatal in file: ',trim(f), &
	      ' at line ',l,' : ',message
	    stop
	  end subroutine fatal
	  
	 



	  !----------------------------------------------------------------------------
	  ! PURPOSE: Computes temperature/Water vapor indexes and interpolation 
	  !          coefficients.
	  !----------------------------------------------------------------------------
	  subroutine settabindx_ir(tavl,pavl,N2,referenceGrid, &
	       indxt_p1,indxt_p2,indxp,ap1,ap2,&
	       at1_p1,at2_p1,at1_p2,at2_p2,adt1_p1,adt2_p1,adt1_p2,adt2_p2)
	    !---Input variables
	    integer                 ,intent(in) :: n2
	    real,    dimension(:)   ,intent(in) :: tavl,pavl
	    logical                 ,intent(in) :: referenceGrid
	    !---Output variables
	    integer, dimension(:),intent(inout) :: indxt_p1,indxt_p2,indxp
	    real,    dimension(:),intent(inout) :: at1_p1,at2_p1,at1_p2,at2_p2,ap1,ap2
	    real,    dimension(:),intent(inout) :: adt1_p1,adt2_p1,adt1_p2,adt2_p2
	    !---Local variables
	    integer               :: l,i,j1,j2,lp1,lp2
	    real                  :: dent1_p1,dent2_p1,dent1_p2,dent2_p2

	    lp2 = 0
	    
	    !---Compute coefficients for temperature interpolation of ODs
	    !!! INDXP: is related to the index of OD pressure layer located just above
	    !!! (by altitude) the current profile layer 
	    do l=1,N2
	       if (referenceGrid) then
		  lp1    = l
		  ap2(l) = 0.
		  ap1(l) = 1.
	       else
		  lp1=1
		  do while(pavlref(lp1).LT.pavl(l))
		     lp1=lp1+1
		
  end do
		  lp2=lp1
		  lp1=lp1-1
		  
		  if (lp1.eq.0) then
		     lp1=1
		     lp2=2
		  end if
		  if (lp2.ge.nlev) then
		     lp1=nlev-2
		     lp2=nlev-1
		  end if
		  !!! if pavl(l) becomes located between lp1 and lp2:
		  ap2(l)=(pavl(l)-pavlref(lp1))/(pavlref(lp2)-pavlref(lp1))
		  ap1(l)=(pavlref(lp2)-pavl(l))/(pavlref(lp2)-pavlref(lp1))
	       end if

               indxp(l)=lp1
	       !!! INDXT_P1: is related to the midst temperature index within the
	       !!!           upper (by altitude) OD pressure layer INDXP
	       !!! INDXT_P2: is related to the midst temperature index within the
	       !!!           lower (by altitude) OD pressure layer INDXP+1 
	       TempLoop1: do i=3,ntmpod-2
		  if(tmptab(i,lp1).GE.tavl(l)) exit TempLoop1
	       end do TempLoop1
	       if (ABS(tavl(l)-tmptab(i-1,lp1)).LE.ABS(tavl(l)-tmptab(i,lp1))) THEN
		  j1=i-1
	       else
		  j1=i
	       end if
	       indxt_p1(l)=j1

	       dent1_p1=(tmptab(j1-1,lp1)-tmptab(j1,lp1))*(tmptab(j1-1,lp1)-tmptab(j1+1,lp1))
	       dent2_p1=(tmptab(j1,lp1)-tmptab(j1-1,lp1))*(tmptab(j1,lp1)-tmptab(j1+1,lp1))
	       at1_p1(l)=(tavl(l)-tmptab(j1,lp1))*(tavl(l)-tmptab(j1+1,lp1))/dent1_p1
	       at2_p1(l)=(tavl(l)-tmptab(j1-1,lp1))*(tavl(l)-tmptab(j1+1,lp1))/dent2_p1
	       adt1_p1(l)=(2*tavl(l)-tmptab(j1,lp1)-tmptab(j1+1,lp1))/dent1_p1
	       adt2_p1(l)=(2*tavl(l)-tmptab(j1-1,lp1)-tmptab(j1+1,lp1))/dent2_p1
	       if (.not.referenceGrid) then
		  TempLoop2: do i=3,ntmpod-2
		     if(tmptab(i,lp2).GE.tavl(l)) exit TempLoop2
		  end do TempLoop2
		  if (ABS(tavl(l)-tmptab(i-1,lp2)).LE.ABS(tavl(l)-tmptab(i,lp2))) THEN
		     j2=i-1
		  else
		     j2=i
		  end if
		  indxt_p2(l)=j2
		  !!! aty_px: interpolation coefficients, which are related to pressure
		  !!!         layer (x=1,or,2, upper, or, lower) and temperature 
		  !!!         interpolation points (y=1,or,2, lower,or, midst) 
		  dent1_p2=(tmptab(j2-1,lp2)-tmptab(j2,lp2))*(tmptab(j2-1,lp2)-tmptab(j2+1,lp2))
		  dent2_p2=(tmptab(j2,lp2)-tmptab(j2-1,lp2))*(tmptab(j2,lp2)-tmptab(j2+1,lp2))
		  at1_p2(l)=(tavl(l)-tmptab(j2,lp2))*(tavl(l)-tmptab(j2+1,lp2))/dent1_p2
		  at2_p2(l)=(tavl(l)-tmptab(j2-1,lp2))*(tavl(l)-tmptab(j2+1,lp2))/dent2_p2
		  adt1_p2(l)=(2*tavl(l)-tmptab(j2,lp2)-tmptab(j2+1,lp2))/dent1_p2
		  adt2_p2(l)=(2*tavl(l)-tmptab(j2-1,lp2)-tmptab(j2+1,lp2))/dent2_p2
	       end if	 
	       !- write(6,*) l,at1_p1(l),at1_p2(l),at2_p1(l),at2_p2(l),adt1_p1(l),adt2_p1(l),adt1_p2(l),adt2_p2(l)
	    end do

          return
	  end subroutine settabindx_ir

	  ! Helper function
	  ! getLastSelectionIndex
	  ! return the last channel selection index loaded from the selection file
	  
	  integer function getLastSelectionIndex()
	    
	    getLastSelectionIndex=lastChSet
	    return
	  end function getLastSelectionIndex
	  
	   function userIndex2chSetIndex(selectionIndex)
	      integer, intent(in) :: selectionIndex
	      integer             :: userIndex2chSetIndex
	      integer             :: kk
	      
	      if (selectionIndex .eq. 0) then
		 userIndex2chSetIndex = 0
		 return
	      else
		 do kk=1, lastChSet
		    if (userIndex2chSetIndexMap(kk) .eq. selectionIndex) then
		       userIndex2chSetIndex = userIndex2chSetIndexMap(kk)
		       return
		    end if
		 end do
	      end if
	      userIndex2chSetIndex = -1

	      return
	   end function userIndex2chSetIndex


	  subroutine setIndex_selfCont(tavl,N2,indxt,at1,at2)
	    !---Input variables
	    integer                 ,intent(in) :: n2
	    real,    dimension(:)   ,intent(in) :: tavl
	    !---Output variables
	    integer, dimension(:),intent(inout) :: indxt
	    real,    dimension(:),intent(inout) :: at1,at2
	    !---Local variables
	    integer               :: l,lt,lt1,lt2

	    do l=1,N2
	       lt=1
	       do while(tmpself(lt) < tavl(l))
		  lt=lt+1
		  if (lt > ntself) exit
	       end do
	       
	       if (lt == 1) then
		  lt1=1
		  lt2=2
	       else if (lt > ntself) then
		  lt1=ntself-1
		  lt2=ntself
	       else
		  lt1=lt-1
		  lt2=lt
	       end if
	       at1(l)=(tmpself(lt2)-tavl(l))/(tmpself(lt2)-tmpself(lt1))
	       at2(l)=1.0-at1(l)
	       indxt(l)=lt1
	    end do

	    return
	  end subroutine setIndex_selfCont
	 


	  subroutine setChanSelect(selectedChSet)
	!Input
	    integer              ,intent(in) :: selectedChSet
	!local variable
	    integer                          :: chSetToSet
	    integer                          :: ich,k, kk
	    if (lastChSet < 0) THEN
		print*, 'Wrn[oss_ir_module:: setChanSelect]: Channel selections have not been loaded ', &
		         ' Default channels selection is used.'
		return
	    end if
	    chSetToSet = userIndex2chSetIndex(selectedChSet)
	    if (chSetToSet .lt. 0) then
		print*, 'Wrn[oss_ir_module:: setChanSelect]: Invalid channel set selection.', &
		              ' Active subset index: ', userIndex2chSetIndexMap(iChSet)
		return
	    end if
	    
	    if (chSetToSet == iChSet) return
	    
	    iChSet=chSetToSet

	    nNodes=nNodes_arr(iChSet)
	    return
	  end subroutine setChanSelect

	  !----------------------------------------------------------------------------
	  ! PURPOSE: Interpolation in wavenumber
	  !----------------------------------------------------------------------------
	  subroutine vinterp(datin,gridin,vn,datout,ip0,coefInt)
	    !---Input variables
	    real, dimension(:,:), intent(in)  :: datin
	    real, dimension(:)  , intent(in)  :: gridin
	    real                , intent(in)  :: vn
	    !---Input/output variables
	    integer           ,intent(inout)  :: ip0
	    !---Output variables
	    real, dimension(:),intent(inout)  :: datout
	    real              ,intent(inout)  :: coefInt
	    !---Local variables
	    integer :: nxdim,ngin,ip,i

	    nxdim=SIZE(datin,2)
	    ngin=SIZE(gridin)
	    if (vn .le. gridin(1)) then
	      DO i=1,nxdim
		 datout(i)=datin(1,i)
	      END DO
	      return
	    else if (vn .gt. gridin(ngin)) then
	      DO i=1,nxdim
		 datout(i)=datin(ngin,i)
	      END DO
	      return
	    end if

	    do ip=ip0,ngin-1
	       if(gridin(ip).GT.vn) EXIT
	    end do 
	    ip0=ip
	    coefInt=(vn-gridin(ip0-1))/(gridin(ip0)-gridin(ip0-1))
	    do i=1,nxdim
	       datout(i)=coefInt*datin(ip0,i) + (1.-coefInt)*datin(ip0-1,i)
	    end do
	    return
	  end subroutine vinterp


	  !----------------------------------------------------------------------------
	  ! PURPOSE: Compute radiances (in mw/m2/str/cm-1) and derivatives of radiances
	  !          with respect to atmospheric and surface parameters
	  !----------------------------------------------------------------------------
	  subroutine ossrad(tautot,abso,dtaudtmp,tavl,tSfc,dtu,dtl,dwqu,dwql, &
				 tempIndex,tSkinIndex,varMolIndex,plogu,plogl,f1,f2, &
				 nmols,imols,xG,surfEmRf,bs_sun,N1,N2,sun,umuLay,umu0Lay,umu0, &
				 lookup,lambertian,rad,xkt,xkEmRf)
	    !Parameters
	    real,    parameter     :: LAMBERTIAN_REFL_SECANT=1.66
	    !---Input variables
	    logical             ,intent(in) :: sun,lookup
	    integer             ,intent(in) :: lambertian
	    integer             ,intent(in) :: N1,N2
	    real                ,intent(in) :: bs_sun,tSfc
	    integer*2           ,intent(in) :: NmolS,ImolS(NmolS)
	    real, dimension(:)  ,intent(in) :: xG,surfEmRf,tautot,dtaudtmp,tavl,dtu,dtl
	    real, dimension(:,:),intent(in) :: abso,dwqu,dwql
	    real, dimension(:)  ,intent(in) :: umuLay,umu0Lay
	    integer             ,intent(in) :: tempIndex,tSkinIndex
	    integer,dimension(:),intent(in) :: varMolIndex
	    !---Output variables:
	    real              ,intent(inout):: rad
	    real, dimension(:),intent(inout):: xkt,xkEmRf
	    !---Local variables: 
	    integer                :: N2prim
	    integer                :: l,ks,k,iXoff
	    real                   :: f1,f2
	    real                   :: em,sun_rsfc,bs,dbs,sumtau_dwn
	    real                   :: tausfc,radsun,draddrsfc,draddemis
	    real                   :: txsun,sumtau0_up,umu0
	    real                   :: draddtskn,sumtau_up,dtran,rsfc
	    real                   :: dbavgdb,dbavgdbdod
	    real                   :: odsec,rlt,plogu,plogl,dtsfcdtu,dtsfcdtl
	    real, dimension(MxLev) :: txdn,txup,bbar,dbbar,draddtmp,draddtau,drdw
	    real, dimension(MxLev) :: blev,dblev,draddtmpdw,draddtmpuw
	    real, dimension(Mxlev) :: secRefl,sec,sec0,draddtau_sun
	    real                   :: localRad

	    rsfc = 0.0
	    dTsfcdtl = 0.0
	    dTsfcdtu = 0.0
	    sumtau_dwn = 0.0
	    tausfc = 1.0
	    sumtau0_up = 0.0

	    em             = surfEmRf(1)
	    sun_rsfc       = surfEmRf(2)
	    sec(1:N2)      = 1.0/umuLay(1:N2)
	    sec0(1:N2)     = 1.0/umu0Lay(1:N2)
	    if (lambertian == 0) then
	       secRefl(1:N2)= sec(1:N2)
	    else
	       secRefl(1:N2)= LAMBERTIAN_REFL_SECANT
	    end if

	  !-----------------------------------------------------------------------      
	  !     Compute Planck function and its derivative wrt temperature
	  !-----------------------------------------------------------------------
	    if (linInTauFlag) then
	       do l=1,n2
		  call planck(f1,f2,xG(l),blev(l),dblev(l))
	       end do
	       call planck(f1,f2,tSfc,blev(n2+1),dblev(n2+1))
	       call planck(f1,f2,xG(tSkinIndex),bs,dbs)
	       dtsfcdtu=tSfc/xG(n2)*plogl/(plogl-plogu)
	       dtsfcdtl=tSfc/xG(n2+1)*plogu/(plogu-plogl)       
	    endif

	    do l=1,n2
	       call planck(f1,f2,tavl(l),bbar(l),dbbar(l))
	    end do
	    call planck(f1,f2,xG(tSkinIndex),bs,dbs)
	
	    !-----------------------------------------------------------------------
	    !     Compute transmittance profile along viewing path down to surface
	    !-----------------------------------------------------------------------
	    if ( .not.lookup) then
	       txdn(1:N1)        = 1.
	       sumtau_dwn        = 0.
	       do l=N1,N2
		  sumtau_dwn     = sumtau_dwn+tautot(l)*sec(l)
		  txdn(l+1)      = EXP(-sumtau_dwn)
	       end do
	       tausfc            = txdn(N2+1)
	    end if
	    !---Initialize radiance and derivative arrays
	    rad               = 0. 
	    radsun            = 0.
	    draddrsfc         = 0.
	    draddtau_sun(1:N2)= 0.
	    draddtmp(1:N2)    = 0. 
	    draddtau(1:N2)    = 0. 
	    draddtmpdw(1:N2+1)= 0.
	    draddtmpuw(1:N2+1)= 0.
	    draddemis         = 0.
	    draddtskn         = 0.

	    !-----------------------------------------------------------------------
	    !     1- Downwelling thermal radiance calculation: 
	    !-----------------------------------------------------------------------
	    if (lookup) then
	       txup(N1:N2+1)  = 1.0
	       sumtau_dwn     = 0.0
	       sumtau_up      = 0.0
	       sumtau0_up     = 0.0
	       N2prim         = min(N1-1,N2)       



	       do l=N2prim,1,-1
		  sumtau_up  = sumtau_up+tautot(l)*sec(l)
		  txup(l)     = EXP(-sumtau_up)
	       end do
	       
	       if (linInTauFlag) then

		  do l=1,N2prim
		     dtran       = txup(l+1)-txup(l)
		     odsec=tautot(l)*sec(l)
		     call rlin(odsec,dbavgdb,dbavgdbdod)

		     if (OSSLININTAU) then
		        rlt=blev(l+1)*(1-dbavgdb)+blev(l)*dbavgdb
		        draddtmpuw(l)=draddtmpuw(l)+dtran*dblev(l)*dbavgdb
		        draddtmpuw(l+1)=dtran*dblev(l+1)*(1.-dbavgdb)
		        draddtau(l)=(txup(l)*rlt - rad + &
		             dtran*dbavgdbdod*(blev(l)-blev(l+1)))*sec(l)
		        draddtmp(l) = 0.0     
		     else   
		        rlt=2.0*bbar(l)*dbavgdb+blev(l+1)*(1.0-2.0*dbavgdb)
		        draddtmpuw(l+1)=dtran*dblev(l+1)*(1.-2.0*dbavgdb)
		        draddtmp(l)=dtran*2.0*dbbar(l)*dbavgdb
		        draddtau(l)=(txup(l)*rlt - rad + &
		             2.0*dtran*dbavgdbdod*(bbar(l)-blev(l+1)))*sec(l)
		     end if  
		     rad=rad+dtran*rlt
		  end do
	       else
		  do l=1,N2prim
		     dtran       = txup(l+1)-txup(l)
		  !---Derivative of upwelling emission wrt to temperature and optical thickness
		     draddtau(l) = (txup(l)*bbar(l)-rad)*sec(l)
		     draddtmp(l) = dtran*dbbar(l)
		     rad         = dtran*bbar(l) + rad
		  end do
	       end if

	       do l=1,N2prim
		  draddtmp(l) = draddtmp(l) + draddtau(l)*dtaudtmp(l)
	       end do       
	    else
		! downlooking
	       if (tausfc.GT.1.e-06) THEN
		  txup(N2+1)     = tausfc
		  sumtau_up     = 0.
		  do l=N2,1,-1
		     sumtau_up  = sumtau_up+tautot(l)*secRefl(l)
		     sumtau0_up  = sumtau0_up+tautot(l)*sec0(l)
		     txup(l)     = EXP(-(sumtau_dwn+sumtau_up))
		  end do
		  if (linInTauFlag) then
		     do l=1,N2
		        dtran       = txup(l+1)-txup(l)
		        odsec=tautot(l)*secRefl(l)
		        call rlin(odsec,dbavgdb,dbavgdbdod)
		        if (OSSLININTAU) then
		            rlt=blev(l+1)+(blev(l)-blev(l+1))*dbavgdb
		            draddtmpdw(l)=draddtmpdw(l)+dtran*dblev(l)*dbavgdb
		            draddtmpdw(l+1)=dtran*dblev(l+1)*(1.-dbavgdb)
		            draddtmp(l) =0.0
		            draddtau(l)=(txup(l)*rlt-rad+dtran*dbavgdbdod* &
		                                 (blev(l)-blev(l+1)))*secRefl(l)
		        else
		            rlt=2.0*bbar(l)*dbavgdb+blev(l+1)*(1.0-2.0*dbavgdb)
		            draddtmpdw(l+1)=dtran*dblev(l+1)*(1.-2.0*dbavgdb)
		            draddtmp(l)=dtran*2.0*dbbar(l)*dbavgdb
		            draddtau(l)=(txup(l)*rlt-rad+dtran*dbavgdbdod* &
		                                2.0*(bbar(l)-blev(l+1)))*secRefl(l)
		        end if
		        rad=rad+dtran*rlt
		     end do
		     draddtmpdw(N2)=draddtmpdw(N2)+draddtmpdw(N2+1)*dtsfcdtu
		     draddtmpdw(N2+1)=draddtmpdw(N2+1)*dtsfcdtl
		  else
		     do l=1,N2
		        dtran       = txup(l+1)-txup(l)
		  !---Derivative of downwelling emission wrt to temperature and optical thickness
		        draddtau(l) = (txup(l)*bbar(l)-rad)*secRefl(l)
		        draddtmp(l) = dtran*dbbar(l)
		        rad         = rad+dtran*bbar(l)
		     end do
		  end if
	       end if

	    !-----------------------------------------------------------------------
	    !      Adjust on surface reflectivity
	    !-----------------------------------------------------------------------
	       rsfc              = (1.-em)
	       do l=1,N2
		  draddtau(l) = draddtau(l)*rsfc
		  draddtmp(l) = draddtmp(l)*rsfc
	       end do
	    !-----------------------------------------------------------------------
	    !     2- Add surface component
	    !-----------------------------------------------------------------------
	    !---Derivatives wrt emissivity and sfc skin temperature: 
	       draddemis         = tausfc*bs-rad
	       draddtskn         = em*tausfc*dbs 
	       localRad          = em*tausfc*bs
	       rad               = rad*rsfc +localRad

	    !-----------------------------------------------------------------------
	    !     3- Add solar component
	    !-----------------------------------------------------------------------
	    if(sun)THEN
	       txsun          = EXP(-(sumtau_dwn+sumtau0_up))
	       !---Derivatives wrt sfc solar reflectance
	       draddrsfc      = txsun*bs_sun*umu0     
	       localRad         = sun_rsfc*draddrsfc
	       do l=1,N2
		  draddtau(l) = draddtau(l) - localRad*sec0(l)
	       end do
	       rad            = rad + localRad
	    end if 
	    
	    !     4- Upwelling thermal radiance calculation
	    !-----------------------------------------------------------------------
	       if (linInTauFlag) then
		  do l=N2,N1,-1
		     dtran=txdn(l)-txdn(l+1)
		     odsec=sec(l)*tautot(l)
		     call rlin(odsec,dbavgdb,dbavgdbdod)
		     if (OSSLININTAU) then
		         rlt=blev(l)*(1-dbavgdb)+blev(l+1)*dbavgdb
		         draddtmpuw(l)=dtran*dblev(l)*(1-dbavgdb)
		         draddtmpuw(l+1)=draddtmpuw(l+1)+dtran*dblev(l+1)*dbavgdb
		         draddtau(l)=draddtau(l) + (txdn(l+1)*rlt-rad + &
		                   dtran*dbavgdbdod*(blev(l+1)-blev(l)))*sec(l)
		     else
		         rlt=2.0*bbar(l)*dbavgdb+blev(l)*(1.-2.*dbavgdb)
		         draddtmpuw(l)=dtran*dblev(l)*(1.-2.*dbavgdb)
		         draddtmp(l)=draddtmp(l) + dtran*2.0*dbbar(l)*dbavgdb
		         draddtau(l)=draddtau(l) + (txdn(l+1)*rlt-rad + &
		                     dtran*dbavgdbdod*2.0*(bbar(l)-blev(l)))*sec(l)
		     end if
		     
		     rad=rad+dtran*rlt
		  end do
		  draddtmpuw(N2)=draddtmpuw(N2)+draddtmpuw(N2+1)*dtsfcdtu
		  draddtmpuw(N2+1)=draddtmpuw(N2+1)*dtsfcdtl
	       else
		  do l=N2,N1,-1
		     dtran          = txdn(l)-txdn(l+1)
		     draddtau(l)    = draddtau(l) + (txdn(l+1)*bbar(l)-rad)*sec(l)
		     draddtmp(l)    = draddtmp(l) + dtran*dbbar(l)
		     rad            = rad+dtran*bbar(l)
		  end do
	       end if
	       do l=1,N2
		   draddtmp(l) = draddtmp(l) + draddtau(l)*dtaudtmp(l)
	       end do
	    end if

 	    !-----------------------------------------------------------------------
	    !     Compute level derivatives and and map to array XKT: 
	    !-----------------------------------------------------------------------
	    xkt=0.
	    !  Air Temperature 
	    xkt(tempIndex+1:tempIndex+N2-1)=draddtmp(2:N2)*dtu(2:N2)+draddtmp(1:N2-1)*dtl(1:N2-1) 
	    !---level 1 
	    xkt(tempIndex)=draddtmp(1)*dtu(1) 
	    !---bottom level (N2+1) 
	    xkt(tempIndex+N2)=draddtmp(N2)*dtl(N2)
	    if (linInTauFlag) then
	       xkt(tempIndex+1:tempIndex+N2-1)=xkt(tempIndex+1:tempIndex+N2-1)+ &
		   draddtmpuw(2:N2)+draddtmpdw(2:N2)*rsfc
	       xkt(tempIndex)=xkt(tempIndex)+draddtmpuw(1)+draddtmpdw(1)*rsfc
	       xkt(tempIndex+N2)=xkt(tempIndex+N2)+draddtmpuw(N2+1)+draddtmpdw(N2+1)*rsfc
	    end if
	    !---Molecular concentrations
	    !   WV  Jacobian correctly
	    ks = 1
		  k=ImolS(kS)  !variable gases indices ImolS(1,...,nmols)
		  iXoff=varMolIndex(k)-1 !Pointer for output xkt
		  drdw(1:N2)=draddtau(1:N2)*abso(kS,1:N2)  !straight set in abso!
		  xkt(iXoff+2:iXoff+N2)=drdw(2:N2)*dwqu(2:N2,k)+drdw(1:N2-1)*dwql(1:N2-1,k) 
		     !---level 1 
		  xkt(iXoff+1)=drdw(1)*dwqu(1,k)
		     !---bottom level (N2+1) 
		  xkt(iXoff+N2+1)=drdw(N2)*dwql(N2,k) 

	    do kS=2,nmolS   !straight set of indices (1,...,nmols) for given wn
		l=ImolS(kS)  !variable gases indices ImolS(1,...,nmols)
		drdw(1:N2)=-draddtau(1:N2)*abso(kS,1:N2)  !straight set in abso!
		xkt(iXoff+2:iXoff+N2)=xkt(iXoff+2:iXoff+N2) + &
		        drdw(2:N2)  *mydwqu(2:N2,l) + &
		        drdw(1:N2-1)*mydwql(1:N2-1,l) 
		!---level 1 
		xkt(iXoff+1)=xkt(iXoff+1)+drdw(1)*mydwqu(1,l)
		!---bottom level (N2+1) 
		xkt(iXoff+N2+1)=xkt(iXoff+N2+1)+drdw(N2)*mydwql(N2,l)
	    end do

	    do kS=2,nmolS   !straight set of indices (1,...,nmols) for given wn 
	       k=ImolS(kS)  !variable gases indices ImolS(1,...,nmols)
	       iXoff=varMolIndex(k)-1 !Pointer for output xkt

	       drdw(1:N2)=draddtau(1:N2)*abso(kS,1:N2)  !straight set in abso!

	       xkt(iXoff+2:iXoff+N2)=drdw(2:N2)*dwqu(2:N2,k)+drdw(1:N2-1)*dwql(1:N2-1,k) 
	       !---level 1 
	       xkt(iXoff+1)=drdw(1)*dwqu(1,k)
	       !---bottom level (N2+1) 
	       xkt(iXoff+N2+1)=drdw(N2)*dwql(N2,k) 

	    end do
	    
	    
	    if ( .not.lookup) then
	       !---Surface terms
	       xkt(tSkinIndex)=draddtskn       !--Tskin
	       xkEmRf(1)=draddemis
	       xkEmRf(2)=draddrsfc
	    else
	       xkt(tSkinIndex)=0.
	       xkEmRf(1:2)=0.
	    end if

	    return
	  
	  contains
				  
		  !----------------------------------------------------------------------------
		  ! PURPOSE: For linear-In-Tau approach, compute the Planck radiance weight
		  !          (hlf) of the far boundary of the layer, and its derivative (dhlf)
		  !----------------------------------------------------------------------------
		  subroutine rlin(od,hlf,dhlf)
		    !Input variables
		    real,    intent(in) :: od
		    !Output variables
		    real, intent(inout) :: hlf,dhlf
		    !Local variables
		    real    tex
		    !Local constants
		    real, parameter :: odMax=20.
		    if (od.lt.pMatchLim) then
		       hlf=.5+od*dbAvg2dtau
		       dhlf=dbAvg2dtau
		    else if (od.gt.odMax) then
		       hlf=1./od
		       dhlf=-1./od**2
		    else
		       tex=exp(od)
		       hlf=1./od-1./(tex-1.)
		       dhlf=-1./od**2+tex/(tex-1)**2
		    end if
		    return
		  end subroutine rlin
 
	  end subroutine ossrad
	  
	  ! Helper function
	  ! getCountChannel
	  ! return channel count
	  
	  integer function getCountChannel()
	    integer  ::    indexLocal
	    
	    !if (present(selectionIndex)) then
	    !  indexLocal = userIndex2chSetIndex(selectionIndex)
	    !  if (indexLocal .eq. -1) then 
	    !	  print*, 'Wrn[oss_ir_module:: getCountChannel]: Invalid channel set selection.', &
		!         ' Channel count for the active subset: ', userIndex2chSetIndexMap(iChSet)
		 ! getCountChannel=nChList_arr(iChSet)
	      !else
		!  getCountChannel=nChList_arr(indexLocal)
	      !end if
	    !else
	    getCountChannel=nChList_arr(iChSet)
	    !end if
	    
	    return
	  end function getCountChannel
	  
	  ! Helper function
	  ! getCountUsedNode
	  ! return used nodes count
	  
	  integer function getCountUsedNode()
	    
	    getCountUsedNode=nNodes
	    return
	  end function getCountUsedNode


	  !----------------------------------------------------------------------------
	  ! PURPOSE: Calculates Planck function and its derivative  wrt temperature
	  !----------------------------------------------------------------------------
	  subroutine planck(f1,f2,t,rad,draddt)
	    !---Input variables
	    real     ,intent(in)      :: t,f1,f2
	    !---Output variables
	    real  ,intent(inout)      :: rad,draddt
	    !---Local variables
	    real            :: f3
	    f3=EXP(-f2/t)
	    rad=f1/(1. - f3)
	    draddt=(f2*f3/f1)*(rad/t)**2
	    rad=rad*f3
	    return
	  end subroutine planck

	  subroutine lzsum4tot(p,t,rhotot)
	    real, parameter     :: ALOSMT=2.6867775E+19  ! (1 / cm^3 )Loschmidt number
	    real, parameter     :: PZERO=1013.25         !  standard pressure ( hPa )
	    real, parameter     :: TZERO=273.15          !  standard temperature ( K )
	    real, intent(in)    :: p,t
	    !real, optional, intent(in)    :: x
	    real, intent(inout) :: RHOTOT                ! total density in molecules/cm^3
	    !real, optional, intent(inout) :: RHOW        ! density of VW in molecules/cm^3
	    real             :: ratio

	    !Local
!	    real, parameter             :: ratio

	    ratio = drymwt/molWt(1)	    
	    RHOTOT=ALOSMT*(p/PZERO)*(TZERO/t)

	    !if (present(x)) THEN
	    !   RHOW  =RHOTOT*x*ratio/(1.0+x*ratio)
	    !end if
	    return
	  end subroutine lzsum4tot

	  subroutine lzsum4T(xu,xl,xint,dxu,dxl)
	    real,   intent(in) :: xu,xl
	    real, intent(inout):: xint,dxu,dxl
	    real               :: eps
	    eps=xl/xu-1
	    if (abs(eps) > 1.e-2) then
	       xint    =(xl-xu)/log(xl/xu)
	       dxu     =xint*(1.0/(xu-xl)+1.0/xu/log(xl/xu))
	       dxl     =xint*(1.0/(xl-xu)+1.0/xl/log(xu/xl))
	    else
	       xint    =0.5*(xl+xu)-xu*eps*eps/12
	       dxu     =0.5+eps/6.0
	       dxl     =0.5-eps/6.0
	    end if
	    return
	  end subroutine lzsum4T

	  subroutine lzsum4W(xu,xl,dz,RHOTOT1,RHOTOT2,xint,dxu,dxl,xint_tot)
	    real,   intent(in) :: xu,xl,dz
	    real,   intent(in) :: RHOTOT1     !total density (molec/cm^3) at upper level
	    real,   intent(in) :: RHOTOT2     !total density (molec/cm^3) at lower level
	    real, intent(inout):: xint        !column amount (molec/cm^2)
	    real, intent(inout):: dxu,dxl
	    real, intent(inout) :: xint_tot !dry column amount (g/cm^2)

	    real               :: RHOAIR_u ! density on the upper level (molec/cm^3)
	    real               :: RHOAIR_l ! density on the lower level (molec/cm^3)
	    real               :: HZ       ! inverse exponential
	    RHOAIR_u=RHOTOT1*xu
	    RHOAIR_l=RHOTOT2*xl

	    HZ      =DZ/log(RHOAIR_u/RHOAIR_l)
	    xint    =HZ*(RHOAIR_u-RHOAIR_l)*1.E+5
	    dxu     =xint*(RHOTOT1/(RHOAIR_u-RHOAIR_l)-1.0/log(RHOAIR_u/RHOAIR_l)/xu)
	    dxl     =xint*(-RHOTOT2/(RHOAIR_u-RHOAIR_l)+1.0/log(RHOAIR_u/RHOAIR_l)/xl)
	    if ( xint_tot /= -999) THEN
	       HZ      =DZ/log(RHOTOT1/RHOTOT2) !inverse exponential for total
	       xint_tot=(HZ*(RHOTOT1-RHOTOT2)*1.E+5)
	    end if
	    return
	  end subroutine lzsum4W

	  
	  !----------------------------------------------------------------------------
	  ! PURPOSE: Computes average layer quantities (or integrated amount) 
	  !          using a linear dependence on P.
	  !----------------------------------------------------------------------------
	  subroutine lpsum(pu,pl,xu,xl,scal,xint,dxu,dxl)
	    !---Input variables
	    real, intent(in)  :: pu,pl,xu,xl,scal
	    !---Output variables
	    real, intent(out) :: xint,dxu,dxl

	    xint = 0.5*(pl-pu)*scal
	    dxu      = xint
	    dxl      = xint
	    xint     = xint*(xu+xl)
	    return
	  end subroutine lpsum

	  !----------------------------------------------------------------------------
	  ! PURPOSE: Computes average layer quantities (or integrated amount) 
	  !          using a log-x dependence on log-p.
	  !----------------------------------------------------------------------------
	  subroutine lpsum_log(pu,pl,xu,xl,scal,xint,dxu,dxl)
	    real, parameter :: epsiln=1.e-12
	    real, intent(in)  :: pu,pl
	    real, intent(in)  :: xu,xl,scal

	    !---Output variables
	    real, intent(out) :: xint,dxu,dxl
	    !---Local variables
	    real              :: hp,x0,zeta,alza,alpha
	    hp       = log(pl/pu)
	    x0       = pl*xl*hp
	    zeta     = pu*xu/(pl*xl)
	    if(ABS(zeta-1.).GT.epsiln)THEN
	       alza  = log(zeta)
	       xint  = x0*(zeta-1.)/alza
	       alpha = zeta/(zeta-1.)-1./alza
	    else
	       xint  = x0*(1.0+(zeta-1.)/2.)
	       alpha = 0.5+(zeta-1.)/12.
	    end if
	    xint     = xint*scal
	    dxu      = xint*alpha/xu
	    dxl      = xint*(1.-alpha)/xl
	    return
	  end subroutine lpsum_log

	  
	  subroutine threePointInterpolationWV(d1, d2, dd1, dd2, z1, z2, z3, dz1, dz2, dz3, y, dy, dq, q)
	     real, intent(in)   :: d1, d2, dd1, dd2
	     real, intent(in)   :: q
	     real, intent(in)   :: z1,z2,z3
	     real, intent(in)   :: dz1,dz2,dz3
	     real, intent(out)  :: y, dy, dq
	     
	     ! local
	     real        :: df1, df2
	     real        :: f1, f2, f3
	     
	     f1 = q*dz1 + z1
	     f2 = q*dz2 + z2
	     f3 = q*dz3 + z3
	     
	      df1 = f1 - f3
	      df2 = f2 - f3
	       
	      y  = d1*df1 + d2*df2 + f3
	      dy = dd1*df1 + dd2*df2
	       
	      dq = d1*(dz1-dz3) + d2*(dz2-dz3) + dz3
	       
	      return
	  end subroutine threePointInterpolationWV
	  
	  subroutine threePointInterpolation(d1, d2, dd1, dd2, z1, z2, z3, y, dy)
	     real, intent(in)   :: d1, d2, dd1, dd2
	     real, intent(in)  :: z1,z2,z3
	     real, intent(out)  :: y, dy
	     
	     ! local
	     real        :: df1, df2
	     
	       df1 = z1 - z3
	       df2 = z2 - z3
	       y  = d1*df1 + d2*df2 + z3
	       dy = dd1*df1 + dd2*df2
	       return
	  end subroutine threePointInterpolation
         !--------------------------------------------------------------- CUSTOM PART
         ! PURPOSE: Create auxiliary array, mapping molID onto a position within 
         !          array of absorption coefficients
         ! ----------------------------------------------------------------------------
         subroutine invertMolID(nmol,invmoldim,molID,invMolID)
          !Input variables
          integer,     intent(in) :: nMol,invmoldim,molID(nMol)
          !Output variables
          integer,    intent(inout) :: invMolID(invmoldim)
          !Local variables
          integer k
          invMolID(:)=0
    
          do k=1,nmol
              invMolID(molID(k))=k
          end do
          return
          end subroutine invertMolID


          !---------------------------------------------------------------------------------------------
          ! The subroutine cum_fix moves some variable species to the fixed part
  
          subroutine cum_fix(nlayod,ntmpod,kfix0,kvar0,imol,nmols,ks,w,choice)
          !---In/Out variables
          integer  ,           intent(in) :: ks,nlayod,ntmpod
          integer  ,           intent(in) :: choice
          integer*2,           intent(in) :: NmolS,Imol
          real, dimension(:,:),intent(in) :: w
          real,                intent(in) :: kvar0(NmolS-1,NlayOD,NtmpOD)
          real,             intent(inout) :: kfix0(NlayOD,NtmpOD)
          !---Local variables
          integer    :: nt

          do nt=1,ntmpOD
             if (choice == 1) then
	             kfix0(1:NLayOD,nt) = kfix0(1:NLayOD,nt) +&
                     kvar0(ks-1,1:NLayOD,nt) * w(imol+2,1:nlayod)
             else
	             kfix0(1:NLayOD,nt) = kfix0(1:NLayOD,nt) * w(1,1:NlayOD)
	      end if	 	
          end do

          return
        end subroutine cum_fix


    !---------------------------------------------------------------------------------------------
    ! The subroutine shrink_var reduce the size of kvar 
    ! and makes a new map of molecular indices 
  
    subroutine shrink_var(nlayod,ntmpod,kvar0,imols0,nmols0,maps,imols_indx,kvar1,imols1,nmols1)
       !---In/Out variables
       integer              , intent(in) :: nlayod,ntmpod
       integer*2            , intent(in) :: nmols0,imols0(nmols0),nmols1
       integer, dimension(:), intent(in) :: maps,imols_indx
       real                 , intent(in) :: kvar0(nmols0-1,NlayOD,NtmpOD)
       real              , intent(inout) :: kvar1(nmols1-1,NlayOD,NtmpOD)
       integer*2         , intent(inout) :: imols1(nmols1)

       kvar1(1:nmols1-1,1:NLayOD,1:NtmpOD) = &
               kvar0(imols_indx(2:nmols1)-1,1:NLayOD,1:NtmpOD)
                     imols1(1:nmols1) = maps(imols0(imols_indx(1:nmols1)))
       return
  end subroutine shrink_var

!---------------------------------------------------------------------------------------------
! The subroutine odthresh finds species (water vapor excluded)
! which contribution to optical depth is less than odfac of total
! and set flag in array iflag

  subroutine odthresh(nlayod,ntmpod,kfix0,kvar0,imols,nmols,w,odfac,iflag)
    !---Input variables
    integer              , intent(in) :: nlayod,ntmpod
    integer*2,                            intent(in)   :: NmolS,ImolS(NmolS)
    real,dimension(NlayOD,NtmpOD),        intent(in)   :: kfix0
    real,dimension(NmolS-1,NlayOD,NtmpOD),intent(in)   :: kvar0
    real,                                 intent(in)   :: odfac
    real, dimension(:,:),                 intent(in)   :: w
    !---Output variables
    integer, dimension(:), intent(inout)               :: iflag
    !---Local variables
    real                                   :: kbuf(50),ktot
    integer                                :: nt,l,ks
    !     Compare optical depth of individual dry constituents to 
    !     total optical depth and flag weakly absorbing molecules
    do nt = 1,NtmpOD
       do l = 1,NlayOD
          ktot = kfix0(l,nt)
          do ks = 2,nmols
             kbuf(ks-1) = kvar0(ks-1,l,nt) * w(imols(ks)+2,l)
             ktot = ktot + kbuf(ks-1)
          end do
          do ks=2,nmols  ! could exit at first occurence of iflag = 1
             if(kbuf(ks-1) > odfac * ktot) iflag(imols(ks)) = 1
          end do
       end do
    end do
    return
  end subroutine odthresh


	  subroutine set_imols(input_nScenemol,nLutmol,scene_molid,lut_molid)
	    implicit none

	    integer , intent(in) :: input_nScenemol, nLutmol
	    integer , intent(in) ,  dimension(input_nScenemol) :: scene_molid
	    integer , intent(in) ,  dimension(nLutmol)          :: lut_molid
	    integer :: k , n

	    nSceneMol = input_nScenemol
	    if ( (nLutmol > mxhmol) .or. (nScenemol > mxhmol)) then
	      call fatal(__FILE__,__LINE__,'Exceeded mxhmol size!')
	    end if
	    if ( scene_molid(1) /= 1 ) then
	      call fatal(__FILE__,__LINE__, &
		'Water vapor must be specified in molecular selection!')
	    end if
	    n = nScenemol
	    do k = 2 , n
	      if ( scene_molid(k) == 0 ) then
		n = k -1
		exit
	      end if
	      if ( scene_molid(k) <= scene_molid(k-1) ) &
		call fatal(__FILE__,__LINE__, &
		  'Molecular selection of HITRAN ids must be in ascending order!')
	    end do
	    nmol = nLutmol
	    allocate(varMolIndex(mxhmol))
	    allocate(varMolID(nScenemol))
            allocate(molid(nmol))
	    varMolID(1:nScenemol)  = scene_molid
            molid(1:nmol)          = lut_molid(1:nmol)

	  end subroutine set_imols


           subroutine set_hitran(rpref,rtmptab,rtmpself,rkself,rcwvn,risels_arr,rnnodes_arr,&
				rcoef_arr,richmap_arr,rnchlist_arr,rnch_arr, rchlist_arr,rchanfreq)
	    implicit none

	    !- Input
	    real(4) , intent(in) , dimension(:)           :: rpref
	    real(4) , intent(in) , dimension(:,:)         :: rtmptab
	    real(4),  intent(in),  dimension(:)           :: rtmpself
	    real,     intent(in),  dimension(:,:)         :: rkself
	    real,     intent(in),  dimension(:)           :: rcWvn
	    real,     intent(in),  dimension(:,:)         :: risels_arr
	    integer,  intent(in),  dimension(:)           :: rnnodes_arr
	    real,     intent(in),   dimension(:,:,:)      :: rcoef_arr
	    integer,  intent(in),  dimension(:,:,:)       :: richMap_arr
	    integer,  intent(in),  dimension(:)           :: rnchlist_arr       
	    integer,  intent(in),  dimension(:,:)         :: rnch_arr
            integer,  intent(in), dimension(:,:)          :: rchList_arr
	    real,     intent(in),   dimension(:)          :: rchanfreq

	    !- Local Variables
	    integer :: i , k, l

	    nChan  = size(rchList_arr,2)
	    mxIndex = size(rcoef_arr,1) - 1
	    nf_sel =  size(rcoef_arr,2)
	    nchmax = size(rcoef_arr,3)

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

	    pref = rpref
	    do k = 1 , nlayod
	      do i = 1 , ntmpod
		tmptab(i,k) = rtmptab(k,i)
	      end do
	    end do
    	    pavlref(1:nlayod) = (pref(2:nlayod+1)-pref(1:nlayod))/log(pref(2:nlayod+1)/pref(1:nlayod))	

	    ntself = size(rtmpself)
	    allocate(tmpself(ntself))
	    if (nfsmp <= 0) then
		nfsmp = size(rkself,1)
	    end if

	    allocate(kself(ntself,nfsmp))
	    tmpself(1:ntself) = rtmpself(1:ntself)
	       
	    do k = 1 , ntself
		kself(k,1:nfsmp) = rkself(1:nfsmp,k)
	    end do

	    !---allocate oss parameters
	    allocate (ichMap_arr(nchmax,nf_sel, 0: mxIndex))
	    allocate (nNodes_arr(0:mxIndex))
	    allocate (iselS_arr(nf_sel,0: mxIndex))
	    allocate (nch_arr(nf_sel, 0:mxIndex))
	    allocate (coef_arr(nchmax,nf_sel,0: mxIndex ))
	    allocate (nChList_arr(0:mxIndex ))
	    allocate (chanList_arr(nChan,0: mxIndex ))
	    allocate (userIndex2chSetIndexMap(0:mxIndex ))
	    allocate (chanIndex(nChan))
	    allocate (chanFreq(nChan))

	    do k =  0, mxIndex
		do l = 1, nchan
			chanList_arr(l,k) = rchList_arr(k + 1,l)
		end do
		do i=1, nf_sel
			do l=1, nchmax
				ichmap_arr(l,i,k) = richmap_arr(k + 1,i,l)
				nch_arr(i,k)      = rnch_arr(k + 1,i)
				iselS_arr(i,k)    = risels_arr(k + 1,i)
				coef_arr(l,i,k)    = rcoef_arr(k + 1,i,l)
			end do		
		end do
	    end do

	    chanIndex(1:nChan)  = chanList_arr(1:nChan,0)
	    chanFreq(1:nChan)   = rchanfreq(1:nChan)

	    nNodes_arr(0:mxIndex)      = rnnodes_arr(1:mxIndex+1)
	    nChlist_arr(0:mxIndex)     = rnchlist_arr(1:mxIndex+1)
	    userIndex2chSetIndexMap(0:mxIndex) = -1
	    userIndex2chSetIndexMap(0)         = 0    

	  end subroutine set_hitran

	  subroutine set_solar_irradiance(rvwvn,rsunrad)
	    implicit none
	    real(4) , intent(in) , dimension(:) :: rvwvn
	    real(4) , intent(in) , dimension(:) :: rsunrad

	    if ( nfsmp > 0 ) then
	      if ( size(vwvn) /= nfsmp ) then
		!if (allocated(vwvn)) deallocate(vwvn)
		!if (allocated(vwvn)) deallocate(sunrad)
		nfsmp = size(rvwvn)
		allocate(vwvn(nfsmp))
		allocate(sunrad(nfsmp))
	      end if
	    else
	      nfsmp = size(rvwvn)
	      allocate(vwvn(nfsmp))
	      allocate(sunrad(nfsmp))
	    end if
	 
	    vwvn(1:nfsmp)   = rvwvn(1:nfsmp)
	    sunrad(1:nfsmp) = rsunrad(1:nfsmp)

	  end subroutine set_solar_irradiance

          subroutine set_hitran_absorption_coefficients(rkfix,rkh2o,rdkh2o, &
			rkvar,rimols_,rnmols_,fmtLUT,xid,fixDMR,wvptab,dflt,molProf)
	    implicit none

	    real(4) , dimension(:,:) , intent(in)   :: rkfix
	    real(4) , dimension(:,:) , intent(in)   :: rkh2o
	    real(4) , dimension(:,:) , intent(in)   :: rdkh2o
	    real(4) , dimension(:,:) , intent(in)   :: rkvar
	    integer , dimension(:,:) , intent(in)   :: rimols_
	    integer , dimension(:) ,   intent(in)   :: rnmols_

	    real(4) , dimension(:)   , intent(in)   :: fixDMR
	    real(4) , dimension(:,:) , intent(in)   :: wvptab,dflt,molProf
	    integer , intent(in)                    :: fmtLUT,xid
         
            real                                    :: dummyReal
	    integer                                 :: nsize1 , nsize2, nsize3, ismp,icnt,invmoldim, m, nSceneMol
	    integer, allocatable                    :: invMolid(:)
	    integer, allocatable                    :: mFix(:)
	    integer, dimension(MxHmol)              :: map,mapS,iflag = 1
	    integer                                 :: kk, ks, k, nFix
	    real                                    :: scale
            integer*2                               :: NmolS_tmp,ImolS_tmp(MxmolS),imol
            integer                                 :: NmolS_,ImolS_(MxmolS)
            real   , allocatable                    :: kvar_tmp(:)
            real, dimension (:,:) , allocatable     :: dkFix_ir
            integer                                 :: imols_indx(MxmolS)

	    nsize1 = size(rkfix,2)
	    nsize2 = size(rkh2o,2)

	    nSceneMol = size(varMolID)
	
            iflag(2:nmol) = 0
	    if ( nlev < 0 ) &
	      call fatal(__FILE__,__LINE__, &
		'OSSTRAN: Set coefficients before setting HITRAN grid.')
	    if ( nfsmp < 0 ) &
	      call fatal(__FILE__,__LINE__, &
		'OSSTRAN: Set coefficients before setting solar irradiance.')
	    if ( nmol < 0 ) &
	      call fatal(__FILE__,__LINE__, &
	      'OSSTRAN: Set coefficients before selecting molecular species.')
	    allocate(kfix_ir(nsize1, nfsmp))
	    allocate(kh2o_ir(nsize2, nfsmp))
	    allocate(dkh2o_ir(nsize1, nfsmp))
	    allocate(kvar_ir(nsize1*mxhmol,nfsmp))	
	    allocate (nmols(nfsmp),imols(mxmols,nfsmp))


	    do ismp = 1 , nfsmp
		  nmols(ismp)		  = nmol
		  kh2o_ir(1:nsize2,ismp)  = rkh2o(ismp,1:nsize2)
		  dkh2o_ir(1:nsize1,ismp) = rdkh2o(ismp,1:nsize1)
		  kfix_ir(1:nsize1,ismp)  = rkfix(ismp,1:nsize1)
		  !kvar_ir(1:nsize1,ismp)  = rkvar(ismp,1:nsize1)
	    end do
    
            !---Build mapping vector for user-selected variable molecules
            kk=0
            map(1:nmol)=0

            do k=1,nmol
                !---map(k) contains new relative index for selected molecule and is zero for others
                if ( ANY(varMolID(1:nSceneMol) == molid(k))) THEN
                   kk=kk+1
                   map(k)=kk
                end if
            end do
            if(kk.NE.nSceneMol)THEN
               PRINT *,'Err[oss_ir_module::GetOD]: '// &
               'This database contains only the molecules with'
               PRINT *,'the following OSS ID as variable molecules: '
               PRINT *,molid(1:nmol)
               call exit(1)
            end if

	    nchanAll = nChan
	    chSetID   = 0
            userIndex2chSetIndexMap(0) = chSetID
	    call setChanSelect(chSetID)
		
            ! Put here kfix <-------------------------

	    allocate (wvpTmp(nmol+2,nlayOD))
    	    wvpTmp = 0.0
	    wvpTmp(1,1:nlayOD)=fixDMR(1:nlayOD)
	    wvpTmp(2,1:nlayOD)=wvpTab(1:nlayOD,1)
	    wvpTmp(3,1:nlayOD)=wvpTab(1:nlayOD,1)
	    !Define wvpTmp in order to use it without changes s/routines of GetOD
	    wvpTmp(4:nmol+2,1:nlayOD)=dflt(1:nlayod,2:nmol)


	    invmoldim = molID(nmol)
	    allocate (invMolID(invmoldim))
	        !in order to count variable species interms of compressed indices, we set
         	!auxiliary array invMolID

	    ! Continue reading defProfFile
	    call invertMolID(nmol,invmoldim,molID(1:nmol),invMolID)

	    nFix = nmol - nSceneMol
	    allocate (mFix(nFix)) !variable species to be transferred into Fix
	    icnt=0
	    !set mFix as those, which are not included into varMolID

	    do k=1,nMol
	       if ( .not. ANY(varMolID(1:nSceneMol) == molid(k)) ) then
		  icnt=icnt+1
		  mFix(icnt)=molid(k)
	       end if
	    end do
            if (icnt .NE. nFix) THEN
                 print*, 'Err[oss_ir_module::GetOD]: Scene variable molecular list is inconsistent with the supplied LUT'
                 call exit(1)
            end if
	    
            do m=1,nFix
	         if (mFix(m)==81) cycle
		 do k=1,nlayOD
			  scale=1./normMolWt(mFix(m))/(pref(k+1) - pref(k))
			  call lpsum_log(pref(k), pref(k+1),molProf(k,invMolid(mfix(m))), molProf(k+1,invMolid(mfix(m))), scale, &
				   wvpTmp(2+invMolID(mFix(m)), k), dummyReal,dummyReal)
		 end do
	    end do
	
	    do ismp=1,nfsmp

     		 nmols_Tmp = rnmols_(ismp)
       		 Imols_Tmp(1:nmols_tmp)=rImols_(1:nmols_tmp,ismp)
	         nsize3 = nsize1*nmols_tmp
            	 allocate(kvar_tmp(nsize3))
		kvar_tmp(1:nsize3) = rkvar(ismp,1:nsize3)

  	        if (fmtLUT < 2) then
		      CALL cum_fix(nlayod,ntmpod,kfix_ir(1,ismp),kvar_tmp,imol,nmols_tmp,ks,wvpTmp,2)
	        end if
                if (xid == 0) THEN
	          	 call cum_fix(nlayod,ntmpod,dkfix_ir(1,ismp),kvar_tmp,imol,nmols_tmp,ks,wvpTmp,2)
	        end if
	        !     Keep only user-selected molecules as variables molecules -
	        !     merge remaining molecules with fixed gases
	        !     First apply OD threshold to identify weakly absorbing constituents
	        ! to eleminate non-contributing species
    	        ! useless since odfac is set to 0

	    	iflag(2:nmol) = 0

		!-checked: kvar_tmp, imols_tmp, nmols_tmp, kfix_ir, wvptmp
	    	call odthresh(nlayod,ntmpod,kfix_ir(1,ismp),kvar_tmp,imols_tmp,nmols_tmp,wvpTmp,odfac,iflag)
	   	maps(1:nmol) = map(1:nmol) * iflag(1:nmol)
		
		!- <---- Check 
	    	kk = 0
	    	do ks = 1,nmols_tmp
			imol = rimols_(ks,ismp)
			if(mapS(imol) > 0)THEN
				kk = kk+1
				imols_indx(kk) = ks
			else
	        		kvar_tmp(1:nsize1) = rkvar(ismp,1:nsize1)
		    		!if (ismp .eq. 1)  write(6,*) shape(kvar_tmp)
				call cum_fix(nlayod,ntmpod,kfix_ir(1,ismp),kvar_tmp,imol,nmols_tmp,ks,wvpTmp,1)
			end if
	    	end do


	    	nmols(ismp) = kk
	    	call shrink_var(nlayod,ntmpod,kvar_tmp,imols_tmp,nmols_tmp,mapS,imols_indx,&
			    kvar_ir(1,ISmp),imols(1,ISmp),nmols(ISmp))
	        if (allocated(kvar_tmp)) deallocate(kvar_tmp)
            end do
           !deallocate (wvpTmp)
	   !deallocate (invmolid)
         
	  end subroutine set_hitran_absorption_coefficients


	  subroutine release()

	    ntmpod = -1
	    nlayod = -1
	    nlev = -1
	    nmol = -1
	    nfsmp = -1
	    nchmax = -1


	    if ( allocated(molid) ) deallocate(molid)
	    if ( allocated(pref) ) deallocate(pref)
	    if ( allocated(pavlref) ) deallocate(pavlref)
	    if ( allocated(tmptab) ) deallocate(tmptab)
	    if ( allocated(sunrad) ) deallocate(sunrad)
	    if ( allocated(wvptmp) ) deallocate(wvptmp)
	    if ( allocated(vwvn) ) deallocate(vwvn)
	    if ( allocated(coef) ) deallocate(coef)


	    if ( allocated(kfix_ir) ) deallocate(kfix_ir)
	    if ( allocated(kh2o_ir) ) deallocate(kh2o_ir)



	    if ( allocated(dkh2o_ir) ) deallocate(dkh2o_ir)
	    if ( allocated(kvar) ) deallocate(kvar)
	    if ( allocated(nmols) ) deallocate(nmols)

	    if ( allocated(imols) ) deallocate(imols)
	    if ( allocated(nch) ) deallocate(nch)

	    if ( allocated(ichmap) ) deallocate(ichmap)
	    if (allocated(tmpself)) deallocate(tmpself)
	    if (allocated(kself)) deallocate(kself)
	    if (allocated(userIndex2chSetIndexMap)) deallocate(userIndex2chSetIndexMap)


	    if ( allocated(VarMolIndex) ) deallocate(VarMolIndex)
	    if ( allocated(VarMolID) ) deallocate(VarMolID)

	    
	    return

	  end subroutine release

          ! Helper function
          ! FindFreeUnit
          ! return a free UNIT
          integer function findFreeUnit()
              logical isOpened
              integer iostat
              do findFreeUnit = 10, 200
                     inquire(UNIT = FindFreeUnit, OPENED = isOpened, iostat = iostat)
                     if (iostat .ne. 0 ) cycle
                     if (.not. isOpened) return
              end do
              STOP '[FindFreeUnit] No Free File Unit Found'
          end function findFreeUnit

end module oss_ir
