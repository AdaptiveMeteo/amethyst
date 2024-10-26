program odcmpr2nc
  use netcdf
  implicit none
   !------------------------------------------------------------------------
   !     public items made available to calling program(s)   
   !------------------------------------------------------------------------

  !- Input vars
  integer :: narg
  integer , parameter :: stderr = 0
  integer , parameter :: stdout = 6
  integer, parameter :: LUT_KIND = KIND(1)
  character(len=256) :: basef
  
  !-------------------------
  !- Netcdf variariables   |
  !-------------------------
  integer :: ierr , iomode
  integer :: ncid , ncstat  , icwvnvar, isunradvar , imolidfixvar , imolidvar
  integer :: iprefvar , itmptabvar ,itmpselfvar,ikselfvar, iwvptabvar , icoefvar, ifmtvar,ixdvar
  integer :: iichmapvar , inchvar ,inchlistvar, iselsvar , ivwvnvar , imolsvar,ichindxvar,imolprofvar
  integer :: ikfixvar , idkh2ovar , ikh2ovar , ikvarvar, innodesvar, ichannelsvar,ifreqvar,inmolsvar, inmols_var
  integer :: ifixdmrvar, idfltvar
  !- dimension ID
  integer :: imolidfixdimid , imoliddimid , ilevdimid , intmpoddimid,imxmolsdimid,ikvarnsize1dimid
  integer :: intselfdimid,ikselfdimid,itmpselfdimid,imxindxdimid,ichdimid,ikh2onsize2dimid, inwvpoddimid
  integer :: inlayoddimid , inmoltabdimid , infdimid , inchmaxdimid,infsmpdimid,ikfixsize1dimid,iscalardimid
	

  integer :: iuo, ius,iuf
  integer :: uid_sel , uid_od
  integer :: nchan , nf_sel , nchmax , nfsmp , nmolx , nn
  integer :: nlayod , ntmpod , nlev , ismp , i , j , n
  integer(2) :: nmolfix
  integer(2) , allocatable , dimension(:) :: molidfix
  real , dimension(:) , allocatable :: cwvn
  real , dimension(:) , allocatable :: pref
  real(8) , dimension(:) , allocatable :: vwvn
  real , dimension(:,:) , allocatable :: tmptab , wvptab , coef
  real , dimension(:,:,:) , allocatable :: kfix , dkh2o , kh2o
  real , dimension(:,:,:,:) , allocatable :: kvar
  integer , dimension(:,:) , allocatable :: ichmap
  integer , dimension(:) , allocatable :: isels
  integer(2) , dimension(:) , allocatable :: nch
  integer                   :: nsize1,nsize2,nsize3,nsize4,nsize5

  CHARACTER(len=100) :: instr_info
  integer , dimension(4) :: dimids

  INTEGER, PARAMETER            :: mxhmol=27
  INTEGER                       :: nprof,k
  REAL, ALLOCATABLE, target     :: xG(:)
  real, allocatable             :: pUser(:),zProf(:)
  REAL, DIMENSION(:), pointer   :: t => null()
  REAL, DIMENSION(:), pointer   :: w => null() 
  REAL, pointer                 :: tskin => null()

  REAL                          :: psfc,obsAngle,solZenith,azAngle
  integer, dimension(mxhmol)    :: varMolID
  INTEGER                       :: tempIndex,tSkinIndex,pSurfIndex,varMolIndex(mxhmol)
  INTEGER                       :: nVarMol 
  REAL, ALLOCATABLE             :: chanFreq(:)
  INTEGER, ALLOCATABLE          :: chanIndex(:)
  REAL                          :: lat
  REAL , ALLOCATABLE            :: surfEmRfGrid(:),surfEmRefl(:,:)
  REAL , ALLOCATABLE            :: xkEmRf(:,:,:)
  REAL, ALLOCATABLE             :: yf(:)
  REAL, ALLOCATABLE             :: xk(:,:)
  REAL, ALLOCATABLE             :: molProf(:,:)

  logical                       :: flatSurfOP = .TRUE.
  integer                       :: lambertian = 0
  INTEGER                :: nRef
  REAL            :: em,  refl, zSurf
  INTEGER(KIND=4) :: iAlt, Iprof, iRefFixVar
  INTEGER(KIND=4) :: obsLevel, nparG, nsf, jj
  CHARACTER(len=120)  :: fn
  integer(kind=4)   :: numberNodes, currentSelection
  CHARACTER  :: ch
  integer :: iset,nSets,nChList
  integer, dimension(:), allocatable :: chanListTmp
  integer              :: indxChSel
  integer              :: fmtLUT_Global, xid_Global


   !------------------------------------------------------------------------
   !     Parameters/Dimensions
   !------------------------------------------------------------------------
   integer, parameter     :: UNDEFINED_INDEX = -999
   integer, parameter     :: mxIndex=10

   integer, parameter :: mxlev=205,mxlay=mxlev-1
   integer, parameter :: mxSf=50,MxmolS=20
   integer, parameter :: mxParG=(mxHmol+1)*mxlev
   real,    parameter :: pi=acos(-1.0), deg2rad = pi/180.
   real,    parameter :: dbAvg2dtau=-1./12!Weighting factor derivative wrt optical depth, for linear-in-tau, in low-OD limit
   real,    parameter :: grav_const_req1= 9.80665 ! m / s^2
   real,    parameter :: grav_const_req2=-0.02586 ! m / s^2
   real,    parameter :: Rgas=8.3144621   ! Universal gas constant (J / mol / K)
   real,    parameter :: drymwt=28.964 ! Dry air molecular mass (g / mol)
   real,    parameter :: AVOGAD = 6.02214199E+23 ! Avogadro ( 1 / mol )
   integer, parameter :: MAXSPC=84
   real, dimension(MAXSPC), parameter  :: molWt=(/ &
                18.015, 44.010, 47.998,44.010, 28.011, 16.043, 31.999, &   ! 7
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
                 
   real,    parameter :: mperkm=1000. ! (m/km)
   real,    parameter :: rairDflt=10.0/grav_const_req1 !10./g ->converts p(mb) in g/cm**2
   real, parameter :: rConst =  Rgas/(drymwt*1.e-3)/grav_const_req1/mperkm
   real,    parameter :: rEarth = 6371.23
   real,    parameter :: pMatchLim = 1.e-3
   integer, parameter :: magicNumber = 123456789

   !------------------------------------------------------------------------
   !     Pressure Grid / Molecules / Geophysical Pointers
   !------------------------------------------------------------------------
   integer   :: nMol
   integer, allocatable                 :: molID(:)

   !------------------------------------------------------------------------
   !     Instrument parameters
   !------------------------------------------------------------------------
   integer,save                         :: iChSet=UNDEFINED_INDEX  ! Current channel set
   integer,save                         :: lastChSet    ! Number of loaded chan sets
   integer,save                         :: nchanAll
   integer,save                         :: nfsmp_ir,nNodes
   integer, allocatable                 :: nChList_arr(:)
   integer, allocatable                 :: chanList_arr(:,:)
   integer, dimension(0:mxIndex)        :: userIndex2chSetIndexMap

   
   integer, dimension(:)  , allocatable :: nNodes_arr
   integer,                 allocatable :: iselS_arr(:,:),isel_Loc(:)
   integer*2, dimension(:), allocatable :: nch_ir,nch_Loc
   integer,                 allocatable :: nch_arr(:,:)           
   real, dimension(:,:)   , allocatable :: coef_ir
   real,                    allocatable :: coef_arr(:,:,:)
   real(KIND=LUT_KIND),     allocatable :: coef_Loc(:,:)
   integer, dimension(:,:), allocatable :: ichmap_ir
   integer,                 allocatable :: ichMap_arr(:,:,:)
   integer,                 allocatable :: ichmap_Loc(:,:)
   real, dimension(:)     , allocatable :: pavlref,sunrad
   real, dimension(:)     , allocatable :: f1_arr,f2_arr
   !------------------------------------------------------------------------
   !     Lookup tables 
   !------------------------------------------------------------------------
   real(KIND = LUT_KIND), allocatable :: TmpSelf(:)

   integer                                  :: NWvpOD,NTself
   real(KIND=LUT_KIND),      dimension (:,:) , allocatable :: kfix_ir,kh2o_ir,kvar_ir
   real,      dimension (:,:) , allocatable :: dkh2o_ir
   real,      dimension (:,:) , allocatable :: kself
   integer, dimension (:,:) , allocatable   :: ImolS_output
   integer,dimension (:)   , allocatable    :: nmols_output	


   ! calculate WV Jacobian correctly
   real, dimension(MxLay,MxHmol)     :: mydwqu, mydwql
   !-------------------------------------------------------------ù
   integer, allocatable      :: mvar(:), chIndx(:), iPol(:)
   integer                               :: iuSol
   CHARACTER*120             :: selFile,lutFile
   integer, allocatable      :: iSelFix(:),mFix(:),invMolID(:)
   integer                   :: nselFix,nFix,nLevF,icnt

    real(KIND=LUT_KIND)       :: V1m,V2m
    real(KIND=LUT_KIND),allocatable :: cWvnLoc(:),prefLoc(:),fix(:)
    real(KIND=LUT_KIND),allocatable :: tmpLoc(:,:),wvpLoc(:,:),dfltLoc(:,:)
    real(KIND=LUT_KIND),allocatable :: kfix_Loc(:),kh2o_Loc(:),kvar_Loc(:)
    real(KIND=LUT_KIND),allocatable :: dkfix_Loc(:),dkh2o_Loc(:)
    real(KIND=LUT_KIND),allocatable :: kself_Loc(:)
    
    real, dimension (:,:) , allocatable :: dkFix_ir
    real, dimension (:)   , allocatable :: fixDMR(:)
    real                                :: scale, dummyReal
    character(len=256)                  :: outfile 
    character(len=256)                  :: defProfFile="OSSdefaultMR101.27.dat"


   !- Take base filename as argument 
   call purpose
   narg = command_argument_count( )
   if ( narg < 1 ) then
     write(stderr,*) 'Not enough arguments. Need basis file name'
     call usage( )
   end if
   call get_command_argument(1,basef)
   write(lutFile,'(a,a)') trim(basef), '.lut'
   write(selFile,'(a,a)') trim(basef), '.sel'
   write(outfile,'(a,a)') trim(basef), '.nc'
   print *,'selFile:             '      ,  trim(selFile)
   print *,'lutFile:             '      ,  trim(lutFile)

   ! Initialize OSS model
   print *, 'OSS model initialization '
   allocate(chanFreq(9000), chanIndex(9000))

   !---read OD and SEL data
   call GetOD(selFile,lutFile,defProfFile)

   !---Size conformity checks (safer interface with calling program/subroutine)
   if (iChSet .eq. UNDEFINED_INDEX) then
        if (SIZE(chanFreq) < nChList_arr(0)) THEN
           print*, 'Err[oss_ir_module::oss_init_ir]: Insufficient Size for chanFreq'
           call exit(1)
        end if
        if (SIZE(chanIndex) < nChList_arr(0)) THEN
           print*, 'Err[oss_ir_module::oss_init_ir]: Insufficient Size for chanIndex'
           call exit(1)
        end if
    else    
        if (SIZE(chanFreq) < nChList_arr(iChSet)) THEN
           print*, 'Err[oss_ir_module::oss_init_ir]: Insufficient Size for chanFreq'
           call exit(1)
        end if
        if (SIZE(chanIndex) < nChList_arr(iChSet)) THEN
           print*, 'Err[oss_ir_module::oss_init_ir]: Insufficient Size for chanIndex'
           call exit(1)
        end if
    end if
    if (SIZE(varMolID)  > mxHmol) THEN
       print*, 'Err[oss_ir_module::oss_init_ir]: MolID Vector too large'
       call exit(1)
    end if

    if (SIZE(pRef) < nlev) then 
       print*, 'Err[oss_ir_module::oss_init_ir] Insufficient Size for pRef'
       call exit(1)
    end if    
 
    allocate (sunrad(nfsmp_ir))
    call InterpSolar(NFSmp_ir,vwvn,sunrad) 


    !------------------------------------------------------------------------
    ! Copy to output arguments
    !------------------------------------------------------------------------
    ! 0th channel assumed
    nRef = nlev
    pRef(1:nlev)                       = Pref(1:nlev)
    nChan                              = nChList_arr(0)
    chanFreq(1:nChan)                  = cWvn(chanList_arr(1:nChan,0))
    chanIndex(1:nChan)                 = chanList_arr(1:nChan,0)
    userIndex2chSetIndexMap(1:mxIndex) = -1
    userIndex2chSetIndexMap(0)         = 0
   !------------------------------------------------------------------         
  
   call write_netcdf(outfile)

  call ossDestroy_ir()
  write(6,*) 'Done!'
  return 
 CONTAINS

	SUBROUTINE InterpSolar(nptout,wvnout,sunout)

	!<f90Subroutine>********************************************************
	!
	! NAME:
	!
	!   InterpSolar
	!
	! PURPOSE:
	!
	!   Interpolate solar source data from input grid to output grid.
	!
	! SYNTAX:
	!
	!   CALL InterpSolar(nptout, wvnout, sunout)
	!
	! ARGUMENTS:
	!
	!   INPUTS:
	!   
	!   nptout  INTEGER  Number of output grid
	!   wvnout  REAL*8   Wavenumbers on output grid
	!   
	!   INPUTS/OUTPUTS:
	!   
	!   sunout  REAL     Interpolated solar contribution on output grid
	!
	!   * OPTIONAL
	!
	! INCLUDES:
	!
	!   None
	!
	!*******************************************************</f90Subroutine>

	!c***********************************************************************
	!c* Function name: InterpSolar
	!c* Purpose: Interpolation routine
	!c* Usage: call InterpSolar(nptout,wvnout,sunout)
	!c*
	!c* Description: Interpolate sunin from input grid to output grid
	!c*              Set sunout to zero if nptsSun equal zero
	!c*
	!c* Inputs:
	!c* Var_name     Type       Description
	!c* --------     ----       -----------
	!c* sunin        real       input data (already in module)
	!c* nptin        integer    dimension of input (already in module)
	!c* nptout       integer    dimension of input
	!c*
	!c* Outputs:
	!c* Var_name     Type       Description
	!c* --------     ----       -----------
	!c* sunout       real       input data
	!c*
	!c* Common blocks: defined the includes
	!c* Includes:
	!c* Name        Description
	!c* ----        -----------
	!c* Externals: none
	!c*
	!c*<copyright_OSS>=====================================================
	!c*
	!c* Developed by Atmospheric and Environmental Research, Inc.                
	!c*
	!c* Copyright: Atmospheric and Environmental Inc., 1997-2007
	!c*                 All Rights Reserved
	!c*
	!c====================================================</copyright_OSS>
	!c***********************************************************************
	  INTEGER,                 INTENT(IN)     :: nptout
	  REAL*8,  DIMENSION(:),   INTENT(IN)     :: wvnout
	  REAL,    DIMENSION(:),   INTENT(INOUT)  :: sunout

          INTEGER                                 :: Lsol
          character(len=256)                      :: solFile="newkur.dat"

          !---Local variables
          CHARACTER*10         :: header  ! change format 10 if this changes
          INTEGER              :: ipt

          INTEGER              :: nptsSun
	  INTEGER              :: i,k
          REAL,    DIMENSION(:), ALLOCATABLE  :: sunwvn
          REAL,    DIMENSION(:), ALLOCATABLE  :: sunrad
          
          !---read solar file
          iuSol = findFreeUnit()
          nptsSun=0

          ! open solar file
          !filename=TRIM(solFile)
          open(unit=Lsol,file=solFile,status='old')

          ! read file if it exists
          read(Lsol,*)  nptsSun
          read(Lsol,10) header
          read(Lsol,10) header
          10    format(a10)

          allocate(sunrad(nptsSun),sunwvn(nptsSun))

          do 100 ipt=1,nptsSun
             read(Lsol,*) sunwvn(ipt),sunrad(ipt)
          100  continue
          sunrad=sunrad*1.0e+7  ! convert from (W cm-2 / cm-1) to (mW m-2 / cm-1)

          close(Lsol)

          ! Start InterpSolar
	  if (nptsSun.eq.0) then
	     sunout=0.
	     return
	  endif
	
	  do i=1,nptout
	     if (wvnout(i).le.sunwvn(1))then
	        sunout(i)=sunrad(1)+(sunrad(2)-sunrad(1))*(wvnout(i)-sunwvn(1))/ &
	             (sunwvn(2)-sunwvn(1))
	     else if (wvnout(i).ge.sunwvn(nptsSun)) then
	        sunout(i)=sunrad(nptsSun-1)+(sunrad(nptsSun)-sunrad(nptsSun-1))* &
	             (wvnout(i)-sunwvn(nptsSun-1)) &
	             /(sunwvn(nptsSun)-sunwvn(nptsSun-1))
	     else
	        do k=1,nptsSun-1
	           if (wvnout(i).ge.sunwvn(k).and. &
	                wvnout(i).lt.sunwvn(k+1))then
	              sunout(i)=sunrad(k)+ &
	                   (sunrad(k+1)-sunrad(k))*(wvnout(i)-sunwvn(k))/ &
	                   (sunwvn(k+1)-sunwvn(k))
	           endif
	        end do
	     end if
	  end do

	  return

          !jump here if error on open
          200  continue
          write(*,*) 'solar function not used'
	  return

	END SUBROUTINE InterpSolar

        subroutine read_mol_prof(iuf,ref_molid,vmol)
          !Input variables
          integer                 ,intent(in) :: iuF
          integer                 ,intent(in) :: ref_molid
          !Output variables
          real,    dimension(:),intent(inout) :: vmol
          !Local variables
          integer                :: k

          read(iuF,*)k 
          if (ref_molid == 1 .AND. k /= 1) THEN
             print*, 'Err[oss_ir_module::GetOD]: Include H2O mixing ratio into defProfFile'
             call exit(1)
          end if
          read(iuF,*)vmol(1:nlayOD+1)
          if (k == ref_molid) return
          do while (k /= ref_molid)
             read(iuF,*)k 
             read(iuF,*)vmol(1:nlayOD+1)
          end do
          return
        end subroutine

	  subroutine write_netcdf(outfile)
                  !- Write output NETCDF

		  integer , parameter		 :: stderr = 0
                  integer , parameter 		 :: stdout = 6
                  character(len=256), intent(in) :: outfile


		  write(stdout,*) 'Writing output file    : ',trim(outfile)
		  iomode = ior(ior(nf90_classic_model,nf90_clobber),nf90_netcdf4)
		  ncstat = nf90_create(outfile,iomode,ncid)

		  call check_ncerror('Cannot create output file '//trim(outfile))

		  ncstat = nf90_put_att(ncid,nf90_global,'id_sel',uid_sel)
		  call check_ncerror('Cannot add global attribute id_sel')
		  ncstat = nf90_put_att(ncid,nf90_global,'id_od',uid_od)
		  call check_ncerror('Cannot add global attribute id_od')
		  ncstat = nf90_put_att(ncid,nf90_global,'info',trim(instr_info))
		  call check_ncerror('Cannot add global attribute info')
		  ncstat = nf90_put_att(ncid,nf90_global,'v1',v1m)
		  call check_ncerror('Cannot add global attribute v1')
		  ncstat = nf90_put_att(ncid,nf90_global,'v2',v2m)
		  call check_ncerror('Cannot add global attribute v2')


		  !----------------------------------------------------
		  !- CREATE Dimensions
		  !----------------------------------------------------

		  !- Create scalar dimension - ID: iscalardimid
		  ncstat = nf90_def_dim(ncid,'scalar',1,iscalardimid)
		  call check_ncerror('Cannot create scalar dimension')

		  !- Create nchan dimension - ID: ichdimid
		  ncstat = nf90_def_dim(ncid,'nchan',nchan,ichdimid)
		  call check_ncerror('Cannot create nchan dimension')
		  
		  !- Create ntmpod dimension - ID: intmpoddimid 
		  ncstat = nf90_def_dim(ncid,'ntmpod',ntmpod,intmpoddimid)
		  call check_ncerror('Cannot create ntmpod dimension')

		  !- Create nlayod dimension - ID: inlayoddimid
		  ncstat = nf90_def_dim(ncid,'nlayod',nlayod,inlayoddimid)
		  call check_ncerror('Cannot create nlayod dimension')

		  !- Create ntself dimension - ID: intselfdimid
		  ncstat = nf90_def_dim(ncid,'ntself',ntself,intselfdimid)
		  call check_ncerror('Cannot create ntself dimension')

		  !- Create nfsmp_ir dimension - ID: infsmpdimid
		  ncstat = nf90_def_dim(ncid,'nfsmp',nfsmp_ir,infsmpdimid)
		  call check_ncerror('Cannot create nfsmpr dimension')
		  
		  !- Create nchmax dimension - ID: inchmaxdimid
		  ncstat = nf90_def_dim(ncid,'nchmax',nchmax,inchmaxdimid)
		  call check_ncerror('Cannot create nchmax dimension')

		  !- Create nf_sel dimension - ID: infdimid
		  ncstat = nf90_def_dim(ncid,'nf',nf_sel,infdimid)
		  call check_ncerror('Cannot create nf dimension')

		  !- Create mxindex dimension - ID: imxindx
		  ncstat = nf90_def_dim(ncid,'mxindex',mxindex + 1,imxindxdimid)
		  call check_ncerror('Cannot create mxindx dimension')

		  !- Create kfix first dimension - ID: ikfixsize1dimid
		  ncstat = nf90_def_dim(ncid,'kfix_nsize1',nlayod*ntmpod,ikfixsize1dimid)
		  call check_ncerror('Cannot create kfix_nsize1 dimension')

		  !- Create kh2o first dimension - ID: ikh2odimid
		  ncstat = nf90_def_dim(ncid,'kh2o_nsize2',nlayod*ntmpod*nwvpod,ikh2onsize2dimid)
		  call check_ncerror('Cannot create kh2o_nsize2 dimension')

		  !- Create kvar first dimension - ID: ikvarnsize1dimid
		  ncstat = nf90_def_dim(ncid,'kvar_nsize1',mxmols*nlayod*ntmpod,ikvarnsize1dimid)
		  call check_ncerror('Cannot create kvar_nsize1 dimension')
		  
		  !- Create nlev first dimension - ID: ilevdimid
		  ncstat = nf90_def_dim(ncid,'nlev',nlayod + 1,ilevdimid)
		  call check_ncerror('Cannot create nlev dimension')

		  !- Create molid dimension - ID: imoliddimid
		  ncstat = nf90_def_dim(ncid,'nmol',nmol,imoliddimid)
		  call check_ncerror('Cannot create nmol dimension')
		
		  ncstat = nf90_def_dim(ncid,'nmoltab',nmol+2,inmoltabdimid)
		  call check_ncerror('Cannot create nmoltab dimension')

		  ncstat = nf90_def_dim(ncid,'mxmols',mxhmol,imxmolsdimid)
		  call check_ncerror('Cannot create mxmols dimension')

		  ncstat = nf90_def_dim(ncid,'nwvpod',max(nwvpod-1,1),inwvpoddimid)
		  call check_ncerror('Cannot create nwvpod dimension')


		  !----------------------------------------------------
		  !- CREATE Variables
		  !----------------------------------------------------
		  
		  dimids(1) = imoliddimid
		  dimids(2) = inlayoddimid
		  ncstat = nf90_def_var(ncid,'dflt',nf90_float,dimids(1:2),idfltvar)
		  call check_ncerror('Cannot create var dflt')

		  dimids(1) = inwvpoddimid
		  dimids(2) = inlayoddimid
		  ncstat = nf90_def_var(ncid,'wvptab',nf90_float,dimids(1:2),iwvptabvar)
		  call check_ncerror('Cannot create var wvptab')
		  
		  dimids(1) = inlayoddimid
		  ncstat = nf90_def_var(ncid,'fixDMR',nf90_float,dimids(1:1),ifixdmrvar)
		  call check_ncerror('Cannot create var fixDMR')
		  !- Create scalar var fmtLUT and xid
		  dimids(1) = iscalardimid
		  ncstat = nf90_def_var(ncid,'fmtLUT',nf90_int,dimids(1:1),ifmtvar)
		  ncstat = nf90_def_var(ncid,'xid',nf90_int,dimids(1:1),ixdvar)
		  call check_ncerror('Cannot create var fmtLUT')

		  !- Create array cWvn - ID: icwvnvar - DIM: nchan
		  dimids(1) = ichdimid
		  ncstat = nf90_def_var(ncid,'cWvn',nf90_float,dimids(1:1),icwvnvar)
		  call check_ncerror('Cannot create var cWvn')

		  dimids(1) = infsmpdimid
		  ncstat = nf90_def_var(ncid,'sunrad',nf90_float,dimids(1:1),isunradvar)
		  call check_ncerror('Cannot create var sunrad')
                 
		  !- Create array chindx - ID: ichindxvar - DIM: nchan
		  dimids(1) = ichdimid
		  ncstat = nf90_def_var(ncid,'chindx',nf90_float,dimids(1:1),ichindxvar)
		  call check_ncerror('Cannot create var chindx')

		  !- Create ndarray tmptab - ID: itmptabvar - DIM: ntmpod x nlayod  
		  dimids(1) = intmpoddimid
		  dimids(2) = inlayoddimid
		  ncstat = nf90_def_var(ncid,'tmptab',nf90_float,dimids(1:2),itmptabvar)
		  call check_ncerror('Cannot create var tmptab')

		  !- Create array tmpself - ID: itmpselfvar - DIM: ntself
		  dimids(1) = intselfdimid
		  ncstat = nf90_def_var(ncid,'tmpself',nf90_float,dimids(1:1),itmpselfvar)
		  call check_ncerror('Cannot create var tmpself')

		  !- Create array kself - ID: ikselfvar - DIM: ntself x nfsmp_ir
		  dimids(1) = intselfdimid
		  dimids(2) = infsmpdimid
		  ncstat = nf90_def_var(ncid,'kself',nf90_float,dimids(1:2),ikselfvar)
		  call check_ncerror('Cannot create var kself')

		  !- Create ndarray coef - ID: icoefvar - DIM: nchmax x nf_sel x (mxIndex + 1)
		  dimids(1) = inchmaxdimid
		  dimids(2) = infdimid
		  dimids(3) = imxindxdimid
		  ncstat = nf90_def_var(ncid,'coef',nf90_float,dimids(1:3),icoefvar)
		  call check_ncerror('Cannot create var coef')
		  ncstat = nf90_put_att(ncid,icoefvar,'_FillValue',-9999.0)
		  call check_ncerror('Cannot add attribute _FillValue for coef variable')

		  !- Create ndarray coef - ID: iichmapvar - DIM: nchmax x nf_sel x (mxIndex + 1)
		  dimids(1) = inchmaxdimid
		  dimids(2) = infdimid
		  dimids(3) = imxindxdimid
		  ncstat = nf90_def_var(ncid,'ichmap',nf90_int,dimids(1:3),iichmapvar)
		  call check_ncerror('Cannot create var ichmap')
		  ncstat = nf90_put_att(ncid,iichmapvar,'_FillValue',-1)
		  call check_ncerror('Cannot add attribute _FillValue for ichmap variable')

		  !- Create ndarray nch - ID: inchvar - DIM: nf_sel x mxIndex
		  ncstat = nf90_def_var(ncid,'nch',nf90_short,dimids(2:3),inchvar)
		  call check_ncerror('Cannot create var nch')

		  !- Create array vwvn - ID: ivwvnvar - DIM: nfsmp_ir
		  dimids(1) =  infsmpdimid
		  ncstat = nf90_def_var(ncid,'vwvn',nf90_double,dimids(1:1),ivwvnvar)
		  call check_ncerror('Cannot create var vwvn')
		  
		  !- Create ndarray isels_arr - ID: iselsvar - DIM: nf_sel x (mxIndex + 1)
		  dimids(1) = infdimid
		  dimids(2) = imxindxdimid
		  ncstat = nf90_def_var(ncid,'isels',nf90_int,dimids(1:2),iselsvar)
		  call check_ncerror('Cannot create var isels')

		  !- Create ndarray nmols - ID: inmols_var - DIM: nfsmp_ir
		  dimids(1) = infsmpdimid
		  ncstat = nf90_def_var(ncid,'nmols_',nf90_float,dimids(1:1),inmols_var)
		  call check_ncerror('Cannot create var nmols_')

		  !- Create ndarray imols_ - ID: imolsvar - DIM: nfsmp_ir X mxhmols
		  dimids(1) = infsmpdimid
		  dimids(2) = imxmolsdimid
		  ncstat = nf90_def_var(ncid,'imols_',nf90_short,dimids(1:2),imolsvar)
		  call check_ncerror('Cannot create var imols_')
		  ncstat = nf90_put_att(ncid,imolsvar,'_FillValue',-1_2)
		  call check_ncerror('Cannot add attribute _FillValue for imols variable')

		  !- Create ndarray kfix - ID: ikfixvar - DIM: kfix_nsize1 x nfsmp_ir
		  dimids(1) = ikfixsize1dimid
		  dimids(2) = infsmpdimid
		  ncstat = nf90_def_var(ncid,'kfix',nf90_real,dimids(1:2),ikfixvar)
		  call check_ncerror('Cannot create var kfix')
		  ncstat = nf90_put_att(ncid,ikfixvar,'_FillValue',-9999.0)
		  call check_ncerror('Cannot add attribute _FillValue for kfix variable')
		  ncstat = nf90_def_var_deflate(ncid,ikfixvar,1,1,9)
		  call check_ncerror('Cannot set deflate level to var kfix')



		  !- Create ndarray dkh2o - ID: idkh2ovar - DIM: kfix_nsize1 x nfsmp_ir
		  ncstat = nf90_def_var(ncid,'dkh2o',nf90_real,dimids(1:2),idkh2ovar)
		  call check_ncerror('Cannot create var dkh2o')
		  ncstat = nf90_put_att(ncid,idkh2ovar,'_FillValue',-9999.0)
		  call check_ncerror('Cannot add attribute _FillValue for dkh2o variable')
		  ncstat = nf90_def_var_deflate(ncid,idkh2ovar,1,1,9)
		  call check_ncerror('Cannot set deflate level to var dkh2o')

		  !- Create ndarray kh2o - ID: ikh2ovar - DIM: dh2o_nsize2 x nfsmp_ir
                  dimids(1) = ikh2onsize2dimid
		  ncstat = nf90_def_var(ncid,'kh2o',nf90_real,dimids(1:2),ikh2ovar)
		  call check_ncerror('Cannot create var kh2o')
		  ncstat = nf90_put_att(ncid,ikh2ovar,'_FillValue',-9999.0)
		  call check_ncerror('Cannot add attribute _FillValue for kh2o variable')
		  ncstat = nf90_def_var_deflate(ncid,ikh2ovar,1,1,9)
		  call check_ncerror('Cannot set deflate level to var kh2o')

		  !- Create ndarray kvar - ID: ikvar - DIM: kvar_nsize1 x nfsmp_ir
		  dimids(1) = ikvarnsize1dimid
		  dimids(2) = infsmpdimid
		  ncstat = nf90_def_var(ncid,'kvar',nf90_real,dimids(1:2),ikvarvar)
		  call check_ncerror('Cannot create var kvar')
		  ncstat = nf90_put_att(ncid,ikvarvar,'_FillValue',-9999.0)
		  call check_ncerror('Cannot add attribute _FillValue for kvar variable')
		  ncstat = nf90_def_var_deflate(ncid,ikvarvar,1,1,9)
		  call check_ncerror('Cannot set deflate level to var kvar')


		  !- Create array nNodes - ID: innodesvar - DIM: (mxindex + 1 )
		  dimids(1) = imxindxdimid
		  ncstat = nf90_def_var(ncid,'nnodes',nf90_real,dimids(1:1),innodesvar)
		  call check_ncerror('Cannot create var nnodes')

		  !- Create array channels - ID: ichannelsvar - DIM: nchan x (mxindex + 1 )
		  dimids(1) = ichdimid
		  dimids(2) = imxindxdimid
		  ncstat = nf90_def_var(ncid,'channels',nf90_real,dimids(1:2),ichannelsvar)
		  call check_ncerror('Cannot create var channels')

		  !- Create array nchlist - ID: inchlistvar - DIM:(mxindex + 1 )
		  dimids(1) = imxindxdimid
		  ncstat = nf90_def_var(ncid,'nchlist',nf90_real,dimids(1:1),inchlistvar)
		  call check_ncerror('Cannot create var nchlist')

		  !- Create array molid - ID: imolidvar - DIM: imoliddimid  
		  dimids(1) = imoliddimid
		  ncstat = nf90_def_var(ncid,'molid',nf90_short,dimids(1:1),imolidvar)
		  call check_ncerror('Cannot create var molid')
		  
		  !- Create array pref - ID: iprefvar - DIM: nlev
		  dimids(1) = ilevdimid
		  ncstat = nf90_def_var(ncid,'pref',nf90_double,dimids(1:1),iprefvar)
		  call check_ncerror('Cannot create var pref')

		  !- Create array chanfreq - ID: ifreqvar - DIM: nchan
		  dimids(1) = ichdimid
		  ncstat = nf90_def_var(ncid,'chanfreq',nf90_float,dimids(1:1),ifreqvar)
		  call check_ncerror('Cannot create var chanfreq')

	
		  !- Create array molProf - ID: imolprofvar - DIM: nmol x nlayod
		  dimids(1) = imoliddimid
		  dimids(2) = ilevdimid
		  ncstat = nf90_def_var(ncid,'molProf',nf90_float,dimids(1:2),imolprofvar)
		  call check_ncerror('Cannot create var molProf')

		  ncstat = nf90_enddef(ncid)
		  call check_ncerror('Cannot finalize output file '//trim(outfile))

		  !----------------------------------------------------
		  !- WRITE Variables
		  !----------------------------------------------------

		  !- Write array fixDMR - DIM: ifixdmrvar
		  ncstat = nf90_put_var(ncid,ifixdmrvar,fix)
		  call check_ncerror('Cannot write variable fixDMR')

		  !- Write scalr vars fmtLUT and xid
		  ncstat = nf90_put_var(ncid,ifmtvar,fmtLUT_Global)
		  ncstat = nf90_put_var(ncid,ixdvar,xid_Global)
		  call check_ncerror('Cannot write variable fmtLUT')

		  !- Write array cWvn - DIM: nchan
		  ncstat = nf90_put_var(ncid,icwvnvar,cwvn)
		  call check_ncerror('Cannot write variable cWvn')

		  !- Write array sunrad - DIM: nfsmp
		  ncstat = nf90_put_var(ncid,isunradvar,sunrad)
		  call check_ncerror('Cannot write variable sunrad')

		  !- Write array chindex - ID: ichindxvar - DIM: nchan         
		  ncstat = nf90_put_var(ncid,ichindxvar,chindx)
		  call check_ncerror('Cannot write variable chindx')
		  
		  !- Write ndarray tmptab - ID: itmptabvar - DIM: ntmpod x nlayod
		  ncstat = nf90_put_var(ncid,itmptabvar,tmptab)
		  call check_ncerror('Cannot write variable tmptab')

		  !- Write array tmpself - ID: itmpselfvar - DIM: ntself
		  ncstat = nf90_put_var(ncid,itmpselfvar,tmpself)
		  call check_ncerror('Cannot write variable tmpself')
		  
		  !- Write ndarray kself - ID: ikselfdimid - DIM: ntself x nfsmp_ir
		  ncstat = nf90_put_var(ncid,ikselfvar,kself)
		  call check_ncerror('Cannot write variable kself')

		  !- Write array coef_arr - ID: icoefvar - DIM: nchmax x nf_sel x (mxIndex + 1)
		  ncstat = nf90_put_var(ncid,icoefvar,coef_arr)
		  call check_ncerror('Cannot write variable coef')
  	

		  !- Write array ichmap_arr - ID: iichmapvar - DIM: nchmax x nf_sel x (mxIndex + 1)
		  ncstat = nf90_put_var(ncid,iichmapvar,ichmap_arr)
		  call check_ncerror('Cannot write variable ichmap')

		  !- Write array nch_arr - ID: inchvar - DIM: nf_sel x ( mxIndex + 1)
		  ncstat = nf90_put_var(ncid,inchvar,nch_arr)
		  call check_ncerror('Cannot write variable nch')

		  !- Write array vwvn - ID: ivwvnvar - DIM: nfsmp_ir
		  ncstat = nf90_put_var(ncid,ivwvnvar,vwvn)
		  call check_ncerror('Cannot write variable vwvn')

		  !- Write ndarray isels - ID: iselsvar - DIM: nf_sel x ( mxIndex + 1)
		  ncstat = nf90_put_var(ncid,iselsvar,isels_arr)
		  call check_ncerror('Cannot write variable isels')

		  !- Write ndarray imols - ID: imolsvar - DIM: mxmols x nfsmp_ir
		  ncstat = nf90_put_var(ncid,imolsvar,imols_output)
		  call check_ncerror('Cannot write variable imols_')

		  !- Write ndarray kfix - ID: ikfixvar - DIM: kfix_nsize1 x nfsmp_ir
		  ncstat = nf90_put_var(ncid,ikfixvar,kfix_ir)
		  call check_ncerror('Cannot write variable kfix')
		  !- Write ndarray dkh2o - ID: idkh2ovar - DIM: kfix_nsize1 x nfsmp_ir
		  ncstat = nf90_put_var(ncid,idkh2ovar,dkh2o_ir)
		  call check_ncerror('Cannot write variable dkh2o')


		  !- Write ndarray kh2o - ID: ikh2ovar - DIM: dh2o_nsize2 x nfsmp_ir
		  ncstat = nf90_put_var(ncid,ikh2ovar,kh2o_ir)
		  call check_ncerror('Cannot write variable kh2o')

		  !- Write ndarray kvar - ID: ikvar - DIM: kvar_nsize1 x nfsmp_ir
		  ncstat = nf90_put_var(ncid,ikvarvar, kvar_ir)
		  call check_ncerror('Cannot write variable kvar')

		  !- Write array nNodes - ID: innodesvar - DIM: (mxindex + 1 )
		  ncstat = nf90_put_var(ncid,innodesvar, nNodes_arr)
		  call check_ncerror('Cannot write var nnodes')

		  !- Write array nNodes - ID: ichannelsvar - DIM: (mxindex + 1 )
		  ncstat = nf90_put_var(ncid,ichannelsvar, chanList_arr)
		  call check_ncerror('Cannot write var channels')

		  !- Write array nchlist - ID: inchlistvar - DIM: (mxindex + 1)
		  ncstat = nf90_put_var(ncid,inchlistvar, nchList_arr)
		  call check_ncerror('Cannot write var nchlist')

		  !- Write array pref - ID: iprefvar - DIM: nlev
		  ncstat = nf90_put_var(ncid,iprefvar,pref)
		  call check_ncerror('Cannot write variable pref')

		  !- Write array chanfreq - ID: ifreqvar - DIM: nchan
		  ncstat = nf90_put_var(ncid,ifreqvar,chanFreq(1:nChan))
		  call check_ncerror('Cannot write variable chanfreq')

		  !- Write array nmols_ - ID: inmolsvar - DIM: nfsmp_ir
		  ncstat = nf90_put_var(ncid,inmols_var,nmols_output)
		  call check_ncerror('Cannot write variable nmols_')

		  !- Write array molProf - ID: imolprofvar - DIM: nfsmp_ir
		  ncstat = nf90_put_var(ncid,imolprofvar,molProf)
		  call check_ncerror('Cannot write var molProf')

		  ncstat = nf90_put_var(ncid,imolidvar,molid)
		  call check_ncerror('Cannot write variable imolid')

		  ncstat = nf90_put_var(ncid,iwvptabvar,wvptab)
		  call check_ncerror('Cannot write variable wvptab')
		  ncstat = nf90_put_var(ncid,idfltvar,dfltLoc)
		  call check_ncerror('Cannot write variable dflt')
		  ncstat = nf90_close(ncid)
		  call check_ncerror('Cannot close output file '//trim(outfile))
		                    
	  end subroutine write_netcdf

	  subroutine GetOD(selfile,odfile,defProfFile)
	    CHARACTER(len=*),      intent(in)  :: selfile,odfile
            character(len=256),    intent(in)  :: defProfFile

	    real , parameter          :: SOfL = 29.9792458 !speed of light in km/sec, scaled by 1e-4
	    real , parameter          :: f1   = 1.19106d4 !stands in numerator in dependence of radiance on TB, units: (mW/(m^2 ster cm^-1)) / [cm^-1)^3, scaled by 1e9
	    real , parameter          :: f2   = 1.43879d3 !stands in exponential dependence of radiance on TB, like f2*wn*1e-2/TB
	    !---Input variables

	    real   , allocatable :: wvpTmp(:,:)

	    !---Local variables
	    integer                   :: ius,iuo
	    CHARACTER(len=100)        :: instr_info
	    CHARACTER(len=24)         :: date_IDSel,date_IDLUT
	    CHARACTER(len=20)         :: Sensor_ID,hdrOpt,Unit_Char_s,Unit_Char_l
	    CHARACTER*12, allocatable :: chID(:)
	    integer                   :: fmtver,nHdrOpt
	    integer                   :: fmtLUT,osstran_opt
	    integer                   :: WMO_Satellite_ID,WMO_Sensor_ID,Sensor_Type
	    integer                   :: chSetID
	    integer*2                 :: NmolS_tmp,ImolS_tmp(MxmolS),imol
	    integer                   :: imols_indx(MxmolS),imols_(MxmolS)
	    integer                   :: nmols_
	    integer                   :: uid_sel,uid_od,k,m,kk,ismp
	    integer                   :: i,nx,ks
	    integer                   :: ihO,Spc_Units_s,Spc_Units_l
	    real                      :: odfac = 0.
	    integer, allocatable      :: xid(:),xid_(:),xid_dpnd(:)
	    real   , allocatable      :: vmol(:),qr(:)
	    integer                   :: nwvpod_

	    !real(KIND=LUT_KIND)       :: V1m,V2m
	    real(KIND=LUT_KIND),allocatable :: cWvnLoc(:),prefLoc(:),fixLoc(:)
	    real(KIND=LUT_KIND),allocatable :: tmpLoc(:,:),wvpLoc(:,:)
	    real(KIND=LUT_KIND),allocatable :: kfix_Loc(:),kh2o_Loc(:),kvar_Loc(:)
	    real(KIND=LUT_KIND),allocatable :: dkfix_Loc(:),dkh2o_Loc(:)
	    real(KIND=LUT_KIND),allocatable :: kself_Loc(:)
	  	    
            real                            :: scale, dummyReal
	    logical                         :: isHDO
	    print *, 'getOD start'

	    !---open files
	    ius = findFreeUnit()
	    OPEN(ius,file=selfile,form='unformatted',status='old')

	    iuo = findFreeUnit()
	    OPEN(iuo,file=odfile,form='unformatted',status='old')

 	    iuf = findFreeUnit() 
	    OPEN (unit=iuf,file=defProfFile,status='old')

	    !- Read VarMolId from ref file
 	    read(iuf,*)nlev, nVarMol
	    read(iuf,*)varMolID(1:nVarMol)    

	    !---Check input molecular selection  
	    if (varMolID(1).NE.1) THEN
		print*, 'Err[oss_ir_module::GetOD]: Water vapor must be specified in molecular selection'
		call exit(1)
	    end if
	    do k=2,nVarMol
	       if (varMolID(k).LE.varMolID(k-1)) THEN
		  print*, 'Err[oss_ir_module::GetOD]: OSS IDs must be in ascending order'
		  call exit(1)
	       end if
	    end do
	    

	    !---read selFile header
	    read(ius)uid_sel  !magic number
	    if (uid_sel /= magicNumber) THEN
	       print*, 'Err[oss_ir_module::GetOD]: wrong Format of binary file: Sel'
	       call exit(1)
	    end if
	    read(ius)date_IDsel,fmtver              !date: YYYYMMDDHHMMSS.sss-0X000, X=4,or,5
	    read(ius)WMO_Satellite_ID,WMO_Sensor_ID !if specifically undefined: -1 -1
	    read(ius)Sensor_Type,Sensor_ID          !Sensor_Type: 1 for IR
	    read(ius)instr_info                     ! any useful information
	    read(ius)nChan,nf_sel,nchmax
	    allocate(chIndx(nChan),chID(nChan),iPol(nChan),cWvn(nChan),cWvnLoc(nChan))
	    read(ius)hdrOpt
	    hdrOpt  = adjustl(hdrOpt)
	    nHdrOpt = len_trim(hdrOpt)
	!defaults:
	    chIndx(1:nChan)=(/(i, i=1,nChan)/)
	    chID(1:nChan)='NA'
	    Spc_Units_s= 0      !default due to IR
	    Unit_Char_s= 'cm-1' !default due to IR
	    iPol(1:nChan)=-1    !decide about default (?)
	    do ihO = 1, nHdrOpt
	       SELECT CASE (HdrOpt(ihO:ihO))
		  CASE ('K')
		     read(ius) chIndx(1:nChan)
		  CASE ('I')
		     read(ius) chID(1:nChan)
		  CASE ('U')   ! for spectral units
		     read(ius) Spc_Units_s, Unit_char_s
		  CASE ('C')
		     read(ius)cWvnLoc(1:nChan)
		     if (Spc_Units_s == 1) cWvnLoc(1:nChan)= cWvnLoc(1:nChan)/SOfL!From GHz to cm^-1 
		     if (Spc_Units_s == 2) cWvnLoc(1:nChan)= 1.e4/cWvnLoc(1:nChan)!From mu  to cm^-1
		     cWvn(:)=cWvnLoc(:)
		  CASE ('P')
		     read(ius)iPol(1:nChan)
		  CASE DEFAULT
		     PRINT *,'Selection file header content flag not recognized:',HdrOpt(iHO:iHO)
	       end SELECT
	    end do
	    deallocate (cWvnLoc)
	    !---Get LUT file header
	    read(iuo)uid_od  !magic number
	    if (uid_od /= magicNumber) THEN
	       print*, 'Err[oss_ir_module::GetOD]: Wrong format of binary file: LUT'
	       call exit(1)
	    end if
	    read(iuo)date_IDlut,fmtLUT !date: YYYYMMDDHHMMSS.sss-0X000 (X: 4, or, 5)
	    fmtLUT_Global = fmtLUT
	    if (date_IDlut /= date_IDsel) THEN
	       print*,'Err[oss_ir_module::GetOD]: IDs for SEL-file and LUT-file are inconsistent '
	       call exit(1)
	    end if
	!    if (fmtLUT .lt. 2) THEN
	!       print*,'Err[oss_ir_module::GetOD]: the current OSS module requires fmtLUT of lutFile to be 2 '
	!       call exit(1)
	!    end if
	    read(iuo)HdrOpt
	    hdrOpt  = adjustl(hdrOpt)
	    nHdrOpt = len_trim(hdrOpt)
	    !Defaults
	    nx = 1 !For dK, number of dependent, e.g., dkFix, dkh2o, dkvar(:). Defualt
		   !requires dkh2o only
	    ALLOCATE (xid(nx), xid_dpnd(nx))
	    xid(1)      = 1
	    xid_dpnd(1) = 1
	    ossTran_Opt = 1
	    Spc_Units_l = 0
	    Unit_Char_l = 'cm-1'

	    do iHO = 1, nHdrOpt
	       SELECT CASE (HdrOpt(iHO:iHO))
		  CASE ('F')
		     read(iuo)ossTran_Opt ! 1 is for IR
		  CASE ('U')   ! for spectral units
		     read(iuo) Spc_Units_l, Unit_char_l
		  CASE ('X')
		     read(iuo)nx
		     if ( ALLOCATED(xid) ) DEALLOCATE(xid)
		     if ( ALLOCATED(xid_dpnd) ) DEALLOCATE(xid_dpnd)
		     if (nx > 0) THEN
		        ALLOCATE (xid(nx), xid_dpnd(nx))
		        read(iuo)xid(:), xid_dpnd(:)
		        if ( .not.ANY(xid(1:nx) == 1)) THEN
		           print*, 'Err[oss_ir_module::GetOD]: '// &
		                   'if nx > 0, WV should be among dependent on...'
		           call exit(1)
		        end if
		        if ( ANY(xid_dpnd(1:nx) /= 1)) THEN
		           print*, 'Err[oss_ir_module::GetOD]: '// &
		           'Currently, implementation is done for dependence on WV only'
		           call exit(1)
		        end if
		     else
		        print*, 'Err[oss_ir_module::GetOD]: '// &
		                'At least, one dependence should exist: WV on WV'
		        call exit(1)
		     end if
		  CASE ('N')
		     read(iuo)
		  CASE DEFAULT
		     PRINT *,'Selection file header content flag not recognized:',HdrOpt(iHO:iHO)
	       end SELECT
	    end do
	    xid_Global = xid(1)
	    read(iuo)V1m, V2m, nfsmp_ir
	    read(iuo)nfix,nMol !integer*4,

	    allocate(molid(nmol))
	    if (nMol .gt. MxHmol) then
	       print *, "Err[oss_ir_module::GetOD]: molID size is less than required by LUT file"
	       call exit(1)
	    end if
	    
	    ALLOCATE (mfix(nfix))
	    read(iuo) mfix(1:nfix),molID(1:nMol)
	    DEALLOCATE(mfix)

            if (any(molID == 81)) then
                 isHDO=.true.
            else
                 isHDO=.false.
            end if

	    read(iuo)nlayod,ntmpod,ntself,nwvpod
	    nwvpod_=MAX(nwvpod-1,1)
	!in this version we use nwvpod=1: 2-points WV ODs plus linear dependence of Keff
	!(Keff=kh2o+dkh2o*q) has been used to determine kh2o and dkh2o
	    if (nwvpod_ > 1) THEN
	       print*, 'Err[oss_ir_module::GetOD]: Versions with NWVPOD_ >= 2 will be released later '
	       call exit(1)
	    end if
	    nlev=nlayod+1
	    allocate (pref(nlev),prefLoc(nlev))
	    allocate (tmpTab(ntmpod,nlayod),wvpTab(nwvpod_,nlayod))
	    allocate (tmpLoc(ntmpod,nlev),wvpLoc(nwvpod_,nlayod))
	    allocate (fixDMR(nlayod),TmpSelf(ntself))
	    allocate (dfltLoc(nmol,nlayod),fixLoc(nlayod))
	    allocate (fix(nlayod))
	    allocate (pavlref(nlayod))
	    allocate (vwvn(nfsmp_ir),f1_arr(nfsmp_ir),f2_arr(nfsmp_ir))
	    !wvpTab, mixing ratios for WV
	    read(iuo)prefLoc(1:nlev),tmpLoc(1:ntmpod,1:nlayod),TmpSelf(1:ntself),wvpLoc(1:nwvpod_,1:nlayod)
	    pref(1:nlev)=prefLoc(1:nlev)
	    tmptab(1:ntmpod,1:nlayod)=tmpLoc(1:ntmpod,1:nlayod)
	    wvptab(1:nwvpod_,1:nlayod)=wvpLoc(1:nwvpod_,1:nlayod)

	    pavlref(1:nlayod)=(pref(2:nlayod+1)-pref(1:nlayod))/log(pref(2:nlayod+1)/pref(1:nlayod))
	    !Transform level reference temperatures onto a layer grid
	    read(iuo)fixLoc(1:nlayod),dfltLoc(1:nmol,1:nlayod)

	    fixDMR(1:nlayod)=fixLoc(1:nlayod)
	    allocate (invMolID(molID(nmol)))
	!in order to count variable species interms of compressed indices, we set
	!auxiliary array invMolID

	    ! Continue reading defProfFile
	    call invertMolID(nmol,molID(1:nmol),invMolID)
		    !go to specific quantities

            allocate (molProf(nmol,nlayOD+1))
	    
            !- Allocate temp array for default profiles
	    allocate(vmol(nlayod+1))

	    k = 0
	    do m=1,nvarmol
	         call read_mol_prof(iuF,varmolid(m),vmol)
		 if ( ANY(molid(1:nmol) == varmolid(m))) then 
			 !- Removed vmol integration on pref
			 k = k + 1
	     		 molProf(k,1:nlayod+1) = vmol
		 end if
	    end do
    	    CLOSE (iuf)
	    deallocate(vmol)

            ! <---- Cut by Psolo S.

	    !---read oss parameters
	    allocate (coef_ir(nchmax,nf_sel))
	    allocate (ichmap_ir(nchmax,nf_sel))
	    allocate (iselS(nf_sel),nch_ir(nf_sel))
	    allocate (coef_Loc(nchmax,nf_sel))
	    allocate (ichmap_Loc(nchmax,nf_sel))
	    allocate (isel_Loc(nf_sel),nch_Loc(nf_sel))
	    do ismp=1,nf_sel
	       read(ius)isel_Loc(ismp),nch_Loc(ismp)
	       read(ius)coef_Loc(1:nch_Loc(ismp),ismp),ichmap_Loc(1:nch_Loc(ismp),ismp)
	    end do
	    CLOSE(ius)
	    allocate (nNodes_arr(0:mxIndex))
	    allocate (iselS_arr(nf_sel, 0:mxIndex))
	    allocate (nch_arr(nf_sel, 0:mxIndex))
	    allocate (ichMap_arr(nchmax,nf_sel, 0:mxIndex))
	    allocate (coef_arr(nchmax,nf_sel, 0:mxIndex))
	    allocate (nChList_arr(0:mxIndex))
	    allocate (chanList_arr(nChan, 0:mxIndex))

	    lastChSet=-1
	    nchanAll = nChan
	    chSetID=0
            
	    call loadChanSelect(chIndx, nChan, chSetID)
	    call setChanSelect(chSetID)
	   
	    !---read absorption coefficient tables
	    nsize1=nlayod*ntmpod
	    nsize2=nsize1*nwvpod
	    allocate (kfix_ir(nsize1,nfsmp_ir),kh2o_ir(nsize2,nfsmp_ir))
	    allocate (kvar_ir(mxmols*nsize1,nfsmp_ir),kself(ntself,nfsmp_ir))
	!xid includes 0 (if fixed gases are used as dependent) and OSS indices for
	!other variable species; xid_ uses compressed indices instead
	!The following cases are under consideration: 1. Fixed gases are in xid, WV is
	!also in, other variable species could, or, could not be in. 2. Fixed gases are
	!not in xid, however, WV is in, other variable species could, or, could not be in.
	    allocate (xid_(nx))
	    allocate (dkfix_ir(MAX(1,nsize1*(1-xid(1))),nfsmp_ir))
	    allocate (dkh2o_ir(nsize1,nfsmp_ir))

	    if (xid(1) == 0) then
	       xid_(1) = 0
	       xid_(2:nx)=invMolID(xid(2:nx))
	       if ( ANY(xid_(2:nx) == 0)) THEN
		  print*, 'Err[oss_ir_module::GetOD]: xid contains indices out of molID'
		  call exit(1)
	       end if
	    end if
	    if (xid(1) == 1) then
	       xid_(1:nx) = invMolID(xid(1:nx))
	       if ( ANY(xid_(1:nx) == 0)) THEN
		  print*, 'Err[oss_ir_module::GetOD]: xid contains indices out of molID'
		  call exit(1)
	       end if
	    end if	
	    nsize4=nsize1*(nx-2+xid(1))
	    nsize5=nsize1*(1-xid(1))
	    allocate (kfix_Loc(nsize1),dkh2o_Loc(nsize1),kh2o_Loc(nsize2),kself_Loc(ntself))
	    allocate (dkfix_Loc(MAX(1,nsize5)))

	    allocate (imols_output(nfsmp_ir,mxhmol))   !- imols_ array to be written inside the netcdf
	    allocate (nmols_output(nfsmp_ir))		!- nmols_ array to be written inside the netcdf
	    nmols_output = 0                            !- init nmols array

	!Transformation of xid from molID indices to compressed indices
	    do ISmp=1,Nfsmp_ir

	       imols_output(ismp,1:mxhmol) = 0	       !- init imols_ array (output)

	       read(iuo)vwvn(ismp),nmols_ !real*8, int*4
	       if (Spc_Units_l == 1) vWvn(ismp)= vWvn(ismp)/SOfL!From GHz 2 wn
	       if (Spc_Units_l == 2) vWvn(ismp)= 1.e4/vWvn(ismp)!From mu 2 wn
	       f1_arr(ismp)=f1*(vwvn(ismp)*1.d-3)**3
	       f2_arr(ismp)=f2*(vwvn(ismp)*1.d-3)
	       read(iuo)Imols_(1:nmols_) !int*4
	       ! transformation to int*2 is done for using nmols_Tmp and imols_tmp
	       ! in odtresh, cum_fix and shrink_var
	       nmols_Tmp = nmols_
	       Imols_Tmp(1:nmols_)=Imols_(1:nmols_)

	       !- Fill output arrays
	       nmols_output(ismp) = nmols_
	       imols_output(ismp,1:nmols_) = imols_(1:nmols_)

	       nsize3=nsize1*(nmols_-1)
	       allocate (kvar_Loc(max(1,nsize3)))
	       read(iuo)kfix_Loc(1:nsize1),kh2o_Loc(1:nsize2),kvar_Loc(1:nsize3),kself_Loc(1:ntself)
	       kfix_ir(1:nsize1,ismp)=kfix_Loc(1:nsize1)

	       kh2o_ir(1:nsize2,ismp)=kh2o_Loc(1:nsize2)
	       kvar_ir(1:nsize3,ismp)=kvar_Loc(1:nsize3)
	       kself(1:ntself,ismp)=kself_Loc(1:ntself)
	       read(iuo)dkfix_Loc(1:nsize5),dkh2o_Loc(1:nsize1)
	       !---code uses kfix_ir*rfix instead of kfix_ir
	       dkfix_ir(1:nsize5,ismp)=dkfix_Loc(1:nsize5)
	       dkh2o_ir(1:nsize1,ismp)=dkh2o_Loc(1:nsize1)
               if (isHDO) then
                     ! skip HDO data
                     read(iuo)kfix_Loc(1:nsize1), dkh2o_Loc(1:nsize1)
               end if
              
               if (allocated(kvar_Loc)) deallocate(kvar_Loc)
            end do	       
	    CLOSE(iuo)
		
	    fix(1:nlayod) = fixLoc(1:nlayOD)
	    deallocate (kfix_Loc,kh2o_Loc,dkfix_Loc,dkh2o_Loc)
	    deallocate (dkfix_ir, fixDMR)
	    deallocate(fixLoc)
	    return
	  end subroutine GetOD
	  
	    subroutine purpose( )
	      implicit none
	      write(stderr,*) 'Converts HITRAN pre-computed Look Up Tables in netCDF'
	    end subroutine purpose

	    subroutine usage( )
	      implicit none
	      character(len=256) :: myname
	      call get_command_argument(0,myname)
	      write(stderr,*) 'Usage: '
	      write(stderr,*) '       ',trim(myname),' basis'
	      write(stderr,*)
	      write(stderr,*) 'Example: ',trim(myname),' leo.airs.0.05'
	      write(stderr,*)
	      stop
	    end subroutine usage

	    subroutine check_ncerror(message)
	      implicit none
	      character(len=*) , intent(in) :: message
	      if ( ncstat /= nf90_noerr ) then
		write(stderr,*) message
		write(stderr,*) nf90_strerror(ncstat)
		stop
	      end if
	    end subroutine check_ncerror
  
	!----------------------------------------------------------------------------
	! PURPOSE: to load channel subset for further use 
	! input
	!   iChList - subset indices (has to be the same format as in selection file)
	!   nChList - length of the subset
	!   chSetID - user supplied ID to the subset
	!----------------------------------------------------------------------------
	subroutine loadChanSelect(iChList, nChList, chSetID)
	   integer, dimension(:), intent(in)    :: iChList
	   integer,               intent(in)    :: nChList   
	   integer,               intent(in)    :: chSetID

	   integer                :: k,n1,ich,ismp,ismpSel, idx, lenBase
	   integer                ::curChSet
	   
	   if (lastChSet .ge. 0) then
	      if (chSetID .le. 0) then
		print*, 'Wrn[oss_ir_module:: loadChanSelect]: chSetID must be a positive integer.',  &
		   ' User subset with index  ', chSetID, ' will no be available '
		return   
	      end if
	       
	      curChSet = userIndex2chSetIndex(chSetID)
	      if (curChSet .eq. -1) then
		 if (lastChSet .ge. mxIndex) then
		    print*, 'Wrn[oss_ir_module:: loadChanSelect]: maximum number of channel subset has been reached.',  &
		       ' Subset with index  ', chSetID, ' will no be available '
		    return   
		 end if
		 lastChSet = lastChSet+1
		 curChSet = lastChSet
	      end if
	      lenBase = nChList_arr(0)
	      do k=1,nChList
		 if ( .not. ANY(iChList(k) .eq. chanList_arr(1:lenBase,0))) THEN
		    print*, 'Err[oss_ir_module:: loadChanSelect]: In subset ',  lastChSet, &
		       ' some channel indices do not present in the SEL file '
		    call exit(1)
		 end if
	      end do 
	      do n1=1,nChList-1
		 if ( ANY(iChList(n1+1:nChList) == iChList(n1))) THEN
		    print*, 'Err[oss_ir_module:: loadChanSelect]: There are duplicates '//&
		                'in selected channels '
		    call exit(1)
		 end if
	      end do
	   else
	      lastChSet = 0
	      curChSet = 0
	   end if
	   userIndex2chSetIndexMap(curChSet) = chSetID
	    

	    ismpSel=0
	    
	    if (nChList == nchanAll) then
	       nNodes_arr(curChSet)=nf_sel
	       iselS_arr(1:nf_sel,curChSet)=isel_Loc(1:nf_sel)
	       nch_arr(1:nf_sel,curChSet)=nch_Loc(1:nf_sel)
	       do ismp=1,nf_sel 
		  do k=1,nch_Loc(ismp)
		     ichMap_arr(k,ismp,curChSet)= &
		         ichMap_Loc(k,ismp)
		     coef_arr(k,ismp,curChSet)= &
		         coef_Loc(k,ismp)
		  end do
	       end do
	    else
	       do ismp=1,nf_sel
		  ismpSel=ismpSel+1
		  ich=0
		  do k=1,nch_Loc(ismp)
		     do idx = 1, nChList
		        if (iChList(idx) .EQ. ichmap_Loc(k,ismp))  THEN
		            ich=ich+1
		            ichMap_arr(ich,ismpSel,curChSet)=idx
		            coef_arr(ich,ismpSel,curChSet)=coef_Loc(k,ismp)
		            EXIT
		         end if   
		     end do    
		  end do
		  if (ich == 0) THEN
		     ismpSel=ismpSel-1
		  else
		     nch_arr(ismpSel,curChSet)=ich
		     iselS_arr(ismpSel,curChSet)=isel_Loc(ismp)
		  end if
		  nNodes_arr(curChSet)=ismpSel
	       end do
	    end if

	    nChList_arr(curChSet)=nChList
	    chanList_arr(1:nChList,curChSet)=iChList(1:nChList)

	    return
	  end subroutine loadChanSelect

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
	    iselS(1:nNodes)=iselS_arr(1:nNodes,iChSet)
	    nch_ir(1:nNodes)=nch_arr(1:nNodes,iChSet)
	    do k=1,nNodes
	       do ich=1,nch_ir(k)
		  ichMap_ir(ich,k)=ichMap_arr(ich,k,iChSet)
		  coef_ir(ich,k)=coef_arr(ich,k,iChSet)
	       end do
	    end do
	    return
	  end subroutine setChanSelect


    ! PURPOSE: Create auxiliary array, mapping molID onto a position within 
    !          array of absorption coefficients
    !----------------------------------------------------------------------------
    subroutine invertMolID(nmol,molID,invMolID)
      !Input variables
      integer,     intent(in) :: nMol,molID(nMol)
      !Output variables
      integer,    intent(inout) :: invMolID(molID(nmol))
       !Local variables
      integer k
      invMolID(:)=0
    
      do k=1,nmol
         invMolID(molID(k))=k
      end do
      return
    end subroutine invertMolID

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
  
!---------------------------------------------------------------------------------------------
! The subroutine odthresh finds species (water vapor excluded)
! which contribution to optical depth is less than odfac of total
! and set flag in array iflag

  subroutine odthresh(kfix0,kvar0,imols,nmols,w,odfac,iflag)
    !---Input variables
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


!---------------------------------------------------------------------------------------------
! The subroutine shrink_var reduce the size of kvar 
! and makes a new map of molecular indices 
  
  subroutine shrink_var(kvar0,imols0,nmols0,maps,imols_indx,kvar1,imols1,nmols1)
    !---In/Out variables
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

  !----------------------------------------------------------------------------
  ! PURPOSE: Computes average layer quantities (or integrated amount) 
  !          using a log-x dependence on log-p.
  !----------------------------------------------------------------------------
  subroutine lpsum_log(pu,pl,xu,xl,scal,xint,dxu,dxl)
    real, parameter :: epsiln=1.e-12
    !---Input variables
    real, intent(in)  :: pu,pl,xu,xl,scal
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


  subroutine ossDestroy_ir()
    if (allocated(nch_ir))        deallocate(nch_ir)
    if (allocated(vWvn))          deallocate(vWvn)
    if (allocated(f1_arr))        deallocate(f1_arr)
    if (allocated(f2_arr))        deallocate(f2_arr)
    if (allocated(coef_ir))       deallocate(coef_ir)
    if (allocated(cWvn))          deallocate(cWvn)
    if (allocated(pref))          deallocate(pref)
    if (allocated(pavlref))       deallocate(pavlref)
    if (allocated(sunrad))        deallocate(sunrad)
    if (allocated(ichmap_ir))     deallocate(ichmap_ir)
    if (allocated(Tmptab))        deallocate(Tmptab)
    
    if (allocated(kfix_ir))       deallocate(kfix_ir)
    if (allocated(kh2o_ir))       deallocate(kh2o_ir)
    if (allocated(dkh2o_ir))      deallocate(dkh2o_ir)
    if (allocated(kvar_ir))       deallocate(kvar_ir)
    if (allocated(kself))         deallocate(kself)
    if (allocated(molid))         deallocate(molid)
    if (allocated(nChList_arr))   deallocate(nChList_arr)
    if (allocated(chanList_arr))  deallocate(chanList_arr)
    if (allocated(iselS))         deallocate(iselS)
    if (allocated(nNodes_arr))    deallocate(nNodes_arr)
    if (allocated(iselS_arr))     deallocate(iselS_arr)
    if (allocated(isel_Loc))      deallocate(isel_Loc)
    if (allocated(nch_Loc))       deallocate(nch_Loc)
    if (allocated(nch_arr))       deallocate(nch_arr)
    if (allocated(coef_arr))      deallocate(coef_arr)
    if (allocated(coef_Loc))      deallocate(coef_Loc)
    if (allocated(ichMap_arr))    deallocate(ichMap_arr)
    if (allocated(ichmap_Loc))    deallocate(ichmap_Loc)
    
   if (allocated(xG))             deallocate (xG)
   if (allocated(surfEmRfGrid))   deallocate (surfEmRfGrid)
   if (allocated(surfEmRefl))     deallocate(surfEmRefl)
   if (allocated(chanFreq))       deallocate(chanFreq)
   if (allocated(chanIndex))      deallocate(chanIndex)
   if (allocated(surfEmRfGrid))   deallocate(surfEmRfGrid)
   if (allocated(molProf))        deallocate(molProf)
   if (allocated(sunrad))         deallocate(sunrad)
   if (allocated(fix))            deallocate(fix)
   if (allocated(wvptab))         deallocate(wvptab)
   if (allocated(dfltLoc))        deallocate(dfltLoc)
   if (allocated(nmols_output))    deallocate(nmols_output)
   if (allocated(imols_output))   deallocate(imols_output)
   if (allocated(TmpSelf))        deallocate(TmpSelf)
   
    iChSet                = UNDEFINED_INDEX
    userIndex2chSetIndexMap(1:mxIndex) = -1
    userIndex2chSetIndexMap(0)         = 0
  end subroutine ossDestroy_ir

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

end program odcmpr2nc

