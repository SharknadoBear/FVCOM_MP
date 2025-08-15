!/===========================================================================/
! Copyright (c) 2007, The University of Massachusetts Dartmouth 
! Produced at the School of Marine Science & Technology 
! Marine Ecosystem Dynamics Modeling group
! All rights reserved.
!
! FVCOM has been developed by the joint UMASSD-WHOI research team. For 
! details of authorship and attribution of credit please see the FVCOM
! technical manual or contact the MEDM group.
!
! 
! This file is part of FVCOM. For details, see http://fvcom.smast.umassd.edu 
! The full copyright notice is contained in the file COPYRIGHT located in the 
! root directory of the FVCOM code. This original header must be maintained
! in all distributed versions.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
! AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO,
! THE IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR
! PURPOSE ARE DISCLAIMED.  
!
!/---------------------------------------------------------------------------/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/

!============================================================================/
! FVCOM Plastic Module
!
! Copyright:    2022(c)
!
! THIS IS A DEMONSTRATION RELEASE. THE AUTHOR(S) MAKE NO REPRESENTATION
! ABOUT THE SUITABILITY OF THIS SOFTWARE FOR ANY OTHER PURPOSE. IT IS
! PROVIDED "AS IS" WITHOUT EXPRESSED OR IMPLIED WARRANTY.
!
! THIS ORIGINAL HEADER MUST BE MAINTAINED IN ALL DISTRIBUTED
! VERSIONS.
!
! Contact:      Y.Huang, Pacific Northwest National Laboratory
!
! Plastic transport model developed, based on the structure of 
! Community Sediment Transport Model (CSTM) as implemented in ROMS 
! by J. Warner (USGS)
!
! Comments:     Plastics Dynamics Module 
!
! FVCOM dependency
!
!   init_sed.F:  - user defined sediment/plastic model initial conditions
!   mod_ncdio.F: - netcdf output includes concentration/bottom/bed fields
!	mod_bbl.F	 - bottom boundary layer module
!   fvcom.F:     - main calls sediment/plastic  setup 
!   internal_step.F:  - calls sediment advance
!
! History
!   Dec 30,2022: merge the whole module from the original mod_sed.F, delete
!               all the redundent part, keep the essential features for a beta 
!				version that could be used for testing the new plastic settling
!				module. Still using the sediment module and most of the 
!				variable names for beta version compiling and testing, will
!				be isolated from the CSTMS coding in the future modification. 
!				Fluid_Mud interaction and concentration depended fluculation
!				and settling velocity variation is completly removed in this
!				version. The wave_current_interaction part is kept, but not 
!				functional. The spatial varying bed feature is not available 
!				as well. The water density in the settling module is depending 
!				on the ambient water density, which could be affected by the 
!				salinity and temperature.
!   Mar 17,2023: Now the module is still compiled as a sediment module, but 
!				the module is almost isolated from the original sediment module. 
!				Things below 'implicit none' are almost completely changed with a plastic 
!				variable name. Although the river forcing is still not isolated from
!				sediment, it is still using N_sed/N_sed_max for setup the forcing. All the 
!               function within this module is adapted for plastic. Wave_current_interaction 
!				part are temporary deleted. Bed variables are completely removed, 
!               asuming the plastic is always at the bottom interface. Additional module 
!				considering the infiltration will be included in the future.
!				Morphology part is completely deleted. Bottom feature is all 
!				about the sediment. Resuspension module is added-in with WS2019 
!				treatment on the critical shear stress, depending on the ambient 
!				density. The plastic on the bottom is resuspended as sediment, 
!				but with the amount related to the sediment bed feature (active-layer-
!				thickness, mud-content, sediment density and porosity). 
!
!==========================================================================/
Module Mod_Sed  

#if defined (SEDIMENT) && (ORIG_SED)
	Use Mod_Par
	Use Mod_Prec 
	Use Mod_Types
	Use Mod_wd
	Use Control, only : seddis
	Use all_vars,only : CNSTNT,UNIFORM

	implicit none 

	!--------------------------------------------------
	!Plastic Type                         
	!
	! pname        => plastic name
	! rtype		   => ambient density field type 	
	! aaxis		   => the longest axis of grain
	! baxis		   => the median axis of grain	   
	! caxis	 	   => the shortest axis of grain
	! cmax         => the maximum value of plastic concentration at each step
	! cmin		   => the minimum value of plastic concentration at each step
	! crms		   => the rms value of plastic concentration at each step
	! dequ         => equivalent grain size 					[mm] 
	!				  (for plastic it will be estimated by the equivalent diameter a/b/caxis)
	! wset		   => spatial varying plastics settling velocity[m/s]
	! tau_ce       => critical shear stress for erosion 		[N/m^2]]
	! Prho         => plastic density 							[kg/m^3]
	! erate        => bed surface erosion mass flux         	[kg m^-1 s^-1] 
	! conc         => plastics concentration in water column 	[kg m^-3 or g/L]
	! cnew         => plastics concentration during update   	[kg m^-3 or g/L]
	! mass         => plastics mass density in at the bottom    [kg m^-2]
	! frac         => plastics fraction in active layer 	    [-]
	! eflx         => suspended plastics erosive flux (+out of bed) [kg m^-2]
	! dflx         => suspended plastics depositional flux (+into bed) [kg m^-2]
	! cdis         => concentration at river source 
	! cflx         => store advective flux at open boundry 
	! cobc         => user specd open bndry concentration 
	! arraysize    => spatial size of sediment arrays
	!--------------------------------------------------

	type plast_type
		character(len=20) :: pname     
		real(sp)		  :: aaxis
		real(sp)		  :: baxis
		real(sp)		  :: caxis
		real(dp)          :: cmax 
		real(dp)          :: cmin 
		real(dp)          :: crms
		real(sp)		  :: dequ 		
		real(sp)		  :: Prho
		real(sp), allocatable :: conc(:,:)
		real(sp), allocatable :: cnew(:,:)
		real(sp), allocatable :: mass(:)
		real(sp), allocatable :: frac(:)
        real(sp), allocatable :: eflx(:)        
		real(sp), allocatable :: dflx(:)        
		real(sp), allocatable :: cdis(:)
		real(sp), allocatable :: cflx(:,:)		
		real(sp), allocatable :: cobc(:)
		real(sp), allocatable :: wset(:,:) 
		real(sp), allocatable :: zactv(:)
		real(sp), allocatable :: erate(:)
		real(sp), allocatable :: tau_ce(:)	
		integer           :: arraysize
	end type
 
	public

	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	!Part II: Public variables
	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	!--------------------------------------------------------------------------
	!Global Model Parameters
	!
	! plastic_model_version : Model version label
	! plastfile   	    : sediment input parameter control file
	! inf_bed   		: switch of infinite plastic supply
	! oned_model		: switch of 1D plastic module
	! plast_nudge		: switch of plastic OBC nudging
	! plast_source		: switch of plastic point sources
	! PLAST_hot_start	: switch of plastic restart
	! nplast      		: number of plastic classes
	! n_report  		: iteration interval for running report
	! Nobc				: number of OBC nodes
	! Nobc_gl           : global number of OBC nodes
	! PLAST_ramp  		: number of iterations over which to ramp plastic nudging
	! PLAST_alpha 		: plastic nudging relaxation factor
	! PLAST_start 		: start interval for plastic model
	! PLAST_its   		: plastic model iteration counter
	! settle_limiter	: settling process FVM limiter type selection (0,1,2)
	!					 0: is using 1st order upwind
	!					 1: is using minmod limiter
	!					 2: is using Van Leer limiter
	! DTsed     		: plastic model time step (seconds)
	! qlim				: limiter parameter
	! T_model   		: model time (seconds)
	! taub      		: array to hold bottom shear stress
	! Wrho0 			: constant water density
	! mass0 			: constant initial plastic mass density in the bottom
	! sd50_0            : constant median grain size in the bottom 
	! dens0				: constant sediment density in the bottom
	! mudp0             : constant mud portion in the bottom
	! poro0             : constant porosity in the bottom
	! rho_water			: seawater density can be constant or spatial varying
	!-------------------------------------------------------------------------

	character(len=120), parameter :: PLAST_model_version = "cptms_v0.3: profiling"
	character(len=120) :: plastfile
	logical  :: inf_bed
	logical  :: oned_model
	logical  :: plast_nudge
	logical  :: PLAST_dumpbed
	logical  :: PLAST_dumpbot 
	logical  :: PLAST_source
	logical  :: PLAST_hot_start
	integer  :: nplast
	integer  :: n_report
	integer  :: Nobc
	integer  :: Nobc_gl
	integer  :: PLAST_ramp  
	integer  :: PLAST_start
	integer  :: PLAST_its
	integer  :: settle_limiter = 2
	real(sp) :: DTsed
	real(sp) :: qlim  = 1.0
	real(sp) :: PLAST_alpha
	real(sp) :: settle_cfl = 1.0
	real(sp) :: T_model
	real(sp) :: Wrho0 = 1025
	real(sp) :: mass0
	real(sp) :: sd50_0
	real(sp) :: dens0 = 2400
	real(sp) :: mudp0 = 0.1
	real(sp) :: poro0 = 0.5
	real(sp), allocatable, target :: rho_water(:,:)
	real(sp), allocatable, target :: taub(:)

	!----------------------------------------------------------------------!
	!  BOTTOM properties indices:                                          !
	!																	   !
	!  data --> bottom(horizontal,prototype)						       !
	!  size --> bottom(0:mt,n_bot_vars)   				           		   !
	!																	   !
	! with proptype indices:											   !
	!   isd50  =>  Median sediment grain diameter.				 (m)       !
	!   idens  =>  Median sediment grain density. 				 (kg/m3)   !
	!   imudp  =>  Mud Potion.  					  			 (-)       !
	!   iporo  =>  Sediment porosity 							 (-)       !
	!																	   !
	! global bed characteristics:										   !
	!   n_bot_vars:  number of bed characteristics						   !
	!	bot_snames:   bottom characteristics name list					   !
	!	bot_lnames:   bottom characteristics long name list				   !
	!	bot_units: 	  bottom characteristics unit list                     !
	!----------------------------------------------------------------------!

	integer, parameter    :: n_bot_vars  = 4
	character(len=80), dimension(n_bot_vars) :: bot_snames = &
	(/"bot_sd50","bot_dens","bot_mudp","bot_poro"/)
	character(len=80), dimension(n_bot_vars) :: bot_lnames = &
	(/"bot_sd50","bot_dens","bot_mudp","bot_poro"/)
	character(len=80), dimension(n_bot_vars) :: bot_units = &
	(/"m","kg/m^3","dimensionless","dimensionless"/)

	integer, parameter    :: isd50  = 1  
	integer, parameter    :: idens  = 2
	integer, parameter    :: imudp  = 3
	integer, parameter    :: iporo  = 4 

	real(sp), allocatable, target :: bottom(:,:)

	!--------------------------------------------------
	!Plastic Array                  
	!--------------------------------------------------

	type(PLAST_type), allocatable, target :: plast(:)

	contains

!==========================================================================
! Allocate Data, Initialize Concentration Fields and Bed Parameters       
!==========================================================================
  
 Subroutine Setup_PLAST(fname,Nin,Nin_gl,N_sed,N_sed_Max,Sed_Names) 
	
	use control, only : msr,casename,input_dir
	use lims,    only : m, kbm1
	use Mod_Utils, only: pstop
	
	implicit none
	
	!Input/Output variables
	integer, intent(in)  :: Nin,Nin_gl   
	character(len=80)    :: fname
	integer, intent(out) :: N_sed
	integer, intent(in ) :: N_sed_max   
	character(len=20)    :: Sed_Names(N_sed_max)
	
	!Local variables
	logical             :: fexist
	integer             :: iplast
	integer             :: i,k

	if(dbg_set(dbg_sbr)) write(ipt,*) "Start: setup_plast" 

	Nobc = Nin
	Nobc_gl = Nin_gl
	
	!----------------------------------------------------------------
	! check arguments and ensure necessary files exist
	!----------------------------------------------------------------

	!initialize plastic model iteration counter
	PLAST_its = 0

	!ensure plastic parameter file exists
	if(fname == "none" .or. fname == "NONE" .or. fname == "None")then
		plastfile = "./"//trim(input_dir)//"/"//trim(casename)//'_plastic.inp'
	else
		plastfile = "./"//trim(input_dir)//"/"//trim(fname)
	endif
    
	inquire(file=trim(plastfile),exist=fexist)
	if(.not.fexist)then
		write(*,*)'plastic parameter file: ',trim(plastfile),' does not exist'
		write(*,*)'stopping'
		call pstop
	end if

	!----------------------------------------------------------------
	! read in plastic parameter file                    
	!----------------------------------------------------------------
		call read_PLAST_params
  
	!----------------------------------------------------------------
	! allocate data                                      
	!----------------------------------------------------------------
		call alloc_PLAST_vars

	!----------------------------------------------------------------
	! initialize fields                                  
	!----------------------------------------------------------------
		call init_Plast

	!----------------------------------------------------------------
	! setup open boundary condition nudging values       
	!----------------------------------------------------------------
		if(plast_nudge) call setup_plast_obc 
		
	!----------------------------------------------------------------
	! report sediment parameters and statistics to screen
	!----------------------------------------------------------------
		call report_PLAST_setup

	!----------------------------------------------------------------
	! convert parameters to MKS unit system 
	!----------------------------------------------------------------
		call convert_PLAST_params 

	!----------------------------------------------------------------
	! set arguments for plastic names and N_sed for FVCOM to use 
	! in setting up river forcing, now the N_sed is still unchanged
	! because the river forcing is not isolated from sediment, will be 
	! updated in the later version, when plastic module is compile 
	! isolatedly by itself.
	!----------------------------------------------------------------
		N_sed = nplast
		do i=1,nplast
			Sed_names(i) = plast(i)%pname
		end do
		
		if(dbg_set(dbg_sbr)) write(ipt,*) "End: setup_plast" 

end Subroutine Setup_PLAST

!=======================================================================
!Convert Sed Params from nonstandard (mm, etc) to mks units
!=======================================================================
Subroutine Convert_PLAST_Params 
	
	Implicit None
	integer :: i

	if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Convert_PLAST_Params" 

	!convert mean grain diameter from mm to m 
	do i=1,nplast
		plast(i)%dequ = plast(i)%dequ*.001
		plast(i)%aaxis= plast(i)%aaxis*.001
		plast(i)%baxis= plast(i)%baxis*.001
		plast(i)%caxis= plast(i)%caxis*.001
	end do
	
	bottom(:,isd50)=bottom(:,isd50)*.001

	if(dbg_set(dbg_sbr)) write(ipt,*) "End: Convert_PLAST_Params" 

End Subroutine Convert_PLAST_Params

!=======================================================================
!Read Sediment Parameters From Sediment Input File
!=======================================================================
Subroutine Read_PLAST_Params
  
  Use Input_Util
  
  Implicit None
  
  integer linenum,i,k1,iscan
  real(sp)           :: ftemp
  character(len=120) :: stemp

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Read_PLAST_Params " 

  !check if we are running 1D model
  Call Get_Val(oned_model,plastfile,'PLAST_ONED',line=linenum,echo=.false.)

  !read in number of plastic classes 
  Call Get_Val(nplast,plastfile,'NPLAST',line=linenum,echo=.false.)  

  !read in start interval for plastic model 
  Call Get_Val(PLAST_start,plastfile,'PLAST_START',line=linenum,echo=.false.)  

  !hot start option
  Call Get_Val(PLAST_hot_start,plastfile,'PLAST_HOT_START',line=linenum,echo=.false.)  

  !read in interation interval for reporting 
  Call Get_Val(n_report,plastfile,'N_REPORT',line=linenum,echo=.false.)  

  !read in logical switch for infinite plastic supply 
  Call Get_Val(inf_bed,plastfile,'INF_PLAST',line=linenum,echo=.false.)

  !read in nudging switch 
  Call Get_Val(plast_nudge,plastfile,'PLAST_NUDGE',line=linenum,echo=.false.)

  !read in params controlling output 
  Call Get_Val(PLAST_dumpbot,plastfile,'PLAST_DUMPBOT',line=linenum,echo=.false.)

  !read in nudging relaxation factor 
  if(plast_nudge) then
  Call Get_Val(PLAST_alpha,plastfile,'PLAST_ALPHA',line=linenum,echo=.false.)
  Call Get_Val(PLAST_ramp ,plastfile,'PLAST_RAMP ',line=linenum,echo=.false.)
  endif

  !settling CFL
  Call Get_Val(settle_cfl ,plastfile,'SETTLE_CFL',line=linenum,echo=.false.)
  if(settle_cfl <= 0 .or. settle_cfl > 1.)then
    write(*,*)'error in mod_sed, settle_cfl must be in range [>0,1]'
    write(*,*)'you entered: ',settle_cfl
    stop
  endif

  !limiter type
  Call Get_Val(settle_limiter,plastfile,'SETTLE_LIMITER',line=linenum,echo=.false.)
  if(settle_limiter < 1 .or. settle_limiter > 3)then
    write(*,*)'error in mod_sed, settle_limiter must be [1,2,3]'
    write(*,*)'you entered: ',settle_limiter
    stop
  endif
  if(settle_limiter==1) qlim = 0.
  if(settle_limiter==2) qlim = 1.
  if(settle_limiter==3) qlim = 2.

  !read in point source switch 
  Call Get_Val(PLAST_source,plastfile,'PLAST_PTSOURCE',line=linenum,echo=.false.)
  
  !read in the constant ambient water density
  Call Get_Val(Wrho0,plastfile,'AMBIENT_WDEN',line=linenum,echo=.false.)
  
  !read in the constant initial mass variation
  Call Get_Val(mass0,plastfile,'INI_MASS',line=linenum,echo=.false.)
  
  !read in the constant median sediment diameter variation
  Call Get_Val(sd50_0,plastfile,'SED_SD50',line=linenum,echo=.false.)
  
  !read in the constant sediment density
  Call Get_Val(dens0,plastfile,'SED_SRHO',line=linenum,echo=.false.)
  
  !read in the constant mud portion
  Call Get_Val(mudp0,plastfile,'SED_MUDP',line=linenum,echo=.false.)
  
  !read in the constant sediment porosity
  Call Get_Val(poro0,plastfile,'SED_PORS',line=linenum,echo=.false.)

  !check values
  if(nplast < 1)then
    write(*,*)'sediment input file must have at least one sediment class'
    write(*,*)'currently nplast = ',nplast
    stop
  endif
  
  !allocate plastic data space
  Allocate(plast(nplast))

	!read in plastic parameters 
	!(this part only read-in the constants, for spatial varying value they are input by mod_input)
	k1 = 1

	do i=1,nplast

		!read PLAST_NAME and mark position
		Call Get_Val(stemp,plastfile,'PLAST_NAME',line=linenum,echo=.false.,start=k1)
		plast(i)%pname = stemp
		k1        = linenum+1

		!read longest plastic axis
		Call Get_Val(ftemp,plastfile,'PLAST_AAXI',line=linenum,echo=.false.,start=k1)
		plast(i)%aaxis = ftemp
		
		!read medium plastic axis
		Call Get_Val(ftemp,plastfile,'PLAST_BAXI',line=linenum,echo=.false.,start=k1)
		plast(i)%baxis = ftemp
		
		!read shortest plastic axis
		Call Get_Val(ftemp,plastfile,'PLAST_CAXI',line=linenum,echo=.false.,start=k1)
		plast(i)%caxis = ftemp

	    !read plastic density 
		Call Get_Val(ftemp,plastfile,'PLAST_PRHO',line=linenum,echo=.false.,start=k1)
		plast(i)%Prho = ftemp 

		!calculate or read mean diameter
		plast(i)%dequ = (plast(i)%aaxis*plast(i)%baxis*plast(i)%caxis)**(0.3333)
		
	end do
 
 if(dbg_set(dbg_sbr)) write(ipt,*) "End: Read_PLAST_Params " 
 
 End Subroutine Read_PLAST_Params 
     
!=======================================================================
! Advance Sediment Model by time step DTSED 
!=======================================================================
Subroutine Advance_plast(DTin,Tin,taub_in,s_in,t_in)  
  
  Use Input_Util
  Use Scalar
  Use ALL_VARS
  Use Control, only : ireport,iint,msr,par
  Use Lims,    only : m,mt,nt,kbm1,numqbc,kb,nprocs,myid
	
# if defined (MULTIPROCESSOR)
  Use Mod_Par
# endif

  Implicit None
  real(sp), intent(in ) :: DTin,Tin
  real(sp), intent(in ) :: taub_in(m),s_in(m,kb),t_in(m,kb)

  integer :: i,k,iplast,l1,l2,ierr,d_cdis,d_cflx
  character(len=4) :: fnum
  real(sp) :: fact,ufact,temp,s_temp(kbm1),t_temp(kbm1),rho_temp(kbm1)

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Advance_Plast "

  !------------------------------------------------------
  !Check for Initialization
  !------------------------------------------------------
  if(iint < PLAST_start) return

  !------------------------------------------------------
  !Increment Plastic Model Counter and Report Plastic Model Init 
  !------------------------------------------------------
  PLAST_its = PLAST_its + 1
  if(msr .and.  PLAST_its == 1)then
    write(*,*)'========Plast Model Initiated======='
  endif
  
  !------------------------------------------------------
  !read-in the density variation from baroclinic time-step 
  !------------------------------------------------------
  if (AMBIENT_DENSITY_TYPE==VRBL) then
	Call AMBI_DENS(m,kbm1,s_in(1:m,1:kbm1),t_in(1:m,1:kbm1),rho_water(1:m,1:kbm1))
  end if

  !------------------------------------------------------
  !Set up Model Forcing  -> time step and bottom shear 
  !Convert taux/tauy to node-based values
  !Convert Model tau/rho to tau [N/m^2]
  !------------------------------------------------------
  DTsed = DTin
  T_model = Tin
  do i=1,m
	taub(i) = taub_in(i)*rho_water(i,kbm1)
  end do

  !------------------------------------------------------
  !Vertical Plastic Dynamics --> Plastic 'Model'
  !------------------------------------------------------

  !Calculate Active Layer Thickness --> bottom(:,iactv)
  call calc_active_layer
    
  !calculated depositional flux from settling (settle flux function is called here)
  call calc_deposition

  !calculate erosive flux (resuspension function is called here)
  call calc_erosion

  !adjust bottom properties
  call update_bottom

  !------------------------------------------------------
  !Horizontal Advection of Sediment Concentration 
  !------------------------------------------------------
  if(.not. oned_model)then
     d_cdis = 1
     if(numqbc > 0) d_cdis = numqbc
     d_cflx = 1
     if(Nobc > 0) d_cflx = Nobc
      
     do i=1,nplast
      
      !set discharge sediment concentrations
      if(numqbc > 0 .and. PLAST_source)then
        plast(i)%cdis(1:numqbc) = seddis(1:numqbc,i)
      endif
  
      !3D advection
      call adv_scal(plast(i)%conc, &
                      plast(i)%cnew, &
                      d_cdis,      &
                      plast(i)%cdis(1:numqbc), &
                      d_cflx,      &
                      plast(i)%cflx, &
                      DTsed,PLAST_source)

	  !This concentration updating function is for sediment, no need to change here.
      call fct_sed(plast(i)%conc,plast(i)%cnew) 
  
    end do
  else
    do i=1,nplast
      plast(i)%cnew = plast(i)%conc
    end do
  endif 

  !------------------------------------------------------
  !Vertical Diffusion of Sediment Concentrations
  !------------------------------------------------------
  do i=1,nplast
    call vdif_scal(plast(i)%cnew,DTsed) 
  end do

  !------------------------------------------------------
  !Exchange Concentration on Processor Boundaries
  !------------------------------------------------------
# if defined (MULTIPROCESSOR)
  if(par)then
    do i=1,nplast
      call aexchange(nc,myid,nprocs,plast(i)%cnew)
    end do
  endif
# endif

  !------------------------------------------------------
  ! Set Open Boundary Conditions (nudging to background)
  !------------------------------------------------------
  if(.not. oned_model)then
    if(PLAST_ramp == 0) PLAST_ramp = 1
    temp = min(float(PLAST_its/PLAST_ramp), 1.0)*PLAST_alpha
    if(.not.plast_nudge) temp = 0.0
    do i=1,nplast 
       call bcond_scal_OBC(plast(i)%conc, &
                           plast(i)%cnew, &
                           plast(i)%cflx, &
                           plast(i)%cobc, &
                           DTsed,       &
                           temp) 
    end do
  endif

  !------------------------------------------------------
  ! Point Source Boundary Conditions 
  !------------------------------------------------------
  if(PLAST_source .and. .not. oned_model)then
    do i=1,nplast
      call bcond_scal_PTsource(plast(i)%conc, &
                               plast(i)%cnew, &
                               plast(i)%cdis) 
    end do
  endif

  !------------------------------------------------------
  ! Update Concentration Variables to next time level
  !------------------------------------------------------
  do i=1,nplast
    plast(i)%conc = plast(i)%cnew
  end do

  !------------------------------------------------------
  ! Limit to Non-Negative (modify to Venkat Limiter) 
  !------------------------------------------------------
  do i=1,nplast
    plast(i)%conc = max(plast(i)%cnew,0.0)
  end do
  
  !------------------------------------------------------
  !Exchange Concentration/tau/rho on Processor Boundaries
  !------------------------------------------------------
# if defined (MULTIPROCESSOR)
  if(par)then
    do i=1,nplast
      call aexchange(nc,myid,nprocs,plast(i)%cnew)
    end do
    call aexchange(nc,myid,nprocs,taub)
	call aexchange(nc,myid,nprocs,rho_water)
  endif
# endif

  !------------------------------------------------------
  ! Report 
  !------------------------------------------------------
  if(mod(PLAST_its,n_report)==0) call PLAST_report

!  call layer_report(1)
  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Advance_Plast "

 End Subroutine Advance_Plast

!=======================================================================
! Allocate Sediment Variables                      
!=======================================================================
Subroutine Alloc_Plast_Vars
  
	use lims,     only: nt,mt,kb,kbm1
  
	implicit none
  
	integer i,tmp1,tmp2

	if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Alloc_PLAST_Vars "

	!ensure arrays have nonzero dimension
	tmp1 = max(Nobc+1,1)
	tmp2 = max(Nobc  ,1)

	!allocate suspended sediment arrays

	do i=1,nplast
		plast(i)%arraysize = mt+1
		allocate(plast(i)%conc(0:mt,kb  ))      ; plast(i)%conc = 0.0
		allocate(plast(i)%cnew(0:mt,kb  ))      ; plast(i)%cnew = 0.0
		allocate(plast(i)%mass(0:mt     ))      ; plast(i)%mass = 0.0
		allocate(plast(i)%frac(0:mt     ))      ; plast(i)%frac = 0.0
		allocate(plast(i)%eflx(0:mt     ))      ; plast(i)%eflx = 0.0
		allocate(plast(i)%dflx(0:mt     ))      ; plast(i)%dflx = 0.0
		allocate(plast(i)%cdis(500      ))      ; plast(i)%cdis = 0.0
		allocate(plast(i)%cflx(tmp1,kbm1))      ; plast(i)%cflx = 0.0
		allocate(plast(i)%cobc(tmp2     ))      ; plast(i)%cobc = 0.0
		allocate(plast(i)%wset(0:mt,kbm1))      ; plast(i)%wset = 0.0
		allocate(plast(i)%zactv(0:mt    ))      ; plast(i)%zactv = 0.0
		allocate(plast(i)%erate(0:mt    ))      ; plast(i)%erate = 0.0
		allocate(plast(i)%tau_ce(0:mt   ))      ; plast(i)%tau_ce = 0.0
	end do

	!allocate bottom shear stress array 
	allocate(taub(0:mt)) ; taub = 0.0

	!allocate bottom array
	allocate(bottom(0:mt , n_bot_vars)) ; bottom = 0.0
	
	!allocate the water density array
	allocate(rho_water(0:mt,kbm1)); rho_water = 0.0

	if(dbg_set(dbg_sbr)) write(ipt,*) "End: Alloc_PLAST_Vars "

 End Subroutine Alloc_PLAST_Vars

!==========================================================================
! Report Sediment Setup To Screen 
!==========================================================================
Subroutine Report_PLAST_Setup
  
	use control, only: msr 
  
	implicit none
  
	integer :: i

	if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Report_PLAST_Setup "

	! echo sediment model parameters to screen
	if(msr)then

		write(*,*)
		write(*,*)'------------------ Sediment Model Setup --------------------'
		if(oned_model)then
			write(*,*)'!  model                 : oneD   '
		endif
		  
		write(*,*)'!  nplast                  :',nplast    
		write(*,*)'!  n_report              :',n_report
		
		if(PLAST_start > 0)then
			write(*,*)'!  PLAST_start             :',PLAST_start
		else
			write(*,*)'!  PLAST_start             : no delay'  
		endif
		
		if(PLAST_hot_start)then
			write(*,*)'!  PLAST_hot_start         :',PLAST_hot_start
		endif
    
		if(plast_nudge)then
			write(*,*)'!  open boundary nudging :  active'
			write(*,*)'!  nudging relax factor  :',PLAST_alpha 
			write(*,*)'!  # nudging ramp its    :',PLAST_ramp  
		else
			write(*,*)'!  open boundary nudging :  not active'
		endif
		
		if(PLAST_source)then
			write(*,*)'!  point source forcing  :  active'
		else
			write(*,*)'!  point source forcing  :  not active'
		endif
  
		if(inf_bed)then
			write(*,*)'!  plastic supply       :  infinite'
		else
			write(*,*)'!  plastic supply       :  finite'
		endif
		
		write(*,*)'!  Settling CFL          :',settle_cfl
		write(*,*)'!'
				
		write(*,*)'|    class       Aaxis   Baxis   Caxis    Prho   |'
	
		do i=1,nplast
			write(*,'(1X,A1,1X,A12,1X,3(F7.4,1X),1X,F8.2,1X,A1)')'|', &
			trim(plast(i)%pname), plast(i)%aaxis,plast(i)%baxis,plast(i)%caxis,plast(i)%Prho,'|'
		end do
	
	end if

  ! echo sediment model initial conditions to screen
  call PLAST_report  

  write(*,*)'------------------------------------------------------------'
  write(*,*)

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Report_PLAST_Setup "

  End Subroutine Report_PLAST_Setup

!==========================================================================
! Monitor Sediment Model: Compute and Write Statistics to Screen
!==========================================================================
Subroutine PLAST_Report  
  
  use control, only: msr ,par, mpi_fvcom_group
  use lims   , only: m,mt,kbm1, msrid
  
  implicit none
  
  real(dp) :: tau_min,tau_max,tmp1(3),tmp2(3)
  integer :: i,dim1,dim2,ierr

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: PLAST_Report "

  !set up limits
  dim1 = mt 
  dim2 = kbm1
  
  !calculate max/min/rms concentrations
  do i=1,nplast
    plast(i)%cmax = maxval(plast(i)%conc(1:dim1,1:dim2))
    plast(i)%cmin = minval(plast(i)%conc(1:dim1,1:dim2))
    plast(i)%crms = sum(dble(abs(plast(i)%conc(1:dim1,1:dim2))))
    plast(i)%crms = plast(i)%crms/(dim1*dim2)
  end do

  tau_max   = maxval(taub(1:m))
  tau_min   = minval(taub(1:m))
 
# if defined (MULTIPROCESSOR)
  IF(PAR)THEN
  do i=1,nplast
    tmp1(1:3) = (/plast(i)%cmax,plast(i)%cmin,plast(i)%crms/)
    CALL MPI_REDUCE(tmp1,tmp2,3,MPI_DP,MPI_MAX,MSRID-1,MPI_FVCOM_GROUP,IERR)
    if(msr)then
      plast(i)%cmax = tmp2(1)
      plast(i)%cmin = tmp2(2)
      plast(i)%crms = tmp2(3)
    endif
  end do

  tmp1(1:2) = (/tau_min,tau_max/)
  CALL MPI_REDUCE(tmp1(1),tmp2(1),1,MPI_DP,MPI_MIN,MSRID-1,MPI_FVCOM_GROUP,IERR)
  CALL MPI_REDUCE(tmp1(2),tmp2(2),1,MPI_DP,MPI_MAX,MSRID-1,MPI_FVCOM_GROUP,IERR)
  tau_min   = tmp2(1) 
  tau_max   = tmp2(2) 

  ENDIF
# endif

  ! write statistics to screen
  if(.not.msr)return
  write(*,*  )'====================Sediment Model Stats======================'
  write(*,*  )'!  quantity              :     avg(g/L)      max(g/L)      '
  do i=1,nplast
    write(*,100)'!',trim(plast(i)%pname),'    :',plast(i)%crms,plast(i)%cmax
  end do
  write(*,*  )'!  max bottom stress     :     ',tau_max
  write(*,*  )'!  min bottom stress     :     ',tau_min 

  100 format(1x,a1,a20,a5,3f12.6)
  101 format(1x,a26,f12.6)

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: PLAST_Report "

  End Subroutine PLAST_Report  

!==========================================================================
! Calculate Erosion Rate (kg/m^2) 
! Calculate erosion rate using the method of Ariathurai and Arulanandan (1978)
!    -Erosion rates of cohesive soils. J. Hydraulics Division, ASCE, 104(2),279-282.-
!    surface_mass_flux[i] = erate[i]*(1-porosity)*bfrac[i]*(tau_w/tau_crit[i] - 1)
!      where
!      iplast            = index of sediment type
!      bfrac[i]        = fraction of bed composed by sediment i
!      erate[i]        = bed erodibility constant
!      porosity        = volume of voids/total volume in the top layer
!      tau_w           = shear stress on the bed
!      tau_ce          = critical shear stress of sediment
!==========================================================================

Subroutine Calc_Erosion
  
  use all_vars
  
  implicit none
  
  integer  :: iplast ,i
  real(sp) :: bottom_por,bottom_frac,tau_w,dep,active,erosion_orig
  real(sp) :: tau_ce(m)
  real(sp) :: Erate(m)
  
  real(sp) :: dx,dy 

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Calc_Erosion "

  !Compute Erosive Flux (kg/m^2), limited by available material in active layer 
  do iplast=1,nplast
	
	if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Critical shear stress computation"
		!Compute the critical stress for erosion
		Call Calc_critial_shear_WS(m,rho_water(1:m,kbm1),bottom(1:m,isd50), &
			plast(iplast)%Prho,plast(iplast)%dequ,tau_ce)
		do i=1,m
			plast(iplast)%tau_ce(i)=tau_ce(i)
		end do
	if(dbg_set(dbg_sbr)) write(ipt,*) "End: Critical shear stress computation"
	
	if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Erosion coefficient"
		!Compute the sediment-like erosion rate based on mud-portion and grain size
		Call Calc_Erosion_Rate(m,tau_ce,bottom(1:m,imudp),bottom(1:m,idens),Erate)
		do i=1,m
			plast(iplast)%erate(i)=Erate(i)
		end do
	if(dbg_set(dbg_sbr)) write(ipt,*) "End: Erosion coefficient"
  
	do i=1,m
# if defined (WET_DRY)
	  if(iswetn(i)/=1)cycle  ! no calculation when dry node
# endif
	  
		bottom_por  = bottom(i,iporo)
		bottom_frac = plast(iplast)%frac(i)
		tau_w    = taub(i)
		!dep      = plast(iplast)%dflx(i)
		!active   = (1.0-bed_por)*Bed_Frac*(tau_w/tau_ce-1.0)*bottom(i,iactv) + dep
		erosion_orig  = MAX(DTsed*Erate(i)*(1.0-bottom_por)*bottom_frac*(tau_w/tau_ce(i)-1.0),0.0)
		plast(iplast)%eflx(i) = erosion_orig 

    end do
  end do

  !Update Concentration in Bottom of Water Column
  do iplast=1,nplast
    do i=1,m
# if defined (WET_DRY)
      if(iswetn(i)==1)then
# endif
        plast(iplast)%conc(i,kbm1) = plast(iplast)%conc(i,kbm1) + plast(iplast)%eflx(i)/(d(i)*dz(i,kbm1))
# if defined (WET_DRY)
      else
        plast(iplast)%conc(i,kbm1) = 0.0_sp
      end if
# endif
    end do
  end do
  
  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Calc_Erosion "

End Subroutine Calc_Erosion

!==========================================================================
! Calculate critical shear stress 
!==========================================================================

Subroutine Calc_critial_shear_WS(m,rho_wb,dsed,rho_p,dequ,tau_ce)

	implicit none
	
	!Input variables
	integer, intent(in)		:: m
	real(sp), intent(inout)   :: tau_ce(m)
	real(sp), intent(in) 	  :: rho_wb(m) !Bottom water density
	real(sp), intent(in)      :: dsed(m)
	real(sp), intent(in)	  :: rho_p
	real(sp), intent(in)      :: dequ
  
	!Local variables
	integer  :: i
	real :: R			   !Relative density
	real :: Dstar		   !Dimensionless diameter
	real :: theta_sc	   !Shields parameter for sediment
	real :: theta_pc       !Shields parameter for plastic
	
	real, parameter :: nu=1E-6	       !Water viscosity at 20deg
	real, parameter :: g=9.81		   !gravitational acceleration
	real, parameter :: HE_coef=-0.503   !Hide-Expose coefficient
	real, parameter :: Pre_coef=0.5588   !Pre-coefficient

	do i=1,m
		
		!Calculate relative density
		R = abs((rho_p-rho_wb(i))/rho_wb(i))
	
		!Calculate Dimensionless diameter
		Dstar=(R*(g/nu**2))**0.3333*dsed(i)
	
		!Calculate Shields parameter for sediment using Van Rijn 1993
		if (Dstar>=1 .and. Dstar<4) then
			theta_sc=0.24*Dstar**(-1)
		elseif (Dstar>=4 .and. Dstar<10) then
			theta_sc=0.14*Dstar**(-0.64)
		elseif (Dstar>=10 .and. Dstar<20) then
			theta_sc=0.04*Dstar**(-0.1)
		elseif (Dstar>=20 .and. Dstar<150) then
			theta_sc=0.013*Dstar**0.39
		elseif (Dstar>=150) then
			theta_sc=0.055
		end if
	
		!Calculate Shields parameter for plastic
		theta_pc=Pre_coef*(dequ/dsed(i))**HE_coef*theta_sc
	
		!Calculate tau_ce
		tau_ce(i)=g*(rho_p-rho_wb(i))*theta_pc*dequ
		
	end do
	
End Subroutine Calc_critial_shear_WS

!==========================================================================
! Calculate the erosion rate base on Wu et al. 
!==========================================================================

Subroutine Calc_Erosion_Rate(m,tau_ce,mudp,dens,Erate)
	
	implicit none
	
	!Input/Output variables
	integer, intent(in)		:: m
	real(sp), intent(inout)   :: Erate(m)
	real(sp), intent(in) 	  :: tau_ce(m)
	real(sp), intent(in)      :: mudp(m)
	real(sp), intent(in)	  :: dens(m)

	!Local variables
	integer :: i
	real, parameter :: ks0=7.6E-5 !pure-sand erosion coefficient
	
	do i=1,m
		Erate(i)=dens(i)*tau_ce(i)*ks0*(1.035*exp(-19.87*mudp(i)**0.86)+0.0981*mudp(i)**0.15)
	end do

End Subroutine Calc_Erosion_Rate

!==========================================================================
! Calculate Depositional Flux and Update Concentration from Settling
!==========================================================================

Subroutine Calc_Deposition
  
  use all_vars
  
  implicit none
  
  integer  :: iplast ,i,mcyc,ncyc,k
  real(sp) :: wset(kbm1),c(kbm1),dx(kbm1),flux,DTmax,eps,DTdep,da,db,dc,rhop

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Calc_Deposition "

  eps = epsilon(eps)

  !Loop over Sediment Class       
  do iplast=1,nplast
    do i=1,m
  
      ! no calculation when dry nodes
# if defined (WET_DRY)
      if (iswetn(i)/=1) cycle  
# endif

      !initialize flux
      flux = 0.0

      !setup 1D concentration and grid arrays
      c(1:kbm1) = plast(iplast)%conc(i,1:kbm1)
      dx(1:kbm1) = dz(i,1:kbm1)*d(i)
      
	  rhop = plast(iplast)%Prho
	  da = plast(iplast)%aaxis
	  db = plast(iplast)%baxis
	  dc = plast(iplast)%caxis
	  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Calc_Wset_Plastic"
		call Calc_Wset_Plastic(kbm1,rho_water(i,1:kbm1),rhop,da,db,dc,wset)
	  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Calc_Wset_Plastic"

	  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Settling advection"
		!set up cycles (use max(CFL) == 1)
		DTmax = Settle_CFL*minval(dx(1:kbm1)/wset(1:kbm1))
		mcyc = int(DTsed/DTmax + 1. - eps)
		DTdep = DTsed/float(mcyc)

		!call flux-limited settling equation
		do ncyc = 1,mcyc
			call settle_flux(kbm1,c,dx,wset,flux,DTdep)
		end do

		!store solution in 3D array
		plast(iplast)%conc(i,1:kbm1) = c(1:kbm1)
	  
		!store wset in 3D array
		plast(iplast)%wset(i,1:kbm1) = wset(1:kbm1)

		!store depositional flux
		plast(iplast)%dflx(i) = flux
	  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Settling advection"	  

      !shutdown erosion if node is inactive 
      !plast(iplast)%dflx(i) = plast(iplast)%dflx(i)
    end do
  end do

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Calc_Deposition "

  End Subroutine Calc_Deposition

!==========================================================================
! Calculate settling velocity of plastic using the method of
!    Mehta, 2021
!==========================================================================
 
Subroutine Calc_Wset_Plastic(n,rhow,rhop,da,db,dc,wset) 
	
	implicit none
	
	!Input variables
	integer , intent(in)  :: n
	real(sp), intent(inout)   :: wset(n)
	real(sp), intent(in) 	  :: rhow(n)
	real(sp), intent(in)      :: rhop
	real(sp), intent(in)	  :: da
	real(sp), intent(in)	  :: db
	real(sp), intent(in)      :: dc
  
	!Local variables
	integer  :: i
	real :: R			   !Relative density
	real :: csf		       !Corey-Shape-Factor csf
	real :: E			   !Jenkes Coefficient E
	real :: nc,C1,C2       !Shape related coefficient (nc is the same as n in the orignal paper)
	real :: Dg,Ds		   !Reference particle size, dimensionless reference particle size 
	real :: Ws		       !Dimensionless settling velocity
	real, parameter :: nu=1E-6	       !Water viscosity at 20deg
	real, parameter :: g=9.81		   !gravitational acceleration
	real, parameter :: coef_a=18      !coefficient for the shape factor C1
	real, parameter :: coef_b=-0.38   !coefficient for the shape factor C1
	real, parameter :: coef_c=0.3708  !coefficient for the shape factor C2
	real, parameter :: coef_d=-0.1602 !coefficient for the shape factor C2
	real, parameter :: coef_e=0.4942  !coefficient for the shape factor nc
	real, parameter :: coef_f=-0.059  !coefficient for the shape factor nc

	csf=dc/sqrt(da*db)
	E=da*((da**2+db**2+dc**2)/3)**(-0.5)
	Dg=da*(csf)**0.34*(db/da)**(0.5)
  
	C1=coef_a*E**coef_b
	C2=coef_c*csf**coef_d
	nc=coef_e*csf**coef_f
  
	do i=1,n
		R=(rhop-rhow(i))/rhow(i)
		Ds=Dg*(g*abs(R)*nu**(-2.0))**(0.3333)
		Ws=Ds**2/(C1+(0.75*C2*Ds**3)**nc)
		wset(i)=sign(Ws*(g*abs(R)*nu)**(0.3333),R)
	enddo
  
 End Subroutine Calc_Wset_Plastic
  
!==========================================================================
! Calculate Settling Flux and Update Sediment Concentration in Water Column 
!   use the SLIP limiter
!   Jameson, A., Analysis and Design of Numerical Schemes for Gas Dynamics I
!      Artificial Diffusion, Upwind Biasing, Limiters, and their effects
!      on Accuracy and Multigrid Convergence, International Journal of 
!      Computational Fluid Dynamics, Vol 4, 171-218, 1995.
!
!      -second order accurate in smooth regions
!      -guaranteed local extremum diminishing
!      
!   conv:  convective flux
!   diss:  dissipative flux
!   
!==========================================================================
Subroutine Settle_Flux(n,c,dx,wset,flux,deltat) 
  implicit none
  integer , intent(in   ) :: n
	  real(sp), intent(inout) ::  c(n)
	  real(sp), intent(in   ) :: dx(n)
	  real(sp), intent(in   ) :: wset(n)
	  real(sp), intent(inout) :: flux 
	  real(sp), intent(in   ) :: deltat 
	  real(sp) :: conv(n+1),diss(n+1)
	  real(sp) :: cin(-1:n+2)
	  real(sp) :: win(-1:n+2)
	  real(sp) :: dis4,wvel
  integer  :: i

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Settle_Flux "

  !transfer to working array
  cin(1:n) = c(1:n)
  win(1:n) = wset(1:n) 

  !surface bcs (no flux)
  cin(0)  =  -cin(1) 
  cin(-1) =  -cin(1)
  win(0)  =  -wset(1)
  win(-1) =  -wset(1)

  !bottom bcs (extrapolate)
  cin(n+1) = cin(n) 
  cin(n+2) = cin(n)
  win(n+1) = win(n)
  win(n+2) = win(n)

  !flux computation
  do i=1,n+1
    wvel    = .5*(win(i)+win(i-1))  !settle velocity at interface
    dis4    = wvel/2.
    conv(i) = wvel*(cin(i)+cin(i-1))/2. 
    diss(i) = dis4*(cin(i)-cin(i-1)-lim(cin(i+1)-cin(i),cin(i-1)-cin(i-2))) 
  end do

  !zero out surface flux
  conv(1) = 0. ; diss(1) = 0.

  !update
  do i=1,n
    c(i) = cin(i) + (deltat/dx(i))*(-conv(i+1)+conv(i) + diss(i+1)-diss(i)) 
  end do

  !set bottom flux (> 0 = deposition)
  flux = flux + deltat*(conv(n+1)-diss(n+1))

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Settle_Flux "

 End Subroutine Settle_Flux

!==========================================================================
! Calculate LED Limiter L(u,v)  
!==========================================================================
Function Lim(a,b)
  real(sp) lim,a,b
  real(sp) q,R
  real(sp) eps
  eps = epsilon(eps)
  
 ! exponent
 ! qlim = 0. !1st order
 ! qlim = 1. !minmod
 ! qlim = 2. !van leer

  R = abs(   (a-b)/(abs(a)+abs(b)+eps) )**qlim
  lim = .5*(1-R)*(a+b)

  End Function Lim

!==========================================================================
! Update Surface Bed Properties Due to SusLoad Flux 
! Change:
!   plast(:)%mass(:,1)   to reflect increased/decrease mass (kg)
!   bed(:,1,itchk)     to reflect changed bed thickness properties
!==========================================================================
Subroutine Update_Bottom
  use lims, only: m
  implicit none
  real(sp) :: dz,accum,flux,eps
  integer  :: i,iplast
  real(sp) :: afrac(m)

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Update_Bottom"

  eps = epsilon(eps)
  
  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Mass depositing"
  do i=1,m
    do iplast=1,nplast
      flux = (plast(iplast)%eflx(i)-plast(iplast)%dflx(i))

      !change top layer bed mass of each plastic type 
      accum  = plast(iplast)%mass(i) - flux 
      plast(iplast)%mass(i) = max(accum,0.0) 

    end do
  end do
  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Mass depositing"
  
  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Fraction computation"
  do iplast=1,Nplast
	 !Compute fraction based on active layer thickness (move to update_bottom)
	 Call Calc_plastic_fraction(m,plast(iplast)%mass,plast(iplast)%zactv, &
			bottom(1:m,iporo),bottom(1:m,idens),afrac)
	 do i=1,m
		plast(iplast)%frac(i)=afrac(i)
	 end do 
   end do
  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Fraction computation"

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Update_Bottom"

End Subroutine Update_Bottom
 
Subroutine Calc_Plastic_Fraction(m,mass,zactv,poro,dens,afrac)

	implicit none
	
	!Input/Output variable
	integer, intent(in)  :: m
    real(sp), intent(inout) :: afrac(m)
	real(sp), intent(in)    :: mass(m) 
	real(sp), intent(in)	:: zactv(m)
	real(sp), intent(in)	:: poro(m)
	real(sp), intent(in)    :: dens(m)
	
	!Local variables
	integer :: i
	
	do i=1,m
		afrac(i)=mass(i)/(mass(i)+dens(i)*(1-poro(i))*zactv(i))
	end do

End Subroutine Calc_Plastic_Fraction
!==========================================================================
! Update Active Layer Thickness 
!   Use method of Harris and Wiberg (1997):
!   -Approaches to quantifying long-term continental shelf sediment transport
!   -with an example from the northern California STRESS mid-shelf site.
!   [Continental Shelf REsearch, 17, 1389-1418]
!
!   z_a = max[k_1*(tau_bot - tau_ce)*rho_0 , 0] + k_2 * D_50
!       where
!       k_1 = .007 (empirical constant)
!       k_2 = 6.0  (empirical constant)
!       tau_bot = shear stress on the bed
!       tau_ce  = critical shear stress for erosion of the bed (averaged over
!                 sediment classes)
!       D_50    = median grain diameter of the surface sediment
!       note: our shear stresses are nondimensionalized by rho, so
!             rho does not appear in actual calculation (below)
!==========================================================================
Subroutine Calc_Active_Layer 
  use lims,     only: m
  implicit none
  integer  :: i,iplast

  !calculate active layer thickness 
  do iplast=1,Nplast
	do i=1,m
		plast(iplast)%zactv(i) = max(0.,.007*(taub(i)-plast(iplast)%tau_ce(i))) &
			+6.0*bottom(i,isd50)
	end do
  end do

End Subroutine Calc_Active_Layer
  
!==========================================================================
!  Setup open boundary nudging for sediment obc
!   - check for existence of sed nudging file
!   - read data (concentration of each sed type at each obc)
!==========================================================================
Subroutine Setup_PLAST_OBC

  Use Control,  only : casename,input_dir,msr,serial,par
!  Use Mod_OBCs, only : iobcn,i_obc_gl,iobcn_gl
  Use Mod_Utils, only: pstop
# if defined (MULTIPROCESSOR)
  Use Mod_Par,  only : nlid
# endif

  integer :: i,j,tmp,ncnt,i1
  integer,  allocatable :: vtmp(:),atemp(:,:)
  real(sp), allocatable :: PLAST_obc_gl(:,:)
  integer, parameter    :: plastobc = 81
  character(len=120)    :: fname
  logical :: fexist
  

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Setup_PLAST_OBC "

  !------------------------------------------------------
  ! make sure sediment nudging data file exists 
  !------------------------------------------------------
  fname = "./"//trim(input_dir)//"/"//trim(casename)//"_plast_nudge.dat"

  inquire(file=trim(fname),exist=fexist)
  if(msr .and. .not.fexist)then
    write(*,*)'sediment nudging data: ',fname,' does not exist'
    write(*,*)'halting.....'
    call pstop
  end if

  open(unit=plastobc,file=fname,form='formatted')
  !------------------------------------------------------
  ! allocate global data
  !------------------------------------------------------
  allocate(PLAST_obc_gl(Nobc_gl,nplast)) ; PLAST_obc_gl = 0.

  !header
  read(plastobc,*)
  read(plastobc,*)
  read(plastobc,*)

  !# obc nodes
  read(plastobc,*) tmp 

  if(tmp /= Nobc_gl .and. msr)then
    write(*,*)'number of open boundary nodes for sediment nudging: ',tmp
    write(*,*)'this is not consistent with number of open boundary nodes in'
    write(*,*)'the model: ',Nobc_gl
    write(*,*)'halting'
    call pstop
  endif

  allocate(vtmp(Nobc_gl))
  do i=1,Nobc_gl  
    read(plastobc,*) i1,vtmp(i),(PLAST_obc_gl(i,j),j=1,nplast)
  end do
  close(plastobc)
    
  !---------------------------------------------------------
  ! shift plastic nudging data to local open boundary nodes
  !---------------------------------------------------------

  if(serial)then
    do i=1,Nobc_gl
      do j=1,nplast
        plast(j)%cobc(i) = PLAST_obc_gl(i,j)
      end do
    end do
  endif

# if defined (MULTIPROCESSOR)
  if(par)then
    allocate(atemp(Nobc_gl,nplast))
    ncnt = 0
    !!set up local open boundary nodes
    do i=1,Nobc_gl
      i1 = nlid( vtmp(i) )
      if(i1 /= 0)then
        ncnt = ncnt + 1
        do j=1,nplast
          atemp(ncnt,j) = PLAST_obc_gl(i,j)
        end do
      end if
    end do

    !transfer plastic nudging data global --> local
    if(ncnt > 0)then
      do j=1,nplast
        plast(j)%cobc(1:ncnt) = atemp(1:ncnt,j)
      end do
    end if

    deallocate(atemp) 
  end if
# endif

  if(msr)write(*,*)'obc setup for plastic model: complete'

  !---------------------------------------------------------
  ! cleanup 
  !---------------------------------------------------------
  deallocate(PLAST_obc_gl)

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Setup_PLAST_OBC "

  End Subroutine Setup_PLAST_OBC
  
  SUBROUTINE AMBI_DENS(mb,kb,MYS,MYT,MYRHO)               

    !==============================================================================|
    IMPLICIT NONE
	INTEGER :: mb,kb
	REAL(SP), INTENT(IN), DIMENSION(:,:) :: MYS,MYT
    REAL(SP), INTENT(INOUT), DIMENSION(:,:):: MYRHO
	
	INTEGER :: I,K,ub1,lb1,ub2,lb2

    !==============================================================================|

    ! SET DIMS USING MY S
    lb1 = lbound(MYS,1)
    ub1 = ubound(MYS,1)

    lb2 = lbound(MYS,2)
    ub2 = ubound(MYS,2)

    ! CHECK MY T
    if(lb1 /= lbound(myt,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")
    if(ub1 /= ubound(myt,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")

    if(lb2 /= lbound(myt,2)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")
    if(ub2 /= ubound(myt,2)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")

    ! CHECK MY RHO
    if(lb1 /= lbound(myrho,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")
    if(ub1 /= ubound(myrho,1)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")

    if(lb2 /= lbound(myrho,2)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")
    if(ub2 /= ubound(myrho,2)) CALL FATAL_ERROR("EQS_OF_STATE:: Dimension mismatch!")
    
    !
    !  CALCULATE DENSITY FROM EQUATION OF STATE
    !
	
	Do I=1,mb
	Do K=1,kb
	
    MYRHO(I,K) =MYS(I,K)*MYS(I,K)*MYS(I,K)*6.76786136E-6_SP &
         & - MYS(I,K)*MYS(I,K)*4.8249614E-4_SP &
         & + MYS(I,K)*8.14876577E-1_SP &
         & - 0.22584586E0_SP

    MYRHO(I,K) = MYRHO(I,K) * &
         & ( MYT(I,K)*MYT(I,K)*MYT(I,K)*1.667E-8_SP &
         & - MYT(I,K)*MYT(I,K)*8.164E-7_SP &
         & + MYT(I,K)*1.803E-5_SP &
         & )
    
    MYRHO(I,K) = MYRHO(I,K)+1.0_SP &
         & - MYT(I,K)*MYT(I,K)*MYT(I,K)*1.0843E-6_SP &
         & + MYT(I,K)*MYT(I,K)*9.8185E-5_SP &
         & - MYT(I,K)*4.786E-3_SP
    
    MYRHO(I,K) = MYRHO(I,K)* &
         & ( MYS(I,K)*MYS(I,K)*MYS(I,K)*6.76786136E-6_SP &
         & - MYS(I,K)*MYS(I,K)*4.8249614E-4_SP &
         & + MYS(I,K)*8.14876577E-1_SP &
         & +3.895414E-2_SP &
         & )
    
    MYRHO(I,K) = MYRHO(I,K) &
         & - (MYT(I,K)-3.98_SP) *(MYT(I,K)-3.98_SP) * (MYT(I,K)+283.0_SP) / (503.57_SP*(MYT(I,K)+67.26_SP))
    
		!
		!  CALCULATE RHO
		!
		MYRHO(I,K) =  MYRHO(I,K)+1000.
		
	End Do
	End Do

  END SUBROUTINE AMBI_DENS
  
# endif
   
End Module Mod_Sed 


