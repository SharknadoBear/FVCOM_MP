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
! Contact:      Yicheng Huang, Pacific Northwest National Laboratory
!
! Micro-Plastic transport model (refer to as FVCOM-MP) developed based on 
! the structure of Community Sediment Transport Model (CSTM) as implemented 
! in ROMS by J. Warner (USGS), with all the transport and fate dynamics adapted
! for the plastics.
!
! Comments:     Micro-Plastics Transport and Fate Module FVCOM-MP
!
! FVCOM dependency
!
!   init_sed.F:  - user defined plastic model initial conditions
!   mod_ncdio.F: - netcdf output includes concentration/bottom/bed fields
!	mod_bbl.F	 - bottom boundary layer module (now it is not related)
!	mod_main.F   - miscellaneous global variables and major namelist structure
!	mod_input.F  - input module (namelist, inital field from ncfiles)
!	mod_ncdio.F  - output module
!	mod_startup.F- set-up all the variables in memory before computation
!   fvcom.F:     - main calls plastic/plastic setup 
!   internal_step.F:  - calls plastic advance
!
! History
!   Dec 30,2022: merge the whole module from the original mod_sed.F, delete
!               all the redundent part, keep the essential features for a beta 
!				version that could be used for testing the new plastic settling
!				module. Still using the plastic module and most of the 
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
!				the module is almost isolated from the original plastic module. 
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
!				density. The plastic on the bottom is resuspended as plastic, 
!				but with the amount related to the sediment bed feature (active-layer-
!				thickness, mud-content, sediment density and porosity). 
!   Apr 3, 2023: Now the module is compiled as a single module, it can be run with 
!				sediment. In this version, the major commission is the toxic-vector 
!				functional. When the toxic-vector module is turned on the bonded toxic 
!				and the free toxic concentration can be output. The toxic releasing 
!				process is modelled by first-order consecutive chemical reaction process.
!				Be careful, there are three substances in the plastic type now, 
!				i.e. sediment, Micro-Plastic, toxic. The sediment is using (kg/m^3 or
!				g/L) as before. The MP is using (g/m^3 or mg/L). The toxic is using
!				(mg/m^3 or microgram/L).
!   Jun 2nd,2024: Concentration units are corrected to #/m^3; the flux units are
!				corrected to #/m^2. Erosion rate is updated to the latest version,
!				separating between 1) tau_w < tau_ce; 2) tau_ce < tau_w < tau_ce_sed
!				3) tau_ce > tau_ce_sed
!   May 30, 2025: Now toxic-vector function is refined specifically on several aspect:
!				1) New toxic releasing functionality is included: now there are a four-modes 
!				 shifting mechanism is added into the namelist, a) FOK; b) SI-FOK; c) SIR-FOK
!														a) simple first order decaying, with constant decaying rate
!               										b) shape-impacted direct-leaching 
!														c) shape-impacted reactive-leaching 
!				2) River sources for sediment/plastic/toxins are now correctly configured
!				3) Merge with latest fixing on the transport equation including five modes on bottom treatment
!				4) Now all the input/output of either MP or toxin data are in the SI unit
!   October 1, 2025 [Work on-going...]: Some preliminary flocs function is added-in:
!				1) Communication between sediment and plastic module, passing sediment concentration in
!				2) Commonication between hydrodynamic part and plastic module, passing the turbulence kinetic energy in
!				3) 2-tracer model, when consider flocs, the MP concentration is computed through 2 parts
!
!==========================================================================/
Module Mod_Plast  

#if defined (PLASTIC)
	Use Mod_Par
	Use Mod_Prec 
	Use Mod_Types
	Use Mod_wd
	Use Control, only : TOXBDIS,TOXFDIS,PLASTDIS    !This is related to river need to be changed for plastic case.
	Use all_vars,only : CNSTNT,UNIFORM

	implicit none 

	!------------------------------------------------------------------------------------------
	!Plastic Type                         
	!
	! pname        => plastic name
	! MP_shape     => model shape name 'sphere'/'cylinder'/'cuboid'/'ellipsoid'/'arbitrary'
	! tname_b	   => bonded toxic name
    ! tname_f 	   => free toxic name 	
	! cmax         => Report variable: the maximum value of plastic concentration at each reporting step			[kg m^-3       ]
	! cmin		   => Report variable: the minimum value of plastic concentration at each reporting step			[kg m^-3       ]
	! crms		   => Report variable: the rms value of plastic concentration at each reporting step				[kg m^-3       ]
	! tcmax_f	   => Report variable: the maximum value of free MP toxic concentration at each reporting step		[kg m^-3       ]
	! tcmin_f 	   => Report variable: the minimum value of free MP toxic concentration at each reporting step		[kg m^-3       ]
	! tcrms_f  	   => Report variable: the rms value of free MP toxic concentration at each reporting step			[kg m^-3       ]
	! tcmax_b 	   => Report variable: the maximum value of bonded MP toxic concentration at each reporting step	[kg m^-3       ]
	! tcmin_b  	   => Report variable: the minimum value of bonded MP toxic concentration at each reporting step	[kg m^-3       ]
	! tcrms_b 	   => Report variable: the minimum value of bonded MP toxic concentration at each reporting step	[kg m^-3       ]
	! aaxis		   => model parameter: the longest axis of grain
	! baxis		   => model parameter: the median axis of grain	   
	! caxis	 	   => model parameter: the shortest axis of grain
	! dequ         => model parameter: equivalent grain size 					          [mm                      ] 
	!				  (for plastic it will be estimated by the equivalent diameter a/b/c axis)
	! Prho         => model parameter: plastic density 							          [kg/m^3		           ]
	! tox_leach_k  => model parameter: bonded toxic releasing rate				          [1/s                     ]
	! tox_decay_k  => model parameter: free tox decaying rate                             [1/s                     ]
	! tox_gen_k    => model parameter: tox generation rate							      [1/s                     ]
	! tox_r	       => model parameter: mass ratio of bottom toxin comparing to MP mass    [-                       ]	
	! theta_s      => model parameter: equivalent boundary layer thickness of MP          [m                       ]
	! A_Hamaker    => model parameter: Hamaker constant                                   [J m^-2                  ]
	! ws_sed       => model parameter: sediment settling/rising velocity with in sediment [m/s                     ]
	! unit_factor  => unit factor for MP resuspension flux unit transition                [various                 ]
    !
	!              Notes units system: on all the input and output value of MP and toxin 'variables' is SI!	
	!		 	   All the computing concentration and fluxes for MP is transformed to mass or equivalent 
	!              mass, in kg m^-3 or kg m^-2. All the computing concentration and fluxes for toxins is 
	!              mass, in kg m^-3 or kg m^-2.
	!
	! conc         => plastic variables: plastics concentration in water column 	      [# m^-3 or kg m^-3       ] 	
	! cnew         => plastic variables: plastics concentration during update   	      [# m^-3 or kg m^-3       ]
	! mass         => plastic variables: plastics mass density in at the bottom           [kg m^-2                 ]
	! conc_sed     => plastic variables: plastics concentration in sediment               [# m^-3 or kg m^-3       ]
	! conc_sed_s   => plastic variables: plastics concentration stopped in sediment       [# m^-3 or kg m^-3       ]
	! conc_sed_p   => plastic variables: plastics concentration can pass in sediment      [# m^-3 or kg m^-3       ]
	! frac         => plastic variables: plastics fraction in active layer 	              [-                       ]
	! eflx         => plastic variables: suspended plastics erosive flux (+out of bed)    [# m^-2 or kg m^-2       ]
	! dflx         => plastic variables: suspended plastics depositional flux (+into bed) [# m^-2 or kg m^-2       ]
	! iflx         => plastic variables: suspended plastics infiltrative flux (+into bed) [# m^-2 or kg m^-2       ]
	! bflx         => plastic variables: suspended plastics burial flux (+into bed)       [# m^-2 or kg m^-2       ]
	! cdis         => plastic variables: concentration of MP at river source 		      [# m^-3 or kg m^-3       ]
	! cflx         => plastic variables: store advective flux at open boundry 		      [# m^-2 or kg m^-2	   ]
	! cobc         => plastic variables: user specd open boundry concentration 	          [# m^-3 or kg m^-3       ]
	! wset		   => plastic variables: spatial varying plastics settling velocity       [m/s                     ]
	! tau_ce       => plastic variables: critical shear stress for erosion 		          [N/m^2                   ]
	! erate        => plastic variables: bed surface erosion mass flux         	          [kg m^-1 s^-1            ]
	! prate 	   => plastic variables: free surface precipitation rate                  [# m^-3 s^-3             ]
	! zactv        => plastic variables: active layer thickness                           [m^-1                    ]
	! D_plast      => plastic variables: solid internal diffusivity of MP                 [m^2/s                   ]
	! D_toxin      => plastic variables: solid internal diffusivity of MP associated tox  [XXX]
	! lambda_plast => plastic variables: solid internal retention rate of MP              [XXX]
	! tconc_b	   => embedded toxic variables: embedded toxic concentration			  [kg m^-3                 ]
	! tmass_b	   => embedded toxic variables: bottom embedded toxic					  [kg m^-2		           ]
	! tconc_b_sed  => embedded toxic variables: total embedded toxic in sediment  		  [kg m^-3                 ]
	! tconc_b_sed_s=> embedded toxic variables: stopped embedded toxic in sediment		  [kg m^-3                 ]
	! tconc_b_sed_p=> embedded toxic variables: passing embedded toxic in sediment  	  [kg m^-3                 ]
	! teflx_b	   => embedded toxic variables: erosion flux of embedded toxic 		      [kg m^-2		           ]
	! tdflx_b	   => embedded toxic variables: deposition flux of embedded toxic 		  [kg m^-2		           ]
	! tiflx_b      => embedded toxic variables: embedded toxic flux infiltrated into bed  [kg m^-2		           ]
	! tiflx_b_deep => free toxic variables: free toxic flux infiltrated into deep bed     [kg m^-2                 ]
	! tcdis_b	   => embedded toxic variables: embedded toxic concentration at river 	  [kg m^-3                 ]
	! tcobc_b	   => embedded toxic variables: open boundary embedded toxic concentration[kg m^-3                 ]
	! tcflx_b      => embedded toxic variables: embedded toxic flux at river			  [kg m^-2		           ]
	! tfrac_b      => embedded toxic variables: fraction of embedded toxic in active layer[-                       ]
	! tconc_f      => free toxic variables: free toxic concentration					  [kg m^-3                 ]
	! tconc_f_sed  => free toxic variables: free toxic concentration in sediment	      [kg m^-3                 ]
	! tcdis_f	   => free toxic variables: free toxic concentration at river		      [kg m^-3                 ]
	! tcnew_f      => free toxic variables: free toxic concentration during update	      [kg m^-3                 ]
	! tcflx_f 	   => free toxic variables: free toxic flux at river				      [kg m^-2 	               ]
	! tcobc_f      => free toxic variables: open boundary free toxic concentration        [kg m^-3                 ]
	! tbdflx_f	   => free toxic variables: bottom free toxic diffusion flux			  [kg m^-2 	               ]
	! tiflx_f      => free toxic variables: free toxic flux infiltrated into bed          [kg m^-2                 ]
	! tiflx_f_deep => free toxic variables: free toxic flux infiltrated into deep bed     [kg m^-2                 ]
	! arraysize    => spatial size of plastic arrays
	!---------------------------------------------------------------------------------------------

	type plast_type
		character(len=20) :: pname
		character(len=20) :: MP_shape
		character(len=20) :: tname_b
		character(len=20) :: tname_f
		real(sp)		  :: aaxis
		real(sp)		  :: baxis
		real(sp)		  :: caxis
		real(dp)          :: cmax 
		real(dp)          :: cmin 
		real(dp)          :: crms
		real(dp)		  :: tcmax_f
		real(dp)		  :: tcmin_f
		real(dp)		  :: tcrms_f
		real(dp)		  :: tcmax_b
		real(dp)		  :: tcmin_b
		real(dp)		  :: tcrms_b
		real(sp)		  :: dequ 		
		real(sp)		  :: Prho
		real(sp)		  :: prate
		real(sp) 		  :: tox_leach_k
		real(sp) 		  :: tox_decay_k
		real(sp)          :: tox_gen_k
		real(sp) 		  :: tox_r
		real(sp)          :: theta_s
		real(sp)		  :: A_Hamaker
		real(sp)          :: ws_sed
		real(sp)          :: unit_factor
		real(sp), allocatable :: conc(:,:)
		real(sp), allocatable :: cnew(:,:)
		real(sp), allocatable :: mass(:)
		real(sp), allocatable :: conc_sed(:,:)
		real(sp), allocatable :: conc_sed_s(:,:)
		real(sp), allocatable :: conc_sed_p(:,:)
		real(sp), allocatable :: frac(:)
        real(sp), allocatable :: eflx(:)        
		real(sp), allocatable :: dflx(:)
		real(sp), allocatable :: iflx(:)
		real(sp), allocatable :: bflx(:)
		real(sp), allocatable :: pflx(:)
		real(sp), allocatable :: fflx(:)        
		real(sp), allocatable :: cdis(:)
		real(sp), allocatable :: cflx(:,:)		
		real(sp), allocatable :: cobc(:)
		real(sp), allocatable :: D_plast(:)
		real(sp), allocatable :: D_toxin(:)
		real(sp), allocatable :: lambda_plast(:,:)
		real(sp), allocatable :: tconc_b(:,:)
		real(sp), allocatable :: tcnew_b(:,:)
		real(sp), allocatable :: tconc_b_sed(:,:)
		real(sp), allocatable :: tconc_b_sed_s(:,:)
		real(sp), allocatable :: tconc_b_sed_p(:,:)
		real(sp), allocatable :: tcdis_b(:)
		real(sp), allocatable :: tcflx_b(:,:)
		real(sp), allocatable :: tmass_b(:)
		real(sp), allocatable :: tcobc_b(:)
		real(sp), allocatable :: teflx_b(:)
		real(sp), allocatable :: tiflx_b(:)
		real(sp), allocatable :: tdflx_b(:)
		real(sp), allocatable :: tfrac_b(:)
		real(sp), allocatable :: tconc_f(:,:)
		real(sp), allocatable :: tcnew_f(:,:)
		real(sp), allocatable :: tconc_f_sed(:,:)
		real(sp), allocatable :: tcdis_f(:)
		real(sp), allocatable :: tcflx_f(:,:)
		real(sp), allocatable :: tcobc_f(:)
		real(sp), allocatable :: tbdflx_f(:)
		real(sp), allocatable :: tiflx_f(:)
		real(sp), allocatable :: tiflx_f_deep(:)
		real(sp), allocatable :: wset(:,:) 
		real(sp), allocatable :: zactv(:)
		real(sp), allocatable :: erate(:)
		real(sp), allocatable :: tau_ce(:)
		integer           :: arraysize
	end type
 
	public

	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	!Public variables
	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	!--------------------------------------------------------------------------
	!Global Model Parameters
	!
	! plastic_model_version : Model version label
	! plastfile   	    : plastic input parameter control file
	! inf_bed   		: switch of infinite plastic supply
	! oned_model		: switch of 1D plastic module
	! plast_nudge		: switch of plastic OBC nudging
	! plast_source		: switch of plastic point sources
	! plast_bedset      : switch of plastic bed-setting
	! plast_infil		: switch of plastic infiltration
	! plast_hot_start	: switch of plastic restart
	! tox_vector		: switch of toxic-vector functional
	! tox_leaching      : switch of toxin decaying model, it can take the value of
    !					  1) SFOK: single first order kinetics
	!					  2) DFOK: double first order kinetics
	!					  3) SI-FOK: Shape-impacted (non-Newtonian) leaching with no reaction
	!					  4) SIR-FOK: Shape-impacted (non-Newtonian) leaching with reaction
	! cshields          : switch of shields formula used
	! plast_unit		: switch of plastic concentration unit
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
	! DTplast     		: plastic model time step (seconds)
	! qlim				: limiter parameter
	! T_model   		: model time (seconds)
	! taub      		: array to hold bottom shear stress
	! tauf				: array to hold free surface shear stress 
	! Wrho0 			: constant water density
	! tmass0            : constant initial bond toxin mass density on the bottom
	! mass0 			: constant initial plastic mass density on the bottom
	! sd50_0            : constant median grain size in the bottom 
	! dens0				: constant plastic density in the bottom
	! mudp0             : constant mud portion in the bottom
	! poro0             : constant porosity in the bottom
	! tort0             : tortuosity factor for uniform sediment [-] 
	! temp0             : temperature for uniform sediment [K]
	! darv0             : pore water velocity for uniform sediment [m/s]
	! rho_water			: seawater density can be constant or spatial varying
	! infil_max_depth	: maximum infiltration depth for uniform sediment [m]
	!-------------------------------------------------------------------------

	character(len=120), parameter :: PLAST_model_version = "FVCOM-MP_v0.5: toxin-leaching"
	character(len=120) :: plastfile
	character(len=120) :: plast_bedset
	character(len=120) :: tox_leaching
	character(len=20)  :: plast_unit
	character(len=20)  :: cshields
	logical  :: inf_bed
	logical  :: oned_model = .false.
	logical  :: plast_nudge = .false.
	logical  :: PLAST_dumpbot = .false.
	logical  :: PLAST_source = .false.
	logical  :: PLAST_infil = .false.
	logical  :: plast_hot_start = .false.
	logical  :: tox_vector = .false.
	integer  :: nplast
	integer  :: n_report
	integer  :: Nobc
	integer  :: Nobc_gl
	integer  :: PLAST_ramp  
	integer  :: PLAST_start
	integer  :: PLAST_its
	integer  :: settle_limiter = 2
	real(sp) :: DTplast
	real(sp) :: qlim  = 1.0
	real(sp) :: PLAST_alpha
	real(sp) :: settle_cfl = 1.0
	real(sp) :: T_model
	real(sp) :: Wrho0 = 1025
	real(sp) :: mass0
	real(sp) :: tmass0 = 0.0
	real(sp) :: sd50_0
	real(sp) :: dens0 = 2400
	real(sp) :: mudp0 = 0.1
	real(sp) :: poro0 = 0.5
	real(sp) :: tort0 = 1.0
	real(sp) :: temp0 = 298.15
	real(sp) :: darv0 = 0.0
	real(sp) :: infil_max_depth = 1.0
	real(sp), allocatable, target :: rho_water(:,:)
	real(sp), allocatable, target :: taub(:)
	real(sp), allocatable, target :: tauf(:)

	!-------------------------------------------------------------------------------!
	!  BOTTOM properties indices:                                          			!
	!																	   			!
	!  data --> bottom(horizontal,prototype)						       			!
	!  size --> bottom(0:mt,n_bot_vars)   				           		   			!
	!																	   			!
	! with proptype indices:											   			!
	!   isd50  =>  Median plastic grain diameter.				 (m)       			!
	!   idens  =>  Median plastic grain density. 				 (kg/m3)   			!
	!   imudp  =>  Mud Potion.  					  			 (-)       			!
	!   iporo  =>  plastic porosity 							 (-)       			!
	!   iTsed  =>  sediment temperature 						 (K)       			!
	!   iwDcy  =>  vertical Darcy speed 						 (m/s)     			!
	!   iLsed  =>  sediment layer thickness 					 (m)       			!
	!																	   			!
	! global bed characteristics:										   			!
	!   n_bot_vars:  number of bed characteristics, will be re-allocated if needed. !
	!   N_sed_layer: number of sediment layers								   		!
	!	bot_snames:   bottom characteristics name list					   			!
	!	bot_lnames:   bottom characteristics long name list				   			!
	!	bot_units: 	  bottom characteristics unit list                     			!
	!-------------------------------------------------------------------------------!

	integer, parameter    :: n_bot_vars  = 7
	integer, parameter    :: N_sed_layer = 3
	character(len=80), dimension(n_bot_vars) :: bot_snames = &
	(/"bot_sd50","bot_dens","bot_mudp","bot_poro","bot_Tsed","bot_wDcy","bot_Lsed"/)
	character(len=120), dimension(n_bot_vars) :: bot_lnames = &
	(/"bottom_sediment_d50", &
	"bottom_density", &
	"bottom_mud_portion", &
	"bottom_porosity", &
	"bottom_temperature", &
	"bottom_vertical_darcy_speed", &
	"bottom_layer_thickness"/)
	character(len=80), dimension(n_bot_vars) :: bot_units = &
	(/"m","kg/m^3","dimensionless","dimensionless","K","m/s","m"/)

	integer, parameter    :: isd50  = 1
	integer, parameter    :: idens  = 2
	integer, parameter    :: imudp  = 3
	integer, parameter    :: iporo  = 4
	integer, parameter    :: iTsed  = 5
	integer, parameter    :: iwDcy  = 6
	integer, parameter    :: iLsed  = 7

	real(sp), allocatable, target :: bottom(:,:)

	!--------------------------------------------------
	!Plastic Array                  
	!--------------------------------------------------

	type(PLAST_type), allocatable, target :: plast(:)

	contains

!==========================================================================
! Main submodules of FVCOM-MP modules:
!
! Setup_PLAST           - Initialize plastic model, allocate arrays, and 
!                         configure parameters
! Convert_PLAST_Params   - Convert plastic parameters from non-standard 
!                         units (mm) to MKS units
! Read_PLAST_Params      - Read plastic model parameters from input 
!                         configuration file
! Advance_plast          - Advance plastic transport model by one time step 
!                         including advection diffusion settling all source
!                          terms and kinetics
! Alloc_Plast_Vars       - Allocate memory for all plastic model variables 
!                         and arrays
! Report_PLAST_Setup     - Report plastic model setup parameters and initial 
!                         conditions to screen
! PLAST_Report           - Monitor and report plastic model statistics during 
!                         simulation
! Setup_PLAST_OBC        - Setup open boundary conditions and nudging for 
!                         plastic concentrations
!==========================================================================
  
Subroutine Setup_PLAST(fname,Nin,Nin_gl,N_plast,N_plast_Max,Plast_Names, &
		Toxb_Names,Toxf_Names,Tox_switch) 
	
	!==========================================================================
	! Allocate Data, Initialize Concentration Fields and Bed Parameters       
	!==========================================================================

	use control, only : msr,casename,input_dir
	use lims,    only : m, kbm1
	use Mod_Utils, only: pstop
	
	implicit none
	
	!Input/Output variables	
	integer, intent(in)  :: Nin,Nin_gl   
	character(len=80)    :: fname
	integer, intent(out) :: N_plast
	integer, intent(in)  :: N_plast_max
	logical, intent(out) :: Tox_switch	
	character(len=20), intent(inout) :: Plast_Names(N_plast_max)
	character(len=20), intent(inout) :: Toxb_Names(N_plast_max)
	character(len=20), intent(inout) :: Toxf_names(N_plast_max)
	
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
		if (plast_nudge) call setup_plast_obc 
		
	!----------------------------------------------------------------
	! report plastic parameters and statistics to screen
	!----------------------------------------------------------------
		call report_PLAST_setup

	!----------------------------------------------------------------
	! convert parameters to MKS unit system 
	!----------------------------------------------------------------
		call convert_PLAST_params 

	!----------------------------------------------------------------
	! pass some global variables to outside
	!----------------------------------------------------------------
	
		Tox_switch = tox_vector
		N_plast = nplast
		do i=1,nplast
			Plast_names(i) = plast(i)%pname
			if(tox_vector)then
				Toxb_names(i) = plast(i)%tname_b
				Toxf_names(i) = plast(i)%tname_f
			end if
		end do
		
	if(dbg_set(dbg_sbr)) write(ipt,*) "End: setup_plast" 

end Subroutine Setup_PLAST

Subroutine Convert_PLAST_Params 
	
	!=======================================================================
	!Convert Sed Params from nonstandard (mm, etc) to mks units
	!=======================================================================

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

Subroutine Read_PLAST_Params
  
  	!=======================================================================
	!Read plastic Parameters From plastic Input File
	!=======================================================================
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
  Call Get_Val(plast_hot_start,plastfile,'PLAST_HOT_START',line=linenum,echo=.false.)  

  !read in interation interval for reporting 
  Call Get_Val(n_report,plastfile,'N_REPORT',line=linenum,echo=.false.)  

  !read in logical switch for infinite plastic supply 
  Call Get_Val(inf_bed,plastfile,'INF_PLAST',line=linenum,echo=.false.)

  !read in nudging switch 
  Call Get_Val(plast_nudge,plastfile,'PLAST_NUDGE',line=linenum,echo=.false.)

  !read in params controlling output 
  Call Get_Val(PLAST_dumpbot,plastfile,'PLAST_DUMPBOT',line=linenum,echo=.false.)
  
  !read in plast bed settings
  Call Get_Val(plast_bedset,plastfile,'PLAST_BEDSET',line=linenum,echo=.false.)

  !read infiltration control switch
  Call Get_Val(plast_infil,plastfile,'PLAST_INFIL',line=linenum,echo=.false.)

  !read infiltration parameters only if enabled
  if (plast_infil) then
    ! Global infiltration parameters (uniform sediment properties)
    Call Get_Val(infil_max_depth,plastfile,'INFIL_MAX_DEPTH',line=linenum,echo=.false.) 
    Call Get_Val(tort0,plastfile,'SED_TORT',line=linenum,echo=.false.)
    Call Get_Val(temp0,plastfile,'SED_TEMP',line=linenum,echo=.false.)
    Call Get_Val(darv0,plastfile,'SED_DARV',line=linenum,echo=.false.)
	Call Get_Val(N_sed_layer,plastfile,'SED_LNUM',line=linenum,echo=.false.)

    ! Validate parameters
    if (infil_max_depth <= 0.0) then
        write(*,*)'error: INFIL_MAX_DEPTH must be > 0'
        write(*,*)'you entered: ',infil_max_depth
        stop
    endif
  endif

  !read in nudging relaxation factor 
  if(plast_nudge) then
  Call Get_Val(PLAST_alpha,plastfile,'PLAST_ALPHA',line=linenum,echo=.false.)
  Call Get_Val(PLAST_ramp ,plastfile,'PLAST_RAMP ',line=linenum,echo=.false.)
  endif
  
  !read in the Tox_vector switch
  Call Get_Val(tox_vector,plastfile,'TOX_VECTOR ',line=linenum,echo=.false.)
  
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
  
  !read in the constant median plastic diameter variation
  Call Get_Val(sd50_0,plastfile,'SED_SD50',line=linenum,echo=.false.)
  
  !read in the constant plastic density
  Call Get_Val(dens0,plastfile,'SED_SRHO',line=linenum,echo=.false.)
  
  !read in the constant mud portion
  Call Get_Val(mudp0,plastfile,'SED_MUDP',line=linenum,echo=.false.)
  
  !read in the constant plastic porosity
  Call Get_Val(poro0,plastfile,'SED_PORS',line=linenum,echo=.false.)
  
  if (tox_vector) then
	!read toxin leaching mode
	Call Get_Val(tox_leaching,plastfile,'TOX_MOD',line=linenum,echo=.false.)
  end if
  
  !check values
  if (nplast < 1) then
    write(*,*)'plastic input file must have at least one plastic class'
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
		
		Call Get_Val(stemp,plastfile,'PLAST_SHAP',line=linenum,echo=.false.,start=k1)
		plast(i)%MP_shape = stemp

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
		
		!read precipitation rate
		Call Get_Val(ftemp,plastfile,'PLAST_PRAT',line=linenum,echo=.false.,start=k1)
		plast(i)%prate = ftemp
		
		select case (plast_unit)
			case ('kg/m^3')
				plast(i)%unit_factor = 1.0
			case ('g/l')
				plast(i)%unit_factor = 1.0
			case ('#/m^3')
				plast(i)%unit_factor = 1.0/(plast(i)%Prho*plast(i)%aaxis*plast(i)%baxis*plast(i)%caxis*(10.0)**(-9.0))
			case ('#/l')
				plast(i)%unit_factor = 1.0/(1000.0*plast(i)%Prho*plast(i)%aaxis*plast(i)%baxis*plast(i)%caxis*(10.0)**(-9.0))
			case DEFAULT
				plast(i)%unit_factor = 1.0
		end select
		
		!Toxic-vector functional
		if (tox_vector) then
			
			!read the name of the toxic
			Call Get_Val(stemp,plastfile,'TOX_NAME',line=linenum,echo=.false.,start=k1)
			plast(i)%tname_b = TRIM(TRIM(stemp) // '_b')
			
			Call Get_Val(stemp,plastfile,'TOX_NAME',line=linenum,echo=.false.,start=k1)
			plast(i)%tname_f = TRIM(TRIM(stemp) // '_f')
  
			  select case (tox_leaching)

				case ('SFOK')
				
					!read constant bonded toxin releasing rates
					Call Get_Val(ftemp,plastfile,'TOX_CKL',line=linenum,echo=.false.,start=k1)
					plast(i)%tox_leach_k = ftemp
					
					!read free toxic decaying rate
					Call Get_Val(ftemp,plastfile,'TOX_CKD',line=linenum,echo=.false.,start=k1)
					plast(i)%tox_decay_k = ftemp
					
					!set bonded toxin generation rates
					Call Get_Val(ftemp,plastfile,'TOX_CKG',line=linenum,echo=.false.,start=k1)
					plast(i)%tox_gen_k = ftemp
				
				case ('SI-FOK')
				
					!Compute shape-impacted leaching factor
					Call Calc_Toxin_Leach_Rate_SI(plast(i)%MP_shape,plast(i)%aaxis,plast(i)%baxis,plast(i)%caxis,plast(i)%theta_s,plast(i)%D_plast,ftemp)
					plast(i)%tox_leach_k = ftemp
					
					!read free toxic decaying rate
					Call Get_Val(ftemp,plastfile,'TOX_CKD',line=linenum,echo=.false.,start=k1)
					plast(i)%tox_decay_k = ftemp
					
					!set bonded toxin generation rates
					Call Get_Val(ftemp,plastfile,'TOX_CKG',line=linenum,echo=.false.,start=k1)
					plast(i)%tox_gen_k = ftemp
				
				case ('SIR-FOK')
				
					!Under-development!!!
					!Compute shape-impacted leaching factor
					Call Calc_Toxin_Leach_Rate_SI(plast(i)%MP_shape,plast(i)%aaxis,plast(i)%baxis,plast(i)%caxis,plast(i)%theta_s,plast(i)%D_plast,ftemp)
					plast(i)%tox_leach_k = ftemp
					
					!read free toxic decaying rate
					Call Get_Val(ftemp,plastfile,'TOX_CKD',line=linenum,echo=.false.,start=k1)
					plast(i)%tox_decay_k = ftemp
					
					!set bonded toxin generation rates
					Call Get_Val(ftemp,plastfile,'TOX_CKG',line=linenum,echo=.false.,start=k1)
					plast(i)%tox_gen_k = ftemp
				
				case DEFAULT

					!Default is set to no leaching
					plast(i)%tox_leach_k = 0.0
					plast(i)%tox_decay_k = 0.0
					plast(i)%tox_gen_k = 0.0
						
			  end select
			
			!read initial ratio of toxic in the MP (mass-rato) on the bed
			Call Get_Val(ftemp,plastfile,'TOX_INR',line=linenum,echo=.false.,start=k1)
			plast(i)%tox_r = ftemp
			
		end if

		!read species-specific infiltration parameters if enabled
		if (plast_infil) then
		    ! Read species-specific infiltration parameters
		    Call Get_Val(ftemp,plastfile,'PLAST_AHAMA',line=linenum,echo=.false.,start=k1)
		    plast(i)%A_Hamaker = ftemp
		    
		    Call Get_Val(ftemp,plastfile,'PLAST_WSSED',line=linenum,echo=.false.,start=k1)
		    plast(i)%ws_sed = ftemp
		    
		end if
		
	end do
 
 if(dbg_set(dbg_sbr)) write(ipt,*) "End: Read_PLAST_Params " 
 
End Subroutine Read_PLAST_Params 
     
Subroutine Advance_plast(DTin,Tin,taub_in,s_in,t_in)  
  
    !=======================================================================
	! Advance plastic Model by time step DTplast 
	!=======================================================================

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
  real(sp), dimension(mt) :: usurf,vsurf

  integer :: i,k,iplast,l1,l2,ierr,d_cdis,d_cflx
  real(sp) :: fact,ufact,temp,s_temp(kbm1),t_temp(kbm1),rho_temp(kbm1),tc_temp(m,kb),dtc_temp

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Advance_Plast "

  !------------------------------------------------------
  !Check for Initialization
  !------------------------------------------------------
  if(iint < PLAST_start) return
  
  ! Compute surface velocity at node
  CALL E2N2D(u(:,1),usurf)
  CALL E2N2D(v(:,1),vsurf)

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
  DTplast = DTin
  T_model = Tin
  
  do i=1,m
	taub(i) = taub_in(i)*rho_water(i,kbm1)
	tauf(i) = (usurf(i)**2+vsurf(i)**2)*rho_water(i,1)
  end do

  !------------------------------------------------------
  !Plastic Dynamics Updating
  !------------------------------------------------------

  !Calculate Active Layer Thickness --> bottom(:,iactv)
  call calc_active_layer
    
  !calculated depositional flux from settling (settle flux function is called here)
  call calc_deposition

  !calculate erosive flux (resuspension function is called here)
  call calc_erosion

  !calculate infiltration and blockage
  if (PLAST_infil) then
	Call Calc_Infiltration
	Call Calc_Blockage
	do iplast=1,nplast
		do i=1,m
			do k=1,N_sed_layer
				plast(iplast)%conc_sed(i,k)=plast(iplast)%conc_sed_p(i,k)+plast(iplast)%conc_sed_s(i,k)
				plast(iplast)%tconc_b_sed(i,k)=plast(iplast)%tconc_b_sed_p(i,k)+plast(iplast)%tconc_b_sed_s(i,k)
			end do
		end do
	end do
  end if

  !adjust bottom properties
  call update_bottom

  !------------------------------------------------------
  !3D Advection of Plastic Concentration
  !------------------------------------------------------
  if(.not. oned_model)then
     d_cdis = 1
     if(numqbc > 0) d_cdis = numqbc
     d_cflx = 1
     if(Nobc > 0) d_cflx = Nobc
      
     do i=1,nplast
      
      !set discharge plastic concentrations
      if(numqbc > 0 .and. PLAST_source)then
        plast(i)%cdis(1:numqbc) = plastdis(1:numqbc,i)
      endif
  
      !3D advection
      call adv_scal(plast(i)%conc, &
                      plast(i)%cnew, &
                      d_cdis,      &
                      plast(i)%cdis(1:numqbc), &
                      d_cflx,      &
                      plast(i)%cflx, &
                      DTplast,PLAST_source)

	  !This concentration updating function is for plastic, no need to change here.
      !call fct_sed(plast(i)%conc,plast(i)%cnew) 
  
    end do
  else
    do i=1,nplast
      plast(i)%cnew = plast(i)%conc
    end do
  endif 

  !------------------------------------------------------
  !Vertical Diffusion of Plastic Concentrations
  !------------------------------------------------------
  do i=1,nplast
    call vdif_scal(plast(i)%cnew,DTplast) 
  end do
  
  if (tox_vector) then
  
	  !------------------------------------------------------
	  !Bonded and free toxic dynamics:
	  !------------------------------------------------------
	  
		!Horizontal Advection of bonded toxic Concentration 
		if(.not. oned_model)then
		  d_cdis = 1
		  if(numqbc > 0) d_cdis = numqbc
		  d_cflx = 1
		  if(Nobc > 0) d_cflx = Nobc
		  
		  do i=1,nplast
		  
		  !set discharge bonded toxic concentrations
		  if(numqbc > 0 .and. PLAST_source)then
			plast(i)%tcdis_b(1:numqbc) = toxbdis(1:numqbc,i)
		  endif
	  
		  !3D advection
		  call adv_scal(plast(i)%tconc_b, &
						  plast(i)%tcnew_b, &
						  d_cdis,      &
						  plast(i)%tcdis_b(1:numqbc), &
						  d_cflx,      &
						  plast(i)%tcflx_b, &
						  DTplast,PLAST_source)

		  !This concentration updating function is for plastic, no need to change here.
		  !call fct_sed(plast(i)%tconc_b,plast(i)%tcnew_b) 
	  
		  end do
		else
		  do i=1,nplast
			plast(i)%tcnew_b = plast(i)%tconc_b
		  end do
		endif 

		!Vertical Diffusion of bonded toxic Concentrations
		do i=1,nplast
			call vdif_scal(plast(i)%tcnew_b,DTplast) 
		end do
	  
		!Horizontal advection of free toxic concentration
		if(.not. oned_model)then
		 d_cdis = 1
		 if(numqbc > 0) d_cdis = numqbc
		 d_cflx = 1
		 if(Nobc > 0) d_cflx = Nobc
		  
			do i=1,nplast
		  
			!set discharge plastic concentrations
			if(numqbc > 0 .and. PLAST_source) then
				plast(i)%tcdis_f(1:numqbc) = toxfdis(1:numqbc,i)
			endif
	  
			!3D advection
			call adv_scal(plast(i)%tconc_f, &
						  plast(i)%tcnew_f, &
						  d_cdis,      &
						  plast(i)%tcdis_f(1:numqbc), &
						  d_cflx,      &
						  plast(i)%tcflx_f, &
						  DTplast,PLAST_source)

			!This concentration updating function is for plastic, no need to change here.
			!call fct_sed(plast(i)%tconc_f,plast(i)%tcnew_f) 
	  
			end do
		else
			do i=1,nplast
				plast(i)%tcnew_f = plast(i)%tconc_f
			end do
		endif 
	  
		!Vertical diffusion of toxic concentration
		do i=1,nplast
			call vdif_scal(plast(i)%tcnew_f,DTplast) 
		end do
		
		!Bonded and free toxic concentration updating
	do iplast=1,nplast
				
			!Leaching kinetic terms for bonded and free toxin in the water column
			do i=1,m
				do k=1,kb
					plast(iplast)%tcnew_f(i,k)=plast(iplast)%tcnew_f(i,k)+DTplast*plast(iplast)%tox_leach_k*plast(iplast)%tcnew_b(i,k)
					plast(iplast)%cnew(i,k)=plast(iplast)%cnew(i,k)-DTplast*plast(iplast)%tox_leach_k*plast(iplast)%tcnew_b(i,k)
					plast(iplast)%tcnew_b(i,k)=plast(iplast)%tcnew_b(i,k)-DTplast*plast(iplast)%tox_leach_k*plast(iplast)%tcnew_b(i,k)
				end do
			end do

			!Leaching kinetic terms for bonded and free toxin in the sediment column
			if (PLAST_infil) then
				do i=1,m
					do k=1,N_sed_layer
						plast(iplast)%tconc_f_sed(i,k)=plast(iplast)%tconc_f_sed(i,k)+DTplast*plast(iplast)%tox_leach_k*plast(iplast)%tconc_b_sed(i,k)
						plast(iplast)%conc_sed(i,k)=plast(iplast)%conc_sed(i,k)-DTplast*plast(iplast)%tox_leach_k*plast(iplast)%tconc_b_sed(i,k)
						plast(iplast)%conc_sed_s(i,k)=plast(iplast)%conc_sed_s(i,k)-DTplast*plast(iplast)%tox_leach_k*plast(iplast)%tconc_b_sed_s(i,k)
						plast(iplast)%conc_sed_p(i,k)=plast(iplast)%conc_sed_p(i,k)-DTplast*plast(iplast)%tox_leach_k*plast(iplast)%tconc_b_sed_p(i,k)
						plast(iplast)%tconc_b_sed(i,k)=plast(iplast)%tconc_b_sed(i,k)-DTplast*plast(iplast)%tox_leach_k*plast(iplast)%tconc_b_sed(i,k)
						plast(iplast)%tconc_b_sed_s(i,k)=plast(iplast)%tconc_b_sed_s(i,k)-DTplast*plast(iplast)%tox_leach_k*plast(iplast)%tconc_b_sed_s(i,k)
						plast(iplast)%tconc_b_sed_p(i,k)=plast(iplast)%tconc_b_sed_p(i,k)-DTplast*plast(iplast)%tox_leach_k*plast(iplast)%tconc_b_sed_p(i,k)
					end do
				end do
			end if
				
			!generation kinetic terms for bonded toxin in the water column
            if (plast(iplast)%tox_gen_k > 0.0) then
				do i=1,m
					do k=1,kb
						plast(iplast)%tcnew_b(i,k)=plast(iplast)%tcnew_b(i,k)*(1.0+DTplast*plast(iplast)%tox_gen_k)
					end do
				end do
			end if

			!generation kinetic terms for bonded toxin in the sediment column
			if (PLAST_infil) then
				if (plast(iplast)%tox_gen_k > 0.0) then
					do i=1,m
						do k=1,N_sed_layer
							plast(iplast)%tconc_f_sed(i,k)=plast(iplast)%tconc_f_sed(i,k)*(1.0+DTplast*plast(iplast)%tox_gen_k)
						end do
					end do
				end if
			end if
				
			!decay kinetic terms for free Toxin in the water column
           if (plast(iplast)%tox_decay_k > 0.0) then
			do i=1,m
				do k=1,kb
					plast(iplast)%tcnew_f(i,k)=plast(iplast)%tcnew_f(i,k)*(1.0-DTplast*plast(iplast)%tox_decay_k)
				end do
			end do
		   end if

			!decay kinetic terms for free Toxin in the sediment column
			if (PLAST_infil) then
				if (plast(iplast)%tox_decay_k > 0.0) then
					do i=1,m
						do k=1,N_sed_layer
							plast(iplast)%tconc_f_sed(i,k)=plast(iplast)%tconc_f_sed(i,k)*(1.0-DTplast*plast(iplast)%tox_decay_k)
						end do
					end do
				end do
			end if
		   
				
	end do

  end if

  !------------------------------------------------------
  !Exchange Concentration on Processor Boundaries
  !------------------------------------------------------
# if defined (MULTIPROCESSOR)
  if(par)then
    do i=1,nplast
      call aexchange(nc,myid,nprocs,plast(i)%cnew)
	  if (tox_vector) then
		call aexchange(nc,myid,nprocs,plast(i)%tcnew_f)
		call aexchange(nc,myid,nprocs,plast(i)%tcnew_b)
	  end if
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
                           DTplast,       &
                           temp) 
    end do
  
	if (tox_vector) then
	   do i=1,nplast 
		  call bcond_scal_OBC(plast(i)%tconc_b, &
                           plast(i)%tcnew_b, &
                           plast(i)%tcflx_b, &
                           plast(i)%tcobc_b, &
                           DTplast,       &
                           temp) 
		  call bcond_scal_OBC(plast(i)%tconc_f, &
                           plast(i)%tcnew_f, &
                           plast(i)%tcflx_f, &
                           plast(i)%tcobc_f, &
                           DTplast,       &
                           temp) 
	   enddo
	endif
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
	
	if (tox_vector) then
		do i=1,nplast
			call bcond_scal_PTsource(plast(i)%tconc_b, &
                               plast(i)%tcnew_b, &
                               plast(i)%tcdis_b)
			call bcond_scal_PTsource(plast(i)%tconc_f, &
                               plast(i)%tcnew_f, &
                               plast(i)%tcdis_f)
		end do
    end if
  endif

  !------------------------------------------------------
  ! Update Concentration Variables to next time level
  !------------------------------------------------------
  do i=1,nplast
    plast(i)%conc = plast(i)%cnew
  end do
  
  if (tox_vector) then
	 do i=1,nplast
		plast(i)%tconc_b = plast(i)%tcnew_b
		plast(i)%tconc_f = plast(i)%tcnew_f
	 end do	
  end if

  !------------------------------------------------------
  ! Limit to Non-Negative (modify to Venkat Limiter) 
  !------------------------------------------------------
  do i=1,nplast
    plast(i)%conc = max(plast(i)%cnew,0.0)
	!plast(i)%conc = min(plast(i)%cnew,10.0)
  end do
  
  if (tox_vector) then
	do i=1,nplast
		plast(i)%tconc_b = max(plast(i)%tcnew_b,0.0)
		plast(i)%tconc_f = max(plast(i)%tcnew_f,0.0)
	end do
  end if
  
  !------------------------------------------------------
  !Exchange Concentration/tau/rho on Processor Boundaries
  !------------------------------------------------------
# if defined (MULTIPROCESSOR)
  if(par)then
    do i=1,nplast
      call aexchange(nc,myid,nprocs,plast(i)%cnew)
	  call aexchange(nc,myid,nprocs,plast(i)%mass)
	  if (tox_vector) then
		call aexchange(nc,myid,nprocs,plast(i)%tcnew_b)
		call aexchange(nc,myid,nprocs,plast(i)%tcnew_f)
		call aexchange(nc,myid,nprocs,plast(i)%tmass_b)
	  end if
    end do
	call aexchange(nc,myid,nprocs,tauf)
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

Subroutine Alloc_Plast_Vars
  
	!=======================================================================
	! Allocate plastic Variables                      
	!=======================================================================

	use lims,     only: nt,mt,kb,kbm1
  
	implicit none
  
	integer i,tmp1,tmp2

	if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Alloc_PLAST_Vars "

	!ensure arrays have nonzero dimension
	tmp1 = max(Nobc+1,1)
	tmp2 = max(Nobc  ,1)

	!allocate suspended plastic arrays

	do i=1,nplast
		plast(i)%arraysize = mt+1

		! Allocate basic plastic transport variables
		allocate(plast(i)%conc(0:mt,kb  ))      ; plast(i)%conc = 0.0
		allocate(plast(i)%cnew(0:mt,kb  ))      ; plast(i)%cnew = 0.0
		allocate(plast(i)%mass(0:mt     ))      ; plast(i)%mass = 0.0
		allocate(plast(i)%frac(0:mt     ))      ; plast(i)%frac = 0.0
		allocate(plast(i)%eflx(0:mt     ))      ; plast(i)%eflx = 0.0
		allocate(plast(i)%dflx(0:mt     ))      ; plast(i)%dflx = 0.0
		allocate(plast(i)%pflx(0:mt     ))      ; plast(i)%pflx = 0.0
		allocate(plast(i)%fflx(0:mt     ))      ; plast(i)%fflx = 0.0
		allocate(plast(i)%cdis(500      ))      ; plast(i)%cdis = 0.0
		allocate(plast(i)%cflx(tmp1,kbm1))      ; plast(i)%cflx = 0.0
		allocate(plast(i)%cobc(tmp2     ))      ; plast(i)%cobc = 0.0
		allocate(plast(i)%wset(0:mt,kbm1))      ; plast(i)%wset = 0.0
		allocate(plast(i)%zactv(0:mt    ))      ; plast(i)%zactv = 0.0
		allocate(plast(i)%erate(0:mt    ))      ; plast(i)%erate = 0.0
		allocate(plast(i)%tau_ce(0:mt   ))      ; plast(i)%tau_ce = 0.0

		! Allocate toxic leaching arrays if enabled
		if (tox_vector) then
		    allocate(plast(i)%tconc_b(0:mt,kb))		; plast(i)%tconc_b = 0.0
			allocate(plast(i)%tcnew_b(0:mt,kb))		; plast(i)%tcnew_b = 0.0
			allocate(plast(i)%tcdis_b(500	))		; plast(i)%tcdis_b = 0.0
			allocate(plast(i)%tcflx_b(tmp1,kbm1))	; plast(i)%tcflx_b = 0.0
			allocate(plast(i)%tmass_b(0:mt	))		; plast(i)%tmass_b = 0.0
			allocate(plast(i)%tcobc_b(tmp2	))		; plast(i)%tcobc_b = 0.0
			allocate(plast(i)%teflx_b(0:mt	))		; plast(i)%teflx_b = 0.0
			allocate(plast(i)%tdflx_b(0:mt	))		; plast(i)%tdflx_b = 0.0
			allocate(plast(i)%tfrac_b(0:mt	))		; plast(i)%tfrac_b = 0.0
			allocate(plast(i)%tconc_f(0:mt,kb))		; plast(i)%tconc_f = 0.0
			allocate(plast(i)%tcnew_f(0:mt,kb))		; plast(i)%tcnew_f = 0.0
			allocate(plast(i)%tcdis_f(500	))		; plast(i)%tcdis_f = 0.0
			allocate(plast(i)%tcflx_f(tmp1,kbm1))	; plast(i)%tcflx_f = 0.0
			allocate(plast(i)%tcobc_f(tmp2  ))	    ;plast(i)%tcobc_f = 0.0
			allocate(plast(i)%tbdflx_f(0:mt ))		; plast(i)%tbdflx_f= 0.0
		endif

		! Allocate infiltration arrays if enabled
		if (plast_infil) then
		    ! Plastic/toxic concentration arrays (species-specific)
		    allocate(plast(i)%conc_sed(0:mt, N_sed_layer)); plast(i)%conc_sed = 0.0
			allocate(plast(i)%conc_sed_s(0:mt, N_sed_layer)); plast(i)%conc_sed_s = 0.0
		    allocate(plast(i)%conc_sed_p(0:mt, N_sed_layer)); plast(i)%conc_sed_p = 0.0
			if (tox_vector) then
		        allocate(plast(i)%tconc_b_sed(0:mt, N_sed_layer)); plast(i)%tconc_b_sed = 0.0
				allocate(plast(i)%tconc_b_sed_s(0:mt, N_sed_layer)); plast(i)%tconc_b_sed_s = 0.0
				allocate(plast(i)%tconc_b_sed_p(0:mt, N_sed_layer)); plast(i)%tconc_b_sed_p = 0.0
		        allocate(plast(i)%tconc_f_sed(0:mt, N_sed_layer)); plast(i)%tconc_f_sed = 0.0
		    endif
		    
		    ! Infiltration flux arrays (nodal)
		    allocate(plast(i)%iflx(0:mt)); plast(i)%iflx = 0.0
			allocate(plast(i)%D_plast(0:mt)); plast(i)%D_plast = 0.0
			allocate(plast(i)%lambda_plast(0:mt, N_sed_layer)); plast(i)%lambda_plast = 0.0
		    if (tox_vector) then
				allocate(plast(i)%D_toxin(0:mt)); plast(i)%D_toxin = 0.0
		        allocate(plast(i)%tiflx_b(0:mt)); plast(i)%tiflx_b = 0.0
		        allocate(plast(i)%tiflx_f(0:mt)); plast(i)%tiflx_f = 0.0
		        allocate(plast(i)%tiflx_f_deep(0:mt)); plast(i)%tiflx_f_deep = 0.0
		    endif
		endif
	end do

		!allocate bottom shear stress array 
		allocate(taub(0:mt)) ; taub = 0.0
		allocate(tauf(0:mt)) ; tauf = 0.0

		!allocate bottom array
		if (plast_infil) then
			n_bot_vars = 7
		else
			n_bot_vars = 4
		end if
		allocate(bottom(0:mt , n_bot_vars)) ; bottom = 0.0
		
		!allocate the water density array
		allocate(rho_water(0:mt,kbm1)); rho_water = 0.0

	if(dbg_set(dbg_sbr)) write(ipt,*) "End: Alloc_PLAST_Vars "

End Subroutine Alloc_PLAST_Vars

Subroutine Report_PLAST_Setup
  
	!==========================================================================
	! Report plastic Setup To Screen 
	!==========================================================================

	use control, only: msr 
  
	implicit none
  
	integer :: i

	if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Report_PLAST_Setup "

	! echo plastic model parameters to screen
	if(msr)then

		write(*,*)
		write(*,*)'------------------ plastic Model Setup --------------------'
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
		
		if(plast_hot_start)then
			write(*,*)'!  plast_hot_start         :',plast_hot_start
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
		
		if (tox_vector) then
			write(*,*) '!	toxic vector		: active'
		else
			write(*,*) '!	toxic vector		: not active'
		end if
				
		write(*,*)'|    class       Aaxis   Baxis   Caxis    Prho   |'
	
		do i=1,nplast
			write(*,'(1X,A1,1X,A12,1X,3(F7.4,1X),1X,F8.2,1X,A1)')'|', &
			trim(plast(i)%pname), plast(i)%aaxis,plast(i)%baxis,plast(i)%caxis,plast(i)%Prho,'|'
		end do
	
	end if

  ! echo plastic model initial conditions to screen
  call PLAST_report  

  write(*,*)'------------------------------------------------------------'
  write(*,*)

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Report_PLAST_Setup "

End Subroutine Report_PLAST_Setup

Subroutine PLAST_Report  
  
	!==========================================================================
	! Monitor plastic Model: Compute and Write Statistics to Screen
	!==========================================================================

  use control, only: msr ,par, mpi_fvcom_group
  use lims   , only: m,mt,kbm1, msrid
  
  implicit none
  
  real(dp) :: tau_min,tau_max,tmp1(3),tmp2(3),tmp3(3),tmp4(3),tmp5(3),tmp6(3)
  integer :: i,dim1,dim2,ierr

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: PLAST_Report "

  !set up limits
  dim1 = mt 
  dim2 = kbm1
  
  !calculate max/min/rms concentrations (this rms is actually mean)
  do i=1,nplast
    plast(i)%cmax = maxval(plast(i)%conc(1:dim1,1:dim2))
    plast(i)%cmin = minval(plast(i)%conc(1:dim1,1:dim2))
    plast(i)%crms = sum(dble(abs(plast(i)%conc(1:dim1,1:dim2))))
    plast(i)%crms = plast(i)%crms/(dim1*dim2)
	if (tox_vector) then
		plast(i)%tcmax_f = maxval(plast(i)%tconc_f(1:dim1,1:dim2))
		plast(i)%tcmin_f = minval(plast(i)%tconc_f(1:dim1,1:dim2))
		plast(i)%tcrms_f = sum(dble(abs(plast(i)%tconc_f(1:dim1,1:dim2))))
		plast(i)%tcrms_f = plast(i)%tcrms_f/(dim1*dim2)
		plast(i)%tcmax_b = maxval(plast(i)%tconc_b(1:dim1,1:dim2))
		plast(i)%tcmin_b = minval(plast(i)%tconc_b(1:dim1,1:dim2))
		plast(i)%tcrms_b = sum(dble(abs(plast(i)%tconc_b(1:dim1,1:dim2))))
		plast(i)%tcrms_b = plast(i)%tcrms_b/(dim1*dim2)
	end if
  end do

  tau_max   = maxval(taub(1:m))
  tau_min   = minval(taub(1:m))
 
# if defined (MULTIPROCESSOR)
  IF(PAR)THEN
  do i=1,nplast
    tmp1(1:3) = (/plast(i)%cmax,plast(i)%cmin,plast(i)%crms/)
    CALL MPI_REDUCE(tmp1,tmp2,3,MPI_DP,MPI_MAX,MSRID-1,MPI_FVCOM_GROUP,IERR)
    if (tox_vector) then
		tmp3(1:3) = (/plast(i)%tcmax_f,plast(i)%tcmin_f,plast(i)%tcrms_f/)
		CALL MPI_REDUCE(tmp3,tmp4,3,MPI_DP,MPI_MAX,MSRID-1,MPI_FVCOM_GROUP,IERR)
		tmp5(1:3) = (/plast(i)%tcmax_b,plast(i)%tcmin_b,plast(i)%tcrms_b/)
		CALL MPI_REDUCE(tmp5,tmp6,3,MPI_DP,MPI_MAX,MSRID-1,MPI_FVCOM_GROUP,IERR)
	end if	
	if(msr)then
      plast(i)%cmax = tmp2(1)
      plast(i)%cmin = tmp2(2)
      plast(i)%crms = tmp2(3)
	  if (tox_vector) then
	  plast(i)%tcmax_f = tmp4(1) 
	  plast(i)%tcmin_f = tmp4(2)
	  plast(i)%tcrms_f = tmp4(3)
	  plast(i)%tcmax_b = tmp6(1)
	  plast(i)%tcmin_b = tmp6(2)
	  plast(i)%tcrms_b = tmp6(3)
	  end if
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
  write(*,*  )'====================plastic Model Stats======================'
  write(*,*  )'!  quantity              :     avg(kg/m^3)      max(kg/m^3)      '
  do i=1,nplast
    write(*,100)'!',trim(plast(i)%pname),'    :',plast(i)%crms,plast(i)%cmax
	if (tox_vector) then
	write(*,100)'!',trim(plast(i)%tname_f),'  :',plast(i)%tcrms_f,plast(i)%tcmax_f
	write(*,100)'!',trim(plast(i)%tname_b),'  :',plast(i)%tcrms_b,plast(i)%tcmax_b
	end if
  end do
  write(*,*  )'!  max bottom stress     :     ',tau_max
  write(*,*  )'!  min bottom stress     :     ',tau_min 

  100 format(1x,a1,a20,a5,3f12.6)
  101 format(1x,a26,f12.6)

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: PLAST_Report "

End Subroutine PLAST_Report  

Subroutine Setup_PLAST_OBC

	!==========================================================================
	!  Setup open boundary nudging for plastic obc
	!   - check for existence of sed nudging file
	!   - read data (concentration of each sed type at each obc)
	!==========================================================================
  
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
  ! make sure plastic nudging data file exists 
  !------------------------------------------------------
  fname = "./"//trim(input_dir)//"/"//trim(casename)//"_plast_nudge.dat"

  inquire(file=trim(fname),exist=fexist)
  if(msr .and. .not.fexist)then
    write(*,*)'plastic nudging data: ',fname,' does not exist'
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
    write(*,*)'number of open boundary nodes for plastic nudging: ',tmp
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

!==========================================================================
! Water column plastic/toxic vertical transport submodules:
!
! Calc_Erosion            - Calculate plastic erosion flux from bed using 
!                         critical shear stress and erodibility parameters
! Calc_critial_shear_WS   - Calculate critical shear stress for plastic and 
!                         sediment erosion using Shields parameter
! Calc_Erosion_Rate       - Calculate erosion rate coefficients based on 
!                         bottom properties and shear stress conditions
! Calc_Deposition         - Calculate plastic depositional flux and update 
!                         concentrations from particle settling
! Calc_Wset_Plastic       - Calculate plastic particle settling velocity 
!                         based on size, density and water properties
! Settle_Flux_new         - Calculate settling flux using SLIP limiter for 
!                         plastic concentration transport in water column
! Settle_tFlux_new        - Calculate toxic settling flux using SLIP limiter 
!                         for bonded toxic concentration transport
! Lim                     - Calculate LED limiter function for numerical 
!                         flux limiting in advection schemes
!==========================================================================

Subroutine Calc_Erosion
    
	!==========================================================================
	! Calculate Erosion 
	! Calculate erosion process using the method of Ariathurai and Arulanandan (1978)
	!    -Erosion rates of cohesive soils. J. Hydraulics Division, ASCE, 104(2),279-282.-
	!    surface_mass_flux[i] = erate[i]*(1-porosity)*bfrac[i]*(tau_w/tau_crit[i] - 1)
	!      where
	!      iplast          = index of plastic type
	!      bfrac[i]        = fraction of bed composed by plastic i
	!      erate[i]        = bed erodibility constant
	!      porosity        = volume of voids/total volume in the top layer
	!      tau_w           = shear stress on the bed
	!      tau_ce          = critical shear stress of plastic
	!==========================================================================

  use all_vars
  
  implicit none
  
  integer  :: iplast ,i
  real(sp) :: bottom_por,bottom_frac,t_bottom_frac
  real(sp) :: bottom_mass,tau_f,tau_w,dep,active,erosion_flux,t_erosion_flux,precip_flux
  real(sp) :: tau_c_p(m)
  real(sp) :: tau_c_s(m)
  real(sp) :: Erate(m)
  
  real(sp) :: dx,dy 

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Calc_Erosion "

  !Compute Erosive Flux (kg/m^2), limited by available material in active layer 
  do iplast=1,nplast
	
	if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Critical shear stress computation"
		!Compute the critical stress for erosion
		Call Calc_critial_shear_WS(m,rho_water(1:m,kbm1),bottom(1:m,isd50),plast(iplast)%dequ, &
			dens0,plast(iplast)%Prho,tau_c_s,tau_c_p)
		do i=1,m
			plast(iplast)%tau_ce(i)=tau_c_p(i)
		end do
	if(dbg_set(dbg_sbr)) write(ipt,*) "End: Critical shear stress computation"
	
	if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Erosion coefficient"
		!Compute the sediment-like erosion rate based on mud-portion and grain size
		Call Calc_Erosion_Rate(m,taub,tau_c_p,tau_c_s,plast(iplast)%Prho,plast(iplast)%unit_factor, &
			plast(iplast)%mass,plast(iplast)%zactv,plast(iplast)%frac,Erate)
		do i=1,m
			plast(iplast)%erate(i)=Erate(i)
		end do
	if(dbg_set(dbg_sbr)) write(ipt,*) "End: Erosion coefficient"
  
	do i=1,m
# if defined (WET_DRY)
	  if(iswetn(i)/=1)cycle  ! no calculation when dry node
# endif
	  
		!Note that fractions of the toxic and plastic are computed inside the Update_Bottom routine
		tau_w = taub(i)
		tau_f = tauf(i)
		bottom_mass = plast(iplast)%mass(i)
		if (bottom_mass>0.0) then
			erosion_flux = DTplast*Erate(i)
			if (erosion_flux>bottom_mass) then
				erosion_flux = bottom_mass
			end if
		else
			erosion_flux = 0.0
		end if 
		precip_flux = plast(iplast)%prate*tau_f
		plast(iplast)%eflx(i) = erosion_flux 
		plast(iplast)%pflx(i) = precip_flux
		
		if (tox_vector) then
			t_bottom_frac = plast(iplast)%tfrac_b(i)
			t_erosion_flux = erosion_flux * t_bottom_frac
			plast(iplast)%teflx_b(i) = t_erosion_flux
		end if

    end do
  end do

  !Update Concentration in Bottom of Water Column
  do iplast=1,nplast
    do i=1,m

# if defined (WET_DRY)
      if(iswetn(i)==1) then
# endif
        
		plast(iplast)%conc(i,kbm1) = plast(iplast)%conc(i,kbm1) + plast(iplast)%eflx(i)/(d(i)*dz(i,kbm1))
		plast(iplast)%conc(i,1) = plast(iplast)%conc(i,1) + plast(iplast)%pflx(i)/(d(i)*dz(i,1))
		
		if (tox_vector) then
			plast(iplast)%tconc_b(i,kbm1) = plast(iplast)%tconc_b(i,kbm1) + plast(iplast)%teflx_b(i)/(d(i)*dz(i,kbm1)) 
			plast(iplast)%tconc_f(i,kbm1) = plast(iplast)%tconc_f(i,kbm1) + plast(iplast)%tbdflx_f(i)/(d(i)*dz(i,kbm1))
		end if		

# if defined (WET_DRY)
      else
        plast(iplast)%conc(i,kbm1) = 0.0_sp
		plast(iplast)%conc(i,1) = 0.0_sp
      end if
# endif

    end do
  end do
  
  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Calc_Erosion "

End Subroutine Calc_Erosion

Subroutine Calc_critial_shear_WS(m,rho_wb,dsed,dequ,rho_s,rho_p,tau_c_s,tau_c_p)

	!==========================================================================
    ! Calculate critical shear stress used for both plastic and sediment
    !==========================================================================

	implicit none
	
	!Input variables
	integer, intent(in)		  :: m
	real(sp), intent(inout)   :: tau_c_p(m)
	real(sp), intent(inout)   :: tau_c_s(m)
	real(sp), intent(in) 	  :: rho_wb(m) !Bottom water density
	real(sp), intent(in)      :: dsed(m)
	real(sp), intent(in)	  :: rho_s
	real(sp), intent(in)      :: rho_p
	real(sp), intent(in)      :: dequ
  
	!Local variables
	integer  :: i
	real :: R_p			   !Relative plastic density
	real :: R_s   		   !Relative sediment density
	real :: Dstar_p		   !Plastic Dimensionless diameter
	real :: Dstar_s        !Sediment Dimensionless diameter
	real :: theta_sc	   !Shields parameter for sediment
	real :: theta_pc       !Shields parameter for plastic
	
	real, parameter :: nu=1E-6	       !Water viscosity at 20deg
	real, parameter :: g=9.81		   !gravitational acceleration
	real, parameter :: HE_coef=-0.503   !Hide-Expose coefficient
	real, parameter :: Pre_coef=0.5588   !Pre-coefficient

	do i=1,m
		
		select case (plast_bedset)
		
			case ('levitation')
				!just give some none zero value here
				tau_c_p(i)=0.05
				tau_c_s(i)=0.1
		
			case ('concrete_bottom')
				!Calculate relative density
				R_p = abs((rho_p-rho_wb(i))/rho_wb(i))
	
				!Calculate Dimensionless diameter
				Dstar_p=(R_p*(g/nu**2))**0.3333*dequ
	
				select case (cshields)
					case ('R93')
						!Calculate Shields parameter for sediment using Van Rijn 1993
						if (Dstar_p>=1 .and. Dstar_p<4) then
							theta_pc=0.24*Dstar_p**(-1)
						elseif (Dstar_p>=4 .and. Dstar_p<10) then
							theta_pc=0.14*Dstar_p**(-0.64)
						elseif (Dstar_p>=10 .and. Dstar_p<20) then
							theta_pc=0.04*Dstar_p**(-0.1)
						elseif (Dstar_p>=20 .and. Dstar_p<150) then
							theta_pc=0.013*Dstar_p**0.39
						elseif (Dstar_p>=150) then
							theta_pc=0.055
						end if
					case ('SW97')
						!Calculate Shields parameter for sediment using Soulsby and Whitehouse 1997
						theta_pc=0.30/(1.0+1.2*Dstar_p)+0.055*(1-exp(-0.020*Dstar_p))
					case DEFAULT
						theta_pc=0.30/(1.0+1.2*Dstar_p)+0.055*(1-exp(-0.020*Dstar_p))
				end select
	
				!Calculate plastic critical shear stress tau_c_p
				tau_c_p(i)=g*(rho_p-rho_wb(i))*theta_pc*dequ
				
				!Just set the sediment critical shear stress tau_c_s take the output place
				tau_c_s(i)=tau_c_p(i)
			
			case ('active_bottom')
				!Calculate relative density
				R_s = abs((rho_s-rho_wb(i))/rho_wb(i))
	
				!Calculate Dimensionless diameter
				Dstar_s=(R_s*(g/nu**2))**0.3333*dsed(i)
	
				select case (cshields)
					case ('R93')
						!Calculate Shields parameter for sediment using Van Rijn 1993
						if (Dstar_s>=1 .and. Dstar_s<4) then
							theta_sc=0.24*Dstar_s**(-1)
						elseif (Dstar_s>=4 .and. Dstar_s<10) then
							theta_sc=0.14*Dstar_s**(-0.64)
						elseif (Dstar_s>=10 .and. Dstar_s<20) then
							theta_sc=0.04*Dstar_s**(-0.1)
						elseif (Dstar_s>=20 .and. Dstar_s<150) then
							theta_sc=0.013*Dstar_s**0.39
						elseif (Dstar_s>=150) then
							theta_sc=0.055
						end if
					case ('SW97')
						!Calculate Shields parameter for sediment using Soulsby and Whitehouse 1997
						theta_sc=0.30/(1.0+1.2*Dstar_s)+0.055*(1-exp(-0.020*Dstar_s))
					case DEFAULT
						theta_pc=0.30/(1.0+1.2*Dstar_p)+0.055*(1-exp(-0.020*Dstar_p))
				end select
	
				!Calculate Shields parameter for plastic
				theta_pc=Pre_coef*(dequ/dsed(i))**HE_coef*theta_sc
	
				!Calculate tau_c_p
				tau_c_p(i)=g*(rho_p-rho_wb(i))*theta_pc*dequ
			
				!Calculate tau_c_s
				tau_c_s(i)=g*(rho_s-rho_wb(i))*theta_sc*dsed(i)
			
			case ('bedload_bottom')
				!UNDER DEVELOPMENT....just give some none zero value here
				tau_c_p(i)=0.05
				tau_c_s(i)=0.1
			
			case DEFAULT
				!just give some none zero value here
				tau_c_p(i)=0.05
				tau_c_s(i)=0.1
			
		end select
		
	end do
	
End Subroutine Calc_critial_shear_WS

Subroutine Calc_Erosion_Rate(m,taub,tau_c_p,tau_c_s,rho_p,ufact,mass,zactv,frac,Erate)
	
	!==========================================================================
	! Calculate the erosion rate base on Perera et al. 2020 
	!==========================================================================

	implicit none
	
	!Input/Output variables
	integer, intent(in)		  :: m
	real(sp), intent(inout)   :: Erate(m)
	real(sp), intent(in)      :: taub(m)
	real(sp), intent(in) 	  :: tau_c_p(m)
	real(sp), intent(in)      :: tau_c_s(m)
	real(sp), intent(in)      :: rho_p
	real(sp), intent(in)      :: ufact
	real(sp), intent(in)      :: mass(m)
	real(sp), intent(in)      :: frac(m)
	real(sp), intent(in)      :: zactv(m)

	!Local variables
	integer :: i
	real(sp) :: rho_pb
    real(sp) :: rho_eq
	real(sp) :: cmud
	real, parameter :: ks0=7.6E-5 !pure-sand erosion coefficient
	
	select case (plast_bedset)
		case ('levitation')
			do i=1,m
				Erate(i)=0.
			end do
		case ('concrete_bottom')
			do i=1,m
				Erate(i)=MAX(rho_p*ks0*1.035*ufact*(taub(i)-tau_c_p(i)),0.0)
			end do
		case ('active_bottom')
			cmud = 1.035*exp(-19.87*mudp0**0.86)+0.0981*mudp0**0.15
			
			do i=1,m
				rho_pb = mass(i)/((1.0-poro0)*zactv(i))
				rho_eq = dens0 - (dens0/rho_p-1)*rho_pb
				if (taub(i)<tau_c_p(i)) then
					Erate(i)=0.
				elseif (taub(i)>=tau_c_p(i) .and. taub(i)<tau_c_s(i)) then
					Erate(i)=rho_eq*ks0*cmud*frac(i)*ufact*(taub(i)**2-tau_c_p(i)*taub(i)-tau_c_p(i)*tau_c_s(i)+tau_c_p(i)**2)/(taub(i)-tau_c_p(i))
				elseif (taub(i)>=tau_c_s(i)) then
					Erate(i)=rho_eq*ks0*cmud*frac(i)*ufact*(taub(i)-tau_c_p(i))
				end if
			end do
		case ('bedload_bottom')
		    ! UNDER DEVELOPMENT....
			do i=1,m
				Erate(i)=0.
			end do
		case DEFAULT
			do i=1,m
				Erate(i)=0.
			end do
	end select

End Subroutine Calc_Erosion_Rate

Subroutine Calc_Deposition
  
  !==========================================================================
  ! Calculate Depositional Flux and Update Concentration from Settling
  !==========================================================================
  
  use all_vars
  
  implicit none
  
  integer  :: iplast ,i,mcyc,ncyc,k
  real(sp) :: wset(kbm1),c(kbm1),ctb(kbm1),dx(kbm1),dflux,fflux,tfluxb,DTmax,eps,DTdep,da,db,dc,rhop

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Calc_Deposition "

  eps = epsilon(eps)

  !Loop over plastic classes and spatial points       
  do iplast=1,nplast
    do i=1,m
  
      ! no calculation when dry nodes

# if defined (WET_DRY)
      if (iswetn(i)/=1) cycle  
# endif

      !initialize flux
      dflux = 0.0
	  fflux = 0.0
	  if (Tox_vector) then
		tfluxb = 0.0
	  end if

      !setup 1D concentration and grid arrays
      c(1:kbm1) = plast(iplast)%conc(i,1:kbm1)
	  if (Tox_vector) then
		ctb(1:kb) = plast(iplast)%tconc_b(i,1:kbm1)
	  end if
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
		mcyc = int(DTplast/DTmax + 1. - eps)
		DTdep = DTplast/float(mcyc)

		!call flux-limited settling equation
		do ncyc = 1,mcyc
			call settle_flux_new(kbm1,c,dx,wset,dflux,fflux,DTdep)
		end do

		!store solution in 3D array
		plast(iplast)%conc(i,1:kbm1) = c(1:kbm1)
	  
		!store wset in 3D array
		plast(iplast)%wset(i,1:kbm1) = wset(1:kbm1)

		!store depositional flux
		plast(iplast)%dflx(i) = dflux
		
		!store floating flux
		plast(iplast)%fflx(i) = fflux
		
		!Bonded Toxic Deposition
		if (Tox_vector) then
			do ncyc = 1,mcyc
				call settle_tflux_new(kbm1,ctb,dx,wset,tfluxb,DTdep)
			end do

			!store solution in 3D array
			plast(iplast)%tconc_b(i,1:kbm1) = ctb(1:kbm1)

			!store depositional flux
			plast(iplast)%tdflx_b(i) = tfluxb
		end if
		
	  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Settling advection"	  
	  
	  !shutdown erosion if node is inactive 
      !plast(iplast)%dflx(i) = plast(iplast)%dflx(i)
    end do
  end do

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Calc_Deposition "

End Subroutine Calc_Deposition

Subroutine Calc_Wset_Plastic(n,rhow,rhop,da,db,dc,wset) 

	!==========================================================================
	! Calculate Settling Flux and Update plastic Concentration in Water Column 
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
  
Subroutine Settle_Flux_new(n,c,dx,wset,dflux,fflux,deltat) 

	!==========================================================================
	! Calculate Settling Flux and Update plastic Concentration in Water Column 
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

	   implicit none
	  integer , intent(in   ) :: n
		  real(sp), intent(inout) ::  c(n)
		  real(sp), intent(in   ) :: dx(n)
		  real(sp), intent(in   ) :: wset(n)
		  real(sp), intent(inout) :: dflux, fflux 
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

	  !surface bcs (extrapolate/zero flux)
	  cin(0)  =  -cin(1) 
	  cin(-1) =  -cin(1)
	  win(0)  =  -wset(1)
	  win(-1) =  -wset(1)

	  !bottom bcs (extrapolate/zero flux)
	  select case (plast_bedset)
		case ('levitation')
			cin(n+1) = -cin(n) 
			cin(n+2) = -cin(n)
			win(n+1) = -win(n)
			win(n+2) = -win(n)
		case ('concrete_bottom')
			if (win(n)>0.0) then
				cin(n+1) = cin(n) 
				cin(n+2) = cin(n)
				win(n+1) = win(n)
				win(n+2) = win(n)
			else
				cin(n+1) = -cin(n) 
				cin(n+2) = -cin(n)
				win(n+1) = -win(n)
				win(n+2) = -win(n)
			end if	
		case ('active_bottom')
			if (win(n)>0.0) then
				cin(n+1) = cin(n) 
				cin(n+2) = cin(n)
				win(n+1) = win(n)
				win(n+2) = win(n)
			else
				cin(n+1) = -cin(n) 
				cin(n+2) = -cin(n)
				win(n+1) = -win(n)
				win(n+2) = -win(n)
			end if	
		case ('bedload_bottom')
			if (win(n)>0.0) then
				cin(n+1) = cin(n) 
				cin(n+2) = cin(n)
				win(n+1) = win(n)
				win(n+2) = win(n)
			else
				cin(n+1) = -cin(n) 
				cin(n+2) = -cin(n)
				win(n+1) = -win(n)
				win(n+2) = -win(n)
			end if	
		case DEFAULT
			if (win(n)>0.0) then
				cin(n+1) = cin(n) 
				cin(n+2) = cin(n)
				win(n+1) = win(n)
				win(n+2) = win(n)
			else
				cin(n+1) = -cin(n) 
				cin(n+2) = -cin(n)
				win(n+1) = -win(n)
				win(n+2) = -win(n)
			end if
	  end select

	  !flux computation
	  do i=1,n+1
		wvel    = .5*(win(i)+win(i-1))  !settle velocity at interface
		dis4    = abs(wvel)/2.
		conv(i) = wvel*(cin(i)+cin(i-1))/2. 
		diss(i) = dis4*(cin(i)-cin(i-1)-lim(cin(i+1)-cin(i),cin(i-1)-cin(i-2))) 
	  end do

	  !re-zero out boundary flux
	  conv(1) = 0. ; diss(1) = 0.
	  
	  select case (plast_bedset)
		case ('levitation')
			conv(n+1) = 0. ; diss(n+1) = 0.
		case ('concrete_bottom')
			if (win(n)>0.0) then
				conv(n+1) = conv(n+1) ; diss(n+1) = diss(n+1)
			else
				conv(n+1) = 0. ; diss(n+1) = 0.
			end if
		case ('active_bottom')
			if (win(n)>0.0) then
				conv(n+1) = conv(n+1) ; diss(n+1) = diss(n+1)
			else
				conv(n+1) = 0. ; diss(n+1) = 0.
			end if
		case ('bedload_bottom')
			if (win(n)>0.0) then
				conv(n+1) = conv(n+1) ; diss(n+1) = diss(n+1)
			else
				conv(n+1) = 0. ; diss(n+1) = 0.
			end if
		case DEFAULT
			if (win(n)>0.0) then
				conv(n+1) = conv(n+1) ; diss(n+1) = diss(n+1)
			else
				conv(n+1) = 0. ; diss(n+1) = 0.
			end if
	  end select

	  !update
	  do i=1,n
		c(i) = cin(i) + (deltat/dx(i))*(-conv(i+1)+conv(i) + diss(i+1)-diss(i)) 
	  end do

	  !set bottom flux (> 0 = deposition)
	  dflux = dflux + deltat*(conv(n+1)-diss(n+1))
	  fflux = fflux + deltat*(conv(1)-diss(1))

	  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Settle_Flux "

End Subroutine Settle_Flux_new
 
Subroutine Settle_tFlux_new(n,c,dx,wset,dflux,deltat) 

	!==========================================================================
	! Calculate Settling Flux and Update plastic Concentration in Water Column 
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

	implicit none
	  integer , intent(in   ) :: n
		  real(sp), intent(inout) ::  c(n)
		  real(sp), intent(in   ) :: dx(n)
		  real(sp), intent(in   ) :: wset(n)
		  real(sp), intent(inout) :: dflux
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

	  !surface bcs (extrapolate/zero flux)
	  cin(0)  =  -cin(1) 
	  cin(-1) =  -cin(1)
	  win(0)  =  -wset(1)
	  win(-1) =  -wset(1)

	  !bottom bcs (extrapolate/zero flux)
	  select case (plast_bedset)
		case ('levitation')
			cin(n+1) = -cin(n) 
			cin(n+2) = -cin(n)
			win(n+1) = -win(n)
			win(n+2) = -win(n)
		case ('concrete_bottom')
			if (win(n)>0.0) then
				cin(n+1) = cin(n) 
				cin(n+2) = cin(n)
				win(n+1) = win(n)
				win(n+2) = win(n)
			else
				cin(n+1) = -cin(n) 
				cin(n+2) = -cin(n)
				win(n+1) = -win(n)
				win(n+2) = -win(n)
			end if	
		case ('active_bottom')
			if (win(n)>0.0) then
				cin(n+1) = cin(n) 
				cin(n+2) = cin(n)
				win(n+1) = win(n)
				win(n+2) = win(n)
			else
				cin(n+1) = -cin(n) 
				cin(n+2) = -cin(n)
				win(n+1) = -win(n)
				win(n+2) = -win(n)
			end if	
		case ('bedload_bottom')
			if (win(n)>0.0) then
				cin(n+1) = cin(n) 
				cin(n+2) = cin(n)
				win(n+1) = win(n)
				win(n+2) = win(n)
			else
				cin(n+1) = -cin(n) 
				cin(n+2) = -cin(n)
				win(n+1) = -win(n)
				win(n+2) = -win(n)
			end if	
		case DEFAULT
			if (win(n)>0.0) then
				cin(n+1) = cin(n) 
				cin(n+2) = cin(n)
				win(n+1) = win(n)
				win(n+2) = win(n)
			else
				cin(n+1) = -cin(n) 
				cin(n+2) = -cin(n)
				win(n+1) = -win(n)
				win(n+2) = -win(n)
			end if
	  end select

	  !flux computation
	  do i=1,n+1
		wvel    = .5*(win(i)+win(i-1))  !settle velocity at interface
		dis4    = abs(wvel)/2.
		conv(i) = wvel*(cin(i)+cin(i-1))/2. 
		diss(i) = dis4*(cin(i)-cin(i-1)-lim(cin(i+1)-cin(i),cin(i-1)-cin(i-2))) 
	  end do

	  !re-zero out boundary flux
	  conv(1) = 0. ; diss(1) = 0.
	  
	  select case (plast_bedset)
		case ('levitation')
			conv(n+1) = 0. ; diss(n+1) = 0.
		case ('concrete_bottom')
			if (win(n)>0.0) then
				conv(n+1) = conv(n+1) ; diss(n+1) = diss(n+1)
			else
				conv(n+1) = 0. ; diss(n+1) = 0.
			end if
		case ('active_bottom')
			if (win(n)>0.0) then
				conv(n+1) = conv(n+1) ; diss(n+1) = diss(n+1)
			else
				conv(n+1) = 0. ; diss(n+1) = 0.
			end if
		case ('bedload_bottom')
			if (win(n)>0.0) then
				conv(n+1) = conv(n+1) ; diss(n+1) = diss(n+1)
			else
				conv(n+1) = 0. ; diss(n+1) = 0.
			end if
		case DEFAULT
			if (win(n)>0.0) then
				conv(n+1) = conv(n+1) ; diss(n+1) = diss(n+1)
			else
				conv(n+1) = 0. ; diss(n+1) = 0.
			end if
	  end select

	  !update
	  do i=1,n
		c(i) = cin(i) + (deltat/dx(i))*(-conv(i+1)+conv(i) + diss(i+1)-diss(i)) 
	  end do

	  !set bottom flux (> 0 = deposition)
	  dflux = dflux + deltat*(conv(n+1)-diss(n+1))

	  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Settle_Flux "

End Subroutine Settle_tFlux_new

Function Lim(a,b)
	  
	!==========================================================================
	! Calculate LED Limiter L(u,v)  
	!==========================================================================
	  
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
! Sediment column plastic/toxic vertical transport submodules:
!
! Calc_Infiltration       - Calculate plastic infiltration into sediment bed 
!                         layers with transport and retention processes
! Calc_Retention          - Calculate plastic retention coefficients and 
!                         filtration processes in sediment layers
! Update_Bottom           - Update surface bed properties due to erosion, 
!                         deposition infiltration and mass flux changes
! Calc_Plastic_Fraction   - Calculate plastic mass fraction in active layer 
!                         for resuspension flux calculations
! Calc_Toxic_Fraction     - Calculate toxic mass fraction relative to plastic 
!                         mass for resuspension processes
! Settle_Flux_sed_new     - Calculate settling flux with infiltration and 
!                         burial terms for sediment bed interactions
! Calc_Active_Layer       - Calculate active layer thickness using Harris and 
!                         Wiberg method based on bottom shear stress
!==========================================================================

Subroutine Calc_Infiltration

	use all_vars
  
  	implicit none

	integer  :: iplast ,i,mcyc,ncyc,k
	real(sp) :: wset_sed_plastic(N_sed_layer), wset_sed_toxic(N_sed_layer), c_sed_p(N_sed_layer), ctb_sed_p(N_sed_layer), ctf_sed(N_sed_layer)
	real(sp) :: dx(N_sed_layer), iflx, bflx, itflxb, itflxf, btflxb, btflxf, DTmax, DTdep
	real(sp) :: eps
	logical(sp) :: isfloat_toxic

	eps = epsilon(eps)

	if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Calc_Infiltration "

	!Loop over plastic classes and spatial points for advection    
  	do iplast=1,nplast
    	do i=1,m
		
		!initialize flux
      	iflx = 0.0
	  	bflx = 0.0
		if (Tox_vector) then
			itflxb = 0.0
			itflxf = 0.0
			btflxb = 0.0
			btflxf = 0.0
		end if

		!setup 1D concentration and grid arrays
		c_sed_p(1:N_sed_layer) = plast(iplast)%conc_sed_p(i,1:N_sed_layer)
		if (Tox_vector) then
			ctb_sed_p(1:N_sed_layer) = plast(iplast)%tconc_b_sed_p(i,1:N_sed_layer)
			ctf_sed(1:N_sed_layer) = plast(iplast)%tconc_f_sed(i,1:N_sed_layer)
		end if
		dx(1:N_sed_layer) = infil_max_depth/N_sed_layer

		!Read in settling velocity 
# if defined (WET_DRY)
      	if (iswetn(i)/=1) then
			wset_sed_plastic(1:N_sed_layer) = ws_sed
		else
			wset_sed_plastic(1:N_sed_layer) = ws_sed+darv0
		end if

		if (iswetn(i)/=1) then
			wset_sed_toxic(1:N_sed_layer) = 0.0
			isfloat_toxic = .true.
		else
			wset_sed_toxic(1:N_sed_layer) = darv0
			isfloat_toxic = .false.
		end if
# else
		wset_sed_plastic(1:N_sed_layer) = ws_sed+darv0
		wset_sed_toxic(1:N_sed_layer) = darv0
		isfloat_toxic = .false.
# endif

		!Set first sediment layer plastic/toxin concentration = bottom water plastic/toxin concentration
		c_sed_p(1) = plast(iplast)%conc_sed_p(i,kbm1)
		if (Tox_vector) then
			ctb_sed_p(1) = plast(iplast)%tconc_b(i,kbm1)
			ctf_sed(1) = plast(iplast)%tconc_f(i,kbm1)
		end if

		if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Settling advection in sediment"
			!set up cycles (use max(CFL) == 1)
			DTmax = Settle_CFL*minval(dx(1:N_sed_layer)/wset_sed_plastic(1:N_sed_layer))
			mcyc = int(DTplast/DTmax + 1. - eps)
			DTdep = DTplast/float(mcyc)

			!call flux-limited settling equation
			do ncyc = 1,mcyc
				call settle_flux_sed_new(N_sed_layer,c_sed_p,dx,wset_sed_plastic,iflx,bflx,DTdep)
			end do

			!store solution in 3D array
			plast(iplast)%conc_sed_p(i,1:N_sed_layer) = c_sed_p(1:N_sed_layer)

			!store infiltration flux
			plast(iplast)%iflx(i) = iflx

			!store burial flux
			plast(iplast)%bflx(i) = bflx

			if (Tox_vector) then
				do ncyc = 1,mcyc
					call settle_flux_sed_new(N_sed_layer,ctb_sed_p,dx,wset_sed_plastic,itflxb,btflxb,DTdep)
					if (.not. isfloat_toxic) then
						! Handle free toxic case
						call settle_flux_sed_new(N_sed_layer,ctf_sed,dx,wset_sed_toxic,itflxf,btflxf,DTdep)
					end if
				end do

				!store toxic solution in 3D array
				plast(iplast)%tconc_b_sed_p(i,1:N_sed_layer) = ctb_sed_p(1:N_sed_layer)
				plast(iplast)%tconc_f_sed(i,1:N_sed_layer) = ctf_sed(1:N_sed_layer)

				!store toxic infiltration flux
				plast(iplast)%tiflx_b(i) = itflxb
				plast(iplast)%tiflx_f(i) = itflxf

				!store toxic burial flux
				plast(iplast)%tiflx_b_deep(i) = btflxb
				plast(iplast)%tiflx_f_deep(i) = btflxf

			end if

		if(dbg_set(dbg_sbr)) write(ipt,*) "End: Settling advection in sediment"

		end do
	end do
			
	if(dbg_set(dbg_sbr)) write(ipt,*) "End: Calc_Infiltration "

End Subroutine Calc_Infiltration

Subroutine Calc_Blockage

	use all_vars
  
  	implicit none

	integer  :: iplast ,i,k,w
	real	 :: 

# if defined (WET_DRY)
      	if (iswetn(i)/=1) then
			w = abs(ws_sed)
		else
			w = abs(ws_sed+darv0)
		end if

		if (iswetn(i)/=1) then
			w = 0.0
		else
			w = abs(darv0)
		end if
# else
		w = abs(ws_sed+darv0)
		w = abs(darv0)
# endif

	if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Calc_Infiltration "

		!Calculation of lambda
		Call Calc_Blockage_Rate

		!Loop over plastic classes and spatial points for kinetics updating    
		do iplast=1,nplast
			do i=1,m
				do k=1,N_sed_layer
					plast(iplast)%conc_sed_p(i,k)=plast(iplast)%conc_sed_p(i,k)-DTplast*plast(iplast)%lambda_plast(i)*w*plast(iplast)%conc_sed_p(i,k)
					plast(iplast)%conc_sed_s(i,k)=plast(iplast)%conc_sed_s(i,k)+DTplast*plast(iplast)%lambda_plast(i)*w*plast(iplast)%conc_sed_p(i,k)
					plast(iplast)%tconc_b_sed_p(i,k)=plast(iplast)%tconc_b_sed_p(i,k)-DTplast*plast(iplast)%lambda_plast(i)*w*plast(iplast)%tconc_b_sed_p(i,k)
					plast(iplast)%tconc_b_sed_s(i,k)=plast(iplast)%tconc_b_sed_s(i,k)+DTplast*plast(iplast)%lambda_plast(i)*w*plast(iplast)%tconc_b_sed_p(i,k)
				end do 
			end do
		end do

	if(dbg_set(dbg_sbr)) write(ipt,*) "End: Calc_Infiltration "

End subroutine Calc_Blockage

Subroutine Update_Bottom
	
	!==========================================================================
	! Update Surface Bed Properties Due to ReSuspension Load Flux
	! Change:
	!   plast(:)%mass(:,1)   to reflect increased/decrease mass (kg)
	!   bed(:,1,itchk)     to reflect changed bed thickness properties
	!==========================================================================

  use lims, only: m
  implicit none
  real(sp) :: dz,accum,flux,tflux_b,eps
  real(sp) :: tmass_b(m),tbdflx_f(m)
  integer  :: i,iplast
  real(sp) :: afrac(m),bfrac(m)

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Update_Bottom"

	  eps = epsilon(eps)
	  
	  !Updating the surface bed mass (mass deposition and infiltrating)
	  if (PLAST_infil) then
	
		if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Mass depositing and infiltrating"

		do i=1,m
			do iplast=1,nplast
			flux = (plast(iplast)%eflx(i)+plast(iplast)%iflx(i)-plast(iplast)%dflx(i))

			!change top layer bed mass of each plastic type 
			accum  = plast(iplast)%mass(i) - flux 
			plast(iplast)%mass(i) = max(accum,0.0) 
			end do
		end do

		if (tox_vector) then
			do i=1,m
				do iplast=1,nplast
					tflux_b = (plast(iplast)%teflx_b(i)+plast(iplast)%tiflx_b(i)-plast(iplast)%tdflx_b(i))
					
					!change top layer bed mass of each plastic type
					accum = plast(iplast)%tmass_b(i) - tflux_b
					plast(iplast)%tmass_b(i) = max(accum,0.0)
				end do
			end do
		endif

		if(dbg_set(dbg_sbr)) write(ipt,*) "End: Mass depositing and infiltrating"

	  else

		if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Mass depositing"
		
		do i=1,m
			do iplast=1,nplast
			flux = (plast(iplast)%eflx(i)-plast(iplast)%dflx(i))

			!change top layer bed mass of each plastic type 
			accum  = plast(iplast)%mass(i) - flux 
			plast(iplast)%mass(i) = max(accum,0.0) 
			end do
		end do
		
		if (tox_vector) then
			do i=1,m
				do iplast=1,nplast
					tflux_b = (plast(iplast)%teflx_b(i)-plast(iplast)%tdflx_b(i))
					
					!change top layer bed mass of each plastic type
					accum = plast(iplast)%tmass_b(i) - tflux_b
					plast(iplast)%tmass_b(i) = max(accum,0.0)
				end do
			end do
		endif

		if(dbg_set(dbg_sbr)) write(ipt,*) "End: Mass depositing"

	  end do
	  
	  !Updating the surface fractions
	  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Fraction computation"
	  
		 !Compute fraction based on active layer thickness
		 Call Calc_plastic_fraction
	  
		 if (tox_vector) then
			!Compute fraction of toxic fraction
			do iplast=1,nplast
				Call Calc_toxic_fraction(plast(iplast)%mass,plast(iplast)%tmass_b,bfrac)
				do i=1,m
					plast(iplast)%tfrac_b(i)=bfrac(i)
				end do
			end do
		 endif
		 
	  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Fraction computation"
  
	  !Updating the bottom mass of Toxic
	  if (tox_vector) then
		if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Toxic mass updating"
				
			!do iplast=1,nplast
					
				!Kinetics on bottom embedded toxin mass
				!  Call Calc_Linear_Bottom_Kinetics(DTplast,m,plast(iplast)%tox_leach_k, &
				!		plast(iplast)%tmass_b,tmass_b)
				!   do i=1,m
				!	  plast(iplast)%tbdflx_f(i)=tmass_b(i)-plast(iplast)%tmass_b(i)
				!	  plast(iplast)%tmass_b(i)=tmass_b(i)
				!   end do
							
				!Kinetics on free toxin mass
				!   Call Calc_Linear_Bottom_Kinetics(DTplast,m,-plast(iplast)%tox_gen_k, &
				!		plast(iplast)%tmass_b,tmass_b)
				!   do i=1,m
				!	  plast(iplast)%tmass_b(i)=tmass_b(i)
				!   end do
					
			!end do

		if(dbg_set(dbg_sbr)) write(ipt,*) "End: Toxic mass updating"
	  end if
  
  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Update_Bottom"

End Subroutine Update_Bottom

Subroutine Calc_Plastic_Fraction
	  
	  !==========================================================================
	! Calculate the Plastic Fraction for resuspension flux:
	!	The method plastic resuspension treat the plastic as a portion of 
	!sediment. The potion is by the plastic mass/total sediment + plastic mass 
	!in the active layer.
	!==========================================================================
	  
	  use lims,     only: m
	  implicit none
	
	  !Local variables
	  integer  :: i, iplast
	  real(sp) :: rho_pb
      real(sp) :: rho_eq
	  
	  select case (plast_bedset)
		case ('levitation')
		
			do iplast=1,nplast
			 do i=1,m
				plast(iplast)%frac(i) = 1.0
			 end do
			end do
		
		case ('concrete_bottom')
			
			do iplast=1,nplast
			 do i=1,m
				plast(iplast)%frac(i) = 1.0
			 end do
			end do
		
		case ('active_bottom')
			   
			do iplast=1,nplast
			 do i=1,m
				rho_pb = plast(iplast)%mass(i)/((1.0-poro0)*plast(iplast)%zactv(i))
				rho_eq = dens0 - (dens0/plast(iplast)%dequ-1)*rho_pb
				plast(iplast)%frac(i) = rho_pb/rho_eq
			 end do
			end do
			  
		case ('bedload_bottom')
		
			do iplast=1,nplast
			 do i=1,m
				plast(iplast)%frac(i) = 1.0
			 end do
			end do
		
		case DEFAULT
		
			do iplast=1,nplast
			 do i=1,m
				plast(iplast)%frac(i) = 1.0
			 end do
			end do
	  
	  end select

End Subroutine Calc_Plastic_Fraction

Subroutine Calc_Toxic_Fraction(mass_p,mass_t,bfrac)
	
	!==========================================================================
	! Calculate the Toxic Fraction for resuspension flux:
	!	The toxic fraction is a mass fraction of toxic vs the plastics.
	!==========================================================================
	
	use lims,     only: m
	implicit none
	
	!Input/Output variable
	real(sp), intent(inout) :: bfrac(m)
	real(sp), intent(in) 	:: mass_p(m)		!Plastic mass
	real(sp), intent(in)	:: mass_t(m)		!Toxic mass
	
	!Local variables
	integer :: i

	do i = 1, m
      !if (mass_p(i) > 0.0) then
      !  bfrac(i) = mass_t(i) / mass_p(i)
      !else
        bfrac(i) = 0.0
      !end if
    end do

End Subroutine Calc_Toxic_Fraction

Subroutine Calc_Active_Layer 

	!==========================================================================
	! Update Active Layer Thickness 
	!   Use method of Harris and Wiberg (1997):
	!   -Approaches to quantifying long-term continental shelf plastic transport
	!   -with an example from the northern California STRESS mid-shelf site.
	!   [Continental Shelf REsearch, 17, 1389-1418]
	!
	!   z_a = max[k_1*(tau_bot - tau_ce)*rho_0 , 0] + k_2 * D_50
	!       where
	!       k_1 = .007 (empirical constant)
	!       k_2 = 6.0  (empirical constant)
	!       tau_bot = shear stress on the bed
	!       tau_ce  = critical shear stress for erosion of the bed (averaged over
	!                 plastic classes)
	!       D_50    = median grain diameter of the surface plastic
	!       note: our shear stresses are nondimensionalized by rho, so
	!             rho does not appear in actual calculation (below)
	!==========================================================================

	use lims,     only: m
	implicit none
	
	!Local variables
	integer  :: i,iplast
	real(sp)     :: A, Ap
	real(sp)     :: B, Bp
	real(sp)     :: C, Cp
	real(sp)     :: Delta, Deltap
	real(sp)     :: zstar
	
	real(sp), parameter :: alpha_e=6.0    !coefficient for the shape factor nc
	real(sp), parameter :: alpha_f=0.007  !coefficient for the shape factor nc
	
	do iplast=1,Nplast
		do i=1,m
			
			select case (plast_bedset)
				case ('levitation')
					plast(iplast)%zactv(i) = 0.0
				case ('concrete_bottom')
					plast(iplast)%zactv(i) = 0.0
				case ('active_bottom')
				
				if (taub(i)<plast(iplast)%tau_ce(i)) then
					zstar = plast(iplast)%mass(i)/(plast(iplast)%Prho*(1-poro0))
					
					A = 1.0 - poro0
					B = -(alpha_e*(1.0 - poro0)*sd50_0)
					C = -plast(iplast)%mass(i)/(plast(iplast)%Prho)*(plast(iplast)%dequ-sd50_0)*alpha_e
					Delta = B**2.0-4.0*A*C
					
					if (Delta>0.0) then
						plast(iplast)%zactv(i)=(-B+Delta**0.5)/(2*A)
					else
						plast(iplast)%zactv(i)=alpha_e*sd50_0/2.0
					end if
					
					if (plast(iplast)%zactv(i)<zstar) then
						plast(iplast)%zactv(i)=alpha_e*plast(iplast)%dequ
					end if
					
				else
					zstar = plast(iplast)%mass(i)/(plast(iplast)%Prho*(1-poro0))	
					
					Ap = 1.0 - poro0
					Bp = -(alpha_e*(1.0 - poro0)*sd50_0+alpha_f*(taub(i)-plast(iplast)%tau_ce(i)))
					Cp = -plast(iplast)%mass(i)/(plast(iplast)%Prho)*(plast(iplast)%dequ-sd50_0)*alpha_e
					Deltap = Bp**2.0-4.0*Ap*Cp
					
					if (Delta>0.0) then
						plast(iplast)%zactv(i)=(-B+Delta**0.5)/(2*A)
					else
						plast(iplast)%zactv(i)=alpha_e*sd50_0/2.0
					end if
					
					if (plast(iplast)%zactv(i)<zstar) then
						plast(iplast)%zactv(i)=alpha_e*plast(iplast)%dequ
					end if
					
				end if
				
				case ('bedload_bottom')
					!UNDER DEVELOPMENT....
					plast(iplast)%zactv(i) = 0.0

				case DEFAULT
					plast(iplast)%zactv(i) = 0.0
			end select
	
		end do
	end do
End Subroutine Calc_Active_Layer

Subroutine Calc_Mechanical_Diffusivity(r_tracer,D_mech)

	use lims,     only: m
  
  	implicit none

	integer  :: i, j
	real(sp), intent(in) :: r_tracer
	real(sp), intent(inout) :: D_mech(:)
	real(sp)     :: w,tau,eta
	real(sp)     :: D_m
	real(sp)     :: T_norm, rho_norm, mu0, delta_mu, sum_term
	real(sp), parameter :: c=0.1   			!Emperical free parameter for mechanical diffusivity in solid matrix
	real(sp), parameter :: k_B=1.380649E-23 !Boltzmann constant
	real(sp), parameter :: PI=3.14159265359 !PI
	real(sp), parameter :: mu_star = 1.0E-6     ! Reference viscosity [Pa⋅s] [IAPWS 2008 viscosity formulation parameters (Huber et al.)]
	real(sp), parameter :: T_star = 647.067     ! Critical temperature [K] [IAPWS 2008 viscosity formulation parameters (Huber et al.)]
	real(sp), parameter :: rho_star = 322.0     ! Critical density [kg/m3] [IAPWS 2008 viscosity formulation parameters (Huber et al.)]
	real(sp), parameter :: H0(0:3) = [1.67752, 2.20462, 0.6366564, -0.241605] ! H0 coefficients [IAPWS 2008 viscosity formulation parameters (Huber et al.)]
	

# if defined (WET_DRY)
      	if (iswetn(i)/=1) then
			w = abs(ws_sed)
		else
			w = abs(ws_sed+darv0)
		end if

		if (iswetn(i)/=1) then
			w = 0.0
		else
			w = abs(darv0)
		end if
# else
		w = abs(ws_sed+darv0)
		w = abs(darv0)
# endif

	if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Calc_Mechanical_Diffusivity "

		do i=1,m

			!Tortuosity computation []
			tau=bottom(i,iporo)**(0.5)

			!Dynamical viscosity of water [IAPWS 2008 viscosity formulation Huber et al.]
			! Assuming 1 atm pressure and using temperature from bottom(i,iTsed)
				
			T_norm = bottom(i,iTsed) / T_star
			rho_norm = 1000.0 / rho_star  ! Water density at 1 atm ≈ 1000 kg/m³
				
			! Calculate mu0 (viscosity in the limit of zero density)
			sum_term = 0.0
			do j = 0, 3
				sum_term = sum_term + H0(j) / (T_norm**j)
			end do
			mu0 = 100.0 * sqrt(T_norm) / sum_term
				
			! For liquid water at 1 atm, delta_mu contribution is small, approximated as 1.0
			delta_mu = 1.0
				
			! Final viscosity calculation [Pa⋅s]
			eta = mu_star * mu0 * delta_mu

			!Molecular diffusivity at given temperature computation [Stokes-Einstein treatment]
			D_m = (k_B*bottom(i,iTsed))/(6*PI*eta*r_tracer)

			!Mechanical diffusivity estimation
			D_mech(i)=c*w*bottom(i,isd50)+D_m/tau

		end do

	if(dbg_set(dbg_sbr)) write(ipt,*) "End: Calc_Mechanical_Diffusivity "

End Subroutine Calc_Mechanical_Diffusivity

Subroutine Calc_Blockage_Rate
	
	use lims,     only: m
  
  	implicit none

	integer  :: iplast ,i, j, k
	real(sp)     :: lambda_DLVO, lambda_strain, alpha_primary, w, eps
	real(sp)     :: eta, eta_D, eta_I, eta_G, gamma, phi, d_g, d_MP
	real(sp)     :: A_S, N_R, N_Pe, N_vdW, N_A, N_G
	real(sp)     :: k_str, psi_str, dz
	real(sp)     :: Z_level(1:N_sed_layer)
	real(sp), parameter :: k_B=1.380649E-23 !Boltzmann constant
	real(sp), parameter :: PI=3.14159265359 !PI
	real(sp), parameter :: mu_star = 1.0E-6     ! Reference viscosity [Pa⋅s] [IAPWS 2008 viscosity formulation parameters (Huber et al.)]
	real(sp), parameter :: T_star = 647.067     ! Critical temperature [K] [IAPWS 2008 viscosity formulation parameters (Huber et al.)]
	real(sp), parameter :: rho_star = 322.0     ! Critical density [kg/m3] [IAPWS 2008 viscosity formulation parameters (Huber et al.)]
	real(sp), parameter :: H0(0:3) = [1.67752, 2.20462, 0.6366564, -0.241605] ! H0 coefficients [IAPWS 2008 viscosity formulation parameters (Huber et al.)]

	eps = epsilon(eps)

# if defined (WET_DRY)
      	if (iswetn(i)/=1) then
			w = abs(ws_sed)
		else
			w = abs(ws_sed+darv0)
		end if

		if (iswetn(i)/=1) then
			w = 0.0
		else
			w = abs(darv0)
		end if
# else
		w = abs(ws_sed+darv0)
		w = abs(darv0)
# endif

	if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Calc_Blockage_Rate "

		do iplast=1,nplast
			
			!Read MP diameter
			d_MP=plast(iplast)%dequ

			!Read Hamaker parameter
			AH=plast(iplast)%A_Hamaker

			!Estimate the alpha_primary
			if (A_Hamaker>0.0) then
				alpha_primary = 1.0
			else
				alpha_primary = 0.0
			end if

			!Compute diffusivity
			r_tracer = plast(iplast)%dequ/2.0
			Call Calc_Mechanical_Diffusivity(r_tracer,plast(iplast)%D_plast)
			
			do i=1,m

				!Establish level vector
				dz=bottom(i,iLsed)/N_sed_layer
				do k=1, N_sed_layer
					Z_level(k)=dz/2.0+(k-1)*dz
				end do

				!Dynamical viscosity of water [IAPWS 2008 viscosity formulation Huber et al.]
				! Assuming 1 atm pressure and using temperature from bottom(i,iTsed)
				T_norm = bottom(i,iTsed) / T_star
				rho_norm = 1000.0 / rho_star  ! Water density at 1 atm ≈ 1000 kg/m³
					
				! Calculate mu0 (viscosity in the limit of zero density)
				sum_term = 0.0
				do j = 0, 3
					sum_term = sum_term + H0(j) / (T_norm**j)
				end do
				mu0 = 100.0 * sqrt(T_norm) / sum_term
					
				! For liquid water at 1 atm, delta_mu contribution is small, approximated as 1.0
				delta_mu = 1.0
					
				! Final viscosity calculation [Pa⋅s]
				eta = mu_star * mu0 * delta_mu

				!Dimensionless parameter computation fro lambda_DLVO
				d_g=bottom(i,isd50)
				phi=bottom(i,iporo)
				gamma=(1-phi)**(0.33333333)
				A_S=2*(1-gamma**5.0)/(2.0-3.0*gamma+3.0*gamma**5.0-2.0*gamma**6.0)
				N_R=d_MP/d_g
				N_Pe=abs(darv0)*d_g/plast(iplast)%D_plast(i)
				N_G=abs(ws_sed)/(abs(darv0)+eps)
				N_vdW=AH/(k_B*bottom(i,iTsed))
				N_A=AH/(12.0*PI*eta*(d_MP/2)**2.0*abs(darv0+eps))

				!Retention coefficient computation
				eta_D=2.4*A_S**(0.33333333)*N_R**(-0.081)*N_Pe**(-0.715)*N_vdW**(0.052)
				eta_I=0.55*A_S*N_R**(1.675)*N_A**(0.125)
				eta_G=0.22*N_R^(-0.24)*N_G**(1.11)*N_vdW**(0.053)

				!lambda_DLVO computation
				lambda_DLVO = 3*(1-phi)/(2*d_g)*(eta_D+eta_I+eta_G)*alpha_primary

				!lambda_strain computation and combined lambda computation
				k_str=269.7*(d_MP/d_g)**(1.42)

				do k=1,N_sed_layer
					psi_str=((d_MP+Z_level(k))/d_MP)**(-0.43)
					lambda_strain=k_str*psi_str/w
					plast(iplast)%lambda_plast(i,k)=lambda_strain+lambda_DLVO
				end do

			end do
		end do

	if(dbg_set(dbg_sbr)) write(ipt,*) "End: Calc_Blockage_Rate "

End Subroutine Calc_Blockage_Rate

Subroutine Settle_Flux_sed_new(n,c,dx,wset,iflux,bflux,deltat) 

	!==========================================================================
	! Calculate Settling Flux and Update plastic Concentration in Sediment Column 
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

	  implicit none
	  integer , intent(in   ) :: n
		  real(sp), intent(inout) ::  c(n)
		  real(sp), intent(in   ) :: dx(n)
		  real(sp), intent(in   ) :: wset(n)
		  real(sp), intent(inout) :: dflux, fflux 
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

	  !surface bcs (extrapolate/zero flux)
	  if (win(1)>0.0) then
			cin(0)  =  cin(1) 
			cin(-1) =  cin(1)
			win(0)  =  wset(1)
			win(-1) =  wset(1)
	  else
			cin(0)  =  -cin(1) 
			cin(-1) =  -cin(1)
			win(0)  =  -wset(1)
			win(-1) =  -wset(1)
	  end if

	  !bottom bcs (extrapolate/zero flux)
	  if (win(n)>0.0) then
			cin(n+1) = cin(n) 
			cin(n+2) = cin(n)
			win(n+1) = win(n)
			win(n+2) = win(n)
	  else
			cin(n+1) = -cin(n) 
			cin(n+2) = -cin(n)
			win(n+1) = -win(n)
			win(n+2) = -win(n)
	  end if

	  !flux computation
	  do i=1,n+1
		wvel    = .5*(win(i)+win(i-1))  !settle velocity at interface
		dis4    = abs(wvel)/2.
		conv(i) = wvel*(cin(i)+cin(i-1))/2. 
		diss(i) = dis4*(cin(i)-cin(i-1)-lim(cin(i+1)-cin(i),cin(i-1)-cin(i-2))) 
	  end do

	  !re-set boundary flux
	  if (win(1)>0.0) then
		  conv(1) = conv(1) ; diss(1) = diss(1)
	  else
		  conv(1) = 0. ; diss(1) = 0.
	  end if

	  if (win(n)>0.0) then
		  conv(n+1) = conv(n+1) ; diss(n+1) = diss(n+1)
	  else
		  conv(n+1) = 0. ; diss(n+1) = 0.
	  end if

	  !update
	  do i=1,n
		c(i) = cin(i) + (deltat/dx(i))*(-conv(i+1)+conv(i) + diss(i+1)-diss(i)) 
	  end do

	  !set bottom flux (> 0 = deposition)
	  iflux = iflux + deltat*(conv(1)-diss(1))
	  bflux = bflux + deltat*(conv(n+1)-diss(n+1))

	  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Settle_Flux "

End Subroutine Settle_Flux_sed_new

!==========================================================================
! Plastic-toxic interaction related submodules:
!
! Calc_Toxin_Leach_Rate_SI    - Calculate shape-impacted toxin leaching rate 
!                              for microplastic particles with geometric factors
! Calc_Linear_Bottom_Kinetics - Calculate first-order kinetic decay or growth 
!                              terms for bottom plastic and toxic processes
!==========================================================================

Subroutine Calc_Toxin_Leach_Rate_SI(sname,a,b,c,theta,D_p,k)

	!==========================================================================
	! Toxin leaching rate: In this subroutine, the toxin leaching rate for the
	! given MP particle is computed with shape factor considered.
	!========================================================================== 

	Implicit None
	
	!Input variables
	real(sp), intent(inout) :: k
	character(len=20), intent(in)   :: sname
	real(sp), intent(in) 	:: a
	real(sp), intent(in) 	:: b
	real(sp), intent(in) 	:: c
	real(sp), intent(in) 	:: theta
	real(sp), intent(in) 	:: D_p
	
	!Local variables
	real			 :: Bim
	real			 :: Bim_crit
	real			 :: L_p
    real		     :: h_flux
	real             :: e_ratio
	real             :: lambda_1_sq
	real			 :: f_L,f_M,f_S
	real			 :: phi_1,phi_2
	real, parameter  :: Pie=3.14159265358979	!pi
	real, parameter  :: D_w=5E-10  				!Water diffusivity
	real, parameter  :: alpha_1=2.4048  		!first zero of Bessel function J0
	real, parameter  :: factor_ellipsoid=3.5	!factor computed from Rayleigh quotient
	real, parameter  :: p=1.6075
	real, parameter  :: A_sphere=2.55
	real, parameter  :: A_cylinder=23.1323
	real, parameter  :: B_cylinder=9.8696
	real, parameter  :: C_cylinder=2.0
	real, parameter  :: D_cylinder=2.0
	real, parameter  :: gamma_cylinder=0.6667
	real, parameter  :: A_cuboid=4.9348
	real, parameter  :: gamma_cuboid=0.3333
	real, parameter  :: A_ellipsoid=1.1667
	real, parameter  :: gamma_ellipsoid=1.3333
	
	L_p = (a*b*c)**(1.0/3.0)
	Bim = (D_w/theta)*(L_p/2.0)/D_p
	
	!Large Bim Range 
	if (Bim >= 1E3) then
	
		select case (sname)
			case ('sphere')
				lambda_1_sq=Pie**2.0
			case ('cylinder')
				lambda_1_sq=4.0*alpha_1**(2.0)*(a/b)**(2.0/3.0)+Pie**(2.0)*(b/a)**(4.0/3.0)
			case ('cuboid')
				lambda_1_sq=Pie**2.0*((L_p/a)**(2.0)+(L_p/b)**(2.0)+(L_p/c)**(2.0))
			case ('ellipsoid')
				lambda_1_sq=3.5*((L_p/a)**(2.0)+(L_p/b)**(2.0)+(L_p/c)**(2.0))
			!case ('arbitrary')				
			case DEFAULT
				lambda_1_sq=Pie**2.0
		end select

	    f_L=lambda_1_sq
		k=D_p/L_p**2.0*f_L

	!Medium Bim Range
	else if (Bim > 1E-3) then
	
		select case (sname)
			case ('sphere')
				lambda_1_sq=Pie**2.0
				Bim_crit=A_sphere
			case ('cylinder')
				lambda_1_sq=4.0*alpha_1**(2.0)*(a/b)**(2.0/3.0)+Pie**(2.0)*(b/a)**(4.0/3.0)
				phi_1=a/b
				Bim_crit=(A_cylinder*phi_1**2.0+B_cylinder)/(C_cylinder*phi_1**gamma_cylinder*(1.0+D_cylinder*phi_1))
			case ('cuboid')
				lambda_1_sq=Pie**2.0*((L_p/a)**(2.0)+(L_p/b)**(2.0)+(L_p/c)**(2.0))
				phi_1=a/b
				phi_2=b/c
				Bim_crit=A_cuboid*(phi_1/phi_2)**(gamma_cuboid)*(1.0/phi_1**2.0+1.0+phi_2**2.0)/(1.0/phi_1+1.0+phi_2)
			case ('ellipsoid')
				lambda_1_sq=3.5*((L_p/a)**(2.0)+(L_p/b)**(2.0)+(L_p/c)**(2.0))
				phi_1=a/b
				phi_2=b/c
				Bim_crit=A_ellipsoid*(phi_1/phi_2)**(gamma_ellipsoid)*(1.0/phi_1**2.0+1.0+phi_2**2.0)/((phi_1**p+phi_2**(-p)+(phi_1/phi_2)**p)/3.0)**(1.0/p)
			!case ('arbitrary')				
			case DEFAULT
				lambda_1_sq=Pie**2.0
				Bim_crit=A_sphere
		end select
		
		Bim_crit=Bim_crit/2.0
		f_M=Bim/(Bim+Bim_crit)*lambda_1_sq
		k=D_p/L_p**2.0*f_M
	
	!Small Bim Range
	else
		
		select case (sname)
			case ('sphere')
				e_ratio=3.0/a
			case ('cylinder')
				e_ratio=2.0/a+4.0/b
			case ('cuboid')
				e_ratio=2.0*(1.0/a+1.0/b+1.0/c)
			case ('ellipsoid')
				e_ratio=3.0*(((a*b)**p+(b*c)**p+(c*a)**p)/3.0)**(1.0/p)/(a*b*c)
			!case ('arbitrary')				
			case DEFAULT
				e_ratio=3.0/a
		end select
		
		h_flux=D_w/theta
		f_S=L_p**2.0/D_p*h_flux*e_ratio
		k=D_p/L_p**2.0*f_S
		
	end if

End Subroutine Calc_Toxin_Leach_Rate_SI

Subroutine Calc_Linear_Bottom_Kinetics(dT,mb,decay_rate,ms,mm)
	
	!==========================================================================
	! Generic first order kinetic term:
	!	This subroutine applies the linear decaying kinetic term to each grid
	! cell, the decay is using the plastic integrating step size. If the decay_rate
	! -input is a positive value (default value for MP is referring to A. Charms 
	! et al. 2020), it represents the exponential decaying processes. If the 
	! decay_rate-input is a negative value, it represents the chemical reaction 
	! source term (exponential growth).
	!=========================================================================

	implicit none
	integer, intent(in)		:: mb
	real(sp), intent(out)   :: mm(:)
	real(sp), intent(in)	:: ms(:)
	real(sp), intent(in)	:: dT
	real(sp), intent(in)	:: decay_rate
	
	integer :: i
	
	do i=1,mb
	   mm(i)=ms(i)
	end do
	
	do i=1,mb
	   mm(i)=mm(i)-dT*decay_rate*ms(i)
	end do

End Subroutine Calc_Linear_Bottom_Kinetics

!==========================================================================
! Utility submodules:
!
! AMBI_DENS              - Ambient water density computing base on salinity
!                        and temperature at each nodes.
!==========================================================================
  
SUBROUTINE AMBI_DENS(mb,kb,MYS,MYT,MYRHO)               

	!==========================================================================
	!  Ambient density calculation using the 'dens2' function:
	!==========================================================================

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
   
End Module Mod_Plast