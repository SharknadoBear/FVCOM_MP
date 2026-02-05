!=======================================================================
! FVCOM Plastic Initialization Routine
!   Configuration:  default    
!=======================================================================
 Subroutine Init_Plast 

# if defined (PLASTIC)
  USE MOD_UTILS
  USE MOD_NCTOOLS
  USE MOD_INPUT
  USE ALL_VARS
  USE EQS_OF_STATE
  Use Mod_Prec 
  Use Mod_Plast
  Use Lims, only: m,kbm1,kb
  
  implicit none 
  
  integer :: i,k,iplast
  real(sp), allocatable :: conc_tmp(:,:,:)
  real(sp), allocatable :: tconc_b_tmp(:,:,:)
  real(sp), allocatable :: tconc_f_tmp(:,:,:)
  real(sp), allocatable :: conc_sed_tmp(:,:,:)
  real(sp), allocatable :: conc_sed_s_tmp(:,:,:)
  real(sp), allocatable :: conc_sed_p_tmp(:,:,:)
  real(sp), allocatable :: tconc_b_sed_tmp(:,:,:)
  real(sp), allocatable :: tconc_b_sed_s_tmp(:,:,:)
  real(sp), allocatable :: tconc_b_sed_p_tmp(:,:,:)
  real(sp), allocatable :: tconc_f_sed_tmp(:,:,:)
  real(sp), allocatable :: mass_tmp(:,:)
  real(sp), allocatable :: tmass_tmp(:,:)
  real(sp), allocatable :: sd50_temp(:)  
  real(sp), allocatable :: dens_temp(:)
  real(sp), allocatable :: mudp_temp(:)
  real(sp), allocatable :: poro_temp(:)
  real(sp), allocatable :: Tsed_temp(:)
  real(sp), allocatable :: wDcy_temp(:)
  real(sp), allocatable :: Lsed_temp(:)
#if defined (SEDIMENT)
  real(sp), allocatable :: conc_a(:,:)
  real(sp), allocatable :: conc_d(:,:)
  real(sp), allocatable :: cnew_a(:,:)
  real(sp), allocatable :: cnew_d(:,:)
  real(sp), allocatable :: mass_a(:)
  real(sp), allocatable :: mass_d(:)
  real(sp), allocatable :: mass_a_tmp(:,:)
  real(sp), allocatable :: mass_d_tmp(:,:)
#endif

   TYPE(NCVAR),  POINTER :: VAR
   TYPE(NCDIM),  POINTER :: DIM
   LOGICAL :: FOUND

  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: Init_Plast" 
  
  !--------------------------------------------------
  !Setup static ambient water density 
  !--------------------------------------------------
  
  if (AMBIENT_DENSITY_TYPE==CNSTNT) then
	do i=1,m
		do k=1,kbm1
			rho_water(i,k)=Wrho0
		end do
	end do
  end if
  
  if (AMBIENT_DENSITY_TYPE==STTC) then 
		call LOAD_AMBIENT_DENSITY(rho_water)
  end if

  !--------------------------------------------------
  !Initialize fraction in the active layer
  !--------------------------------------------------

  Do i=1,m
     Do iplast=1,Nplast
        plast(iplast)%frac(i) = 0.0
		If (tox_vector) then
			plast(iplast)%tfrac_b(i) = 0.0
		End if 		
     End Do
  End Do
  
  !--------------------------------------------------
  !Initialize bottom_Mass properties
  !-------------------------------------------------- 
  
  if (INIT_PLAST_MASS_TYPE==STTC) then
	
		allocate(mass_tmp(0:mt,nplast)) ; mass_tmp = 0.0
		If (tox_vector) then
			allocate(tmass_tmp(0:mt,nplast)) ; tmass_tmp = 0.0 
		End if
  
		call LOAD_INIT_PLAST_MASS(mass_tmp)
		If (tox_vector) then
			call LOAD_INIT_TOXIN_MASS(tmass_tmp)
	    End if
  
	!Setup initial sediment mass on the bottom
	Do i=1,m
		Do iplast=1,nplast
			plast(iplast)%mass(i) = mass_tmp(i,iplast)
		End do
	End do
	
	if (tox_vector) then
		!Setup initial bonded toxic mass on the bottom
		Do i=1,m
			Do iplast=1,Nplast
				plast(iplast)%tmass_b(i) = tmass_tmp(i,iplast)
			End do 
		End do
	End if

#if defined (SEDIMENT)
	if (floc_intract) then
		allocate(mass_a_tmp(0:mt,nplast)) ; mass_a_tmp = 0.0
		allocate(mass_d_tmp(0:mt,nplast)) ; mass_d_tmp = 0.0
		
		!Setup initial aggregated/dis-aggregated plastic mass on the bottom
		Do i=1,m
			Do iplast=1,Nplast
				plast(iplast)%mass_a(i) = mass_a_tmp(i,iplast)
				plast(iplast)%mass_d(i) = mass_d_tmp(i,iplast)
			End do 
		End do
	end if
#endif

  End if
  
  if (INIT_PLAST_MASS_TYPE==CNSTNT) then
  
    !Setup initial sediment mass on the bottom
	Do i=1,m
		Do iplast=1,nplast
			plast(iplast)%mass(i) = mass0
		End do
	End do
	
	if (tox_vector) then
		!Setup initial bonded toxic mass on the bottom using TOX_INR ratio
		Do i=1,m
			Do iplast=1,Nplast
				plast(iplast)%tmass_b(i) = plast(iplast)%tox_r * plast(iplast)%mass(i)
			End do 
		End do
	End if

#if defined (SEDIMENT)
	if (floc_intract) then
		Do i=1,m
			Do iplast=1,nplast
				plast(iplast)%mass_a(i) = mass0*Arate_m0
				plast(iplast)%mass_d(i) = mass0*(1.0-Arate_m0)
			End do
		End do
	end if
#endif
	
  End if
  
  !--------------------------------------------------
  !Initialize Plastic Concentrations
  !--------------------------------------------------

  If (INIT_PLAST_CONC_TYPE==CNSTNT) then
  
	  Do k=1,Kb
		Do i=1,m
		  Do iplast=1,Nplast
			plast(iplast)%conc(i,k) = 0.0 
			plast(iplast)%cnew(i,k) = plast(iplast)%conc(i,k) 
		  End Do
		End do
	  End do
  
  End if
  
  If (INIT_PLAST_CONC_TYPE==STTC) then

	allocate(conc_tmp(m,kb,nplast)) ; conc_tmp = 0.0

	call LOAD_INIT_PLAST_CONC(conc_tmp)
	
	Do k=1,Kb
		Do i=1,m
			Do iplast=1,nplast
				plast(iplast)%conc(i,k) = conc_tmp(i,k,iplast)
				plast(iplast)%cnew(i,k) = plast(iplast)%conc(i,k) 
			End do
		End do
	End do
  
  End if

#if defined (SEDIMENT)
  !--------------------------------------------------
  !Split concentrations into aggregated and disaggregated fractions
  !--------------------------------------------------
  if (floc_intract) then
	Do k=1,Kb
		Do i=1,m
			Do iplast=1,nplast
				! Split initial concentration based on Arate_c0
				plast(iplast)%conc_a(i,k) = plast(iplast)%conc(i,k) * Arate_c0
				plast(iplast)%conc_d(i,k) = plast(iplast)%conc(i,k) * (1.0 - Arate_c0)
				plast(iplast)%cnew_a(i,k) = plast(iplast)%conc_a(i,k)
				plast(iplast)%cnew_d(i,k) = plast(iplast)%conc_d(i,k)
			End do
		End do
	End do
  End if
#endif

	If (tox_vector) then

		If (INIT_TOXIN_CONC_TYPE==CNSTNT) then

			Do k=1,Kb
			Do i=1,m
				Do iplast=1,Nplast
				plast(iplast)%tconc_b(i,k) = 0.0 
				plast(iplast)%tcnew_b(i,k) = plast(iplast)%tconc_b(i,k)
				plast(iplast)%tconc_f(i,k) = 0.0 
				plast(iplast)%tcnew_f(i,k) = plast(iplast)%tconc_f(i,k) 				
				End Do
			End do
			End do

		End If

		If (INIT_TOXIN_CONC_TYPE==STTC) then

			allocate(tconc_b_tmp(m,kb,nplast)) ; tconc_b_tmp = 0.0
			allocate(tconc_f_tmp(m,kb,nplast)) ; tconc_f_tmp = 0.0

			call LOAD_INIT_TOXIN_CONC(tconc_b_tmp,tconc_f_tmp)

			Do k=1,Kb
				Do i=1,m
					Do iplast=1,nplast
						plast(iplast)%tconc_b(i,k) = tconc_b_tmp(i,k,iplast)
						plast(iplast)%tcnew_b(i,k) = plast(iplast)%tconc_b(i,k)
						plast(iplast)%tconc_f(i,k) = tconc_f_tmp(i,k,iplast)
						plast(iplast)%tcnew_f(i,k) = plast(iplast)%tconc_f(i,k) 
					End do
				End do
			End do

		End if
	
	End if
  
  !--------------------------------------------------
  !Loading bottom (sediment) properties
  !--------------------------------------------------
 
  if (BOTTOM_FEATURES_TYPE==CNSTNT) then
	Do i=1,m
		bottom(i,isd50) = sd50_0
		bottom(i,idens) = dens0
		bottom(i,imudp) = mudp0
		bottom(i,iporo) = poro0
		if (plast_infil) then
			bottom(i,iTsed) = temp0
			bottom(i,iwDcy) = darv0
			bottom(i,iLsed) = infil_max_depth
		end if
	End do	

    ! Clamp to safe defaults if any constant inputs are non-physical
    Do i=1,m
      if (bottom(i,isd50) <= 0.0_sp) bottom(i,isd50) = sd50_0
      if (bottom(i,iporo) <= 0.0_sp) bottom(i,iporo) = poro0
      if (plast_infil) then
        if (bottom(i,iTsed) <= 0.0_sp) bottom(i,iTsed) = temp0
        if (bottom(i,iLsed) <= 0.0_sp) bottom(i,iLsed) = infil_max_depth
      end if
    End do
  end if

  if (BOTTOM_FEATURES_TYPE==STTC) then

	allocate(sd50_temp(m))
	allocate(dens_temp(m))
	allocate(mudp_temp(m))
	allocate(poro_temp(m))
	if (plast_infil) then
		allocate(Tsed_temp(m))
		allocate(wDcy_temp(m))
		allocate(Lsed_temp(m))
	end if

	call LOAD_BOTTOM_FEATURES(sd50_temp,dens_temp,mudp_temp,poro_temp,Tsed_temp,wDcy_temp,Lsed_temp)
	Do i=1,m
		bottom(i,isd50) = sd50_temp(i)
		bottom(i,idens) = dens_temp(i)
		bottom(i,imudp) = mudp_temp(i)
		bottom(i,iporo) = poro_temp(i)
		if (plast_infil) then
			bottom(i,iTsed) = Tsed_temp(i)
			bottom(i,iwDcy) = wDcy_temp(i)
			bottom(i,iLsed) = Lsed_temp(i)
		end if
	End Do
  End if

  !--------------------------------------------------------------------------
  !Initialization of plastics and toxin concentration within sediment column
  !--------------------------------------------------------------------------
  
  if (plast_infil) then

	if (BOTTOM_FEATURES_TYPE==CNSTNT) then

		do iplast=1,nplast
			do i=1,m
				do k=1,N_sed_layer
					plast(iplast)%conc_sed(i,k)=0.0
					plast(iplast)%conc_sed_s(i,k)=0.0
					plast(iplast)%conc_sed_p(i,k)=0.0
					if (tox_vector) then
						plast(iplast)%tconc_b_sed(i,k)=0.0
						plast(iplast)%tconc_b_sed_s(i,k)=0.0
						plast(iplast)%tconc_b_sed_p(i,k)=0.0
						plast(iplast)%tconc_f_sed(i,k)=0.0
					end if
				end do
			end do
		end do

	end if

	if (BOTTOM_FEATURES_TYPE==STTC) then

		allocate(conc_sed_tmp(m,N_sed_layer,nplast)) ; conc_sed_tmp = 0.0
		allocate(conc_sed_p_tmp(m,N_sed_layer,nplast)) ; conc_sed_p_tmp = 0.0
		allocate(conc_sed_s_tmp(m,N_sed_layer,nplast)) ; conc_sed_s_tmp = 0.0
		If (tox_vector) then
			allocate(tconc_b_sed_tmp(m,N_sed_layer,nplast)) ; tconc_b_sed_tmp = 0.0
			allocate(tconc_b_sed_p_tmp(m,N_sed_layer,nplast)) ; tconc_b_sed_p_tmp = 0.0
			allocate(tconc_b_sed_s_tmp(m,N_sed_layer,nplast)) ; tconc_b_sed_s_tmp = 0.0
			allocate(tconc_f_sed_tmp(m,N_sed_layer,nplast)) ; tconc_f_sed_tmp = 0.0
		End if

		call LOAD_INIT_PLAST_SED_CONC(conc_sed_tmp,conc_sed_s_tmp,conc_sed_p_tmp)
		If (tox_vector) then
			call LOAD_INIT_TOXIN_SED_CONC(tconc_b_sed_tmp,tconc_b_sed_s_tmp,tconc_b_sed_p_tmp,tconc_f_sed_tmp)
		End if

		do iplast=1,nplast
			do i=1,m
				do k=1,N_sed_layer
					plast(iplast)%conc_sed(i,k) = conc_sed_tmp(i,k,iplast)
					plast(iplast)%conc_sed_p(i,k) = conc_sed_p_tmp(i,k,iplast)
					plast(iplast)%conc_sed_s(i,k) = conc_sed_s_tmp(i,k,iplast)
				end do
			end do
		end do

		If (tox_vector) then
			do iplast=1,nplast
				do i=1,m
					do k=1,N_sed_layer
						plast(iplast)%tconc_b_sed(i,k) = tconc_b_sed_tmp(i,k,iplast)
						plast(iplast)%tconc_b_sed_p(i,k) = tconc_b_sed_p_tmp(i,k,iplast)
						plast(iplast)%tconc_b_sed_s(i,k) = tconc_b_sed_s_tmp(i,k,iplast)
						plast(iplast)%tconc_f_sed(i,k) = tconc_f_sed_tmp(i,k,iplast)
					end do
				end do
			end do
		end if

	end if

  end if

  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Init_Plast" 

  Return

# endif

  End Subroutine Init_Plast
!=======================================================================
