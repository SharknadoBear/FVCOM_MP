!=======================================================================
! FVCOM Sediment Initialization Routine
!   Configuration:  default    
!=======================================================================
 Subroutine Init_Sed 

# if defined (SEDIMENT)
  USE MOD_UTILS
  USE MOD_NCTOOLS
  USE MOD_INPUT
  USE ALL_VARS
  USE EQS_OF_STATE
  
# if defined (ORIG_SED)
  Use Mod_Prec 
  Use Mod_Sed
  Use Lims, only: m,kbm1,kb
  
  implicit none 
  
  integer :: i,k,iplast
  real(sp), allocatable :: mass_tmp(:,:)
  real(sp), allocatable :: sd50_temp(:)  
  real(sp), allocatable :: dens_temp(:)
  real(sp), allocatable :: mudp_temp(:)
  real(sp), allocatable :: poro_temp(:)

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
        plast(iplast)%frac(i)=0.0
     End Do
  End Do
  
  !--------------------------------------------------
  !Initialize bottom_Mass properties
  !-------------------------------------------------- 
  
  if (INIT_PLAST_MASS_TYPE==CNSTNT) then
	Do i=1,m
		Do iplast=1,nplast
			plast(iplast)%mass(i) = mass0
		End do
	End do
  end if
  
  if (INIT_PLAST_MASS_TYPE==STTC) then
  
	call LOAD_INIT_PLAST_MASS(mass_tmp) 
  
	Do i=1,m
		Do iplast=1,nplast
			plast(iplast)%mass(i) = mass_tmp(i,iplast)
		End do
	End do

  End if

  !--------------------------------------------------
  !Initialize Plastic Concentrations
  !--------------------------------------------------
  
  Do k=1,Kb
    Do i=1,m
      Do iplast=1,Nplast
        plast(iplast)%conc(i,k) = 0.0 
        plast(iplast)%cnew(i,k) = plast(iplast)%conc(i,k) 
      End Do
    End do
  End do
  
  !--------------------------------------------------
  !Initialize Toxic Concentration
  !--------------------------------------------------
	
   if (tox_vector) then
    Do k=1,Kb
		Do i=1,m
			Do iplast=1,Nplast
				plast(iplast)%tconc_b(i,k) = plast(iplast)%tox_r0*plast(iplast)%conc(i,k)
				plast(iplast)%tconc_f(i,k) = 0.0
				plast(iplast)%tcnew_f(i,k) = plast(iplast)%tconc_f(i,k) 
			End Do
		End do
	End do
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
	End do	
  end if

  if (BOTTOM_FEATURES_TYPE==STTC) then
	call LOAD_BOTTOM_FEATURES(sd50_temp,dens_temp,mudp_temp,poro_temp)
	Do i=1,m
		bottom(i,isd50) = sd50_temp(i)
		bottom(i,idens) = dens_temp(i)
		bottom(i,imudp) = mudp_temp(i)
		bottom(i,iporo) = poro_temp(i)
	End Do
  End if
  
  if(dbg_set(dbg_sbr)) write(ipt,*) "End: Init_Plast" 

  Return

# endif

# endif

  End Subroutine Init_Sed
!=======================================================================
