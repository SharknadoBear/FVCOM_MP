!------------------------------------------------------------------------------!
!----This is the Marine Hydro-Kinetic Module by T.W. in Dec., 2010 for EERE----!
!----Added to V3.2.2 in 2021 for seedling project------------------------------!
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
!---Migrated to FVCOM v431, TPW--------------!

MODULE MOD_MHK
# if defined (MHK)
   USE MOD_TYPES
   USE MOD_PREC
# if defined (MULTIPROCESSOR)
   USE MOD_PAR
# endif
   IMPLICIT NONE
   SAVE

!--File Unit Specifiers -------------------------------------------------------!
   INTEGER :: MHK_MAP

! --> above variables are saved under each time step <--!

   CONTAINS !------------------------------------------------------------------|
            ! ALLOC_VARS_MHK  : Allocate MHK arrays                            | 
            ! IOFILES_MHK  : Open Input Files for MHK Calculation              |
            ! MS_MHK  : Calculate momentum sink for MHK Calculation            |
!==============================================================================|
!   Input Parameters Which Control the Calculation of Air Pressure             |
!==============================================================================|

!==============================================================================!
!   Open & Read Input Files for MHK Calculation                                !
!==============================================================================!
   SUBROUTINE ALLOC_VARS_MHK
!==============================================================================!
   USE ALL_VARS
   IMPLICIT NONE

!==============================================================================!
!  ALLOCATE:                                                                   !
!==============================================================================!
   ALLOCATE(rn_of_turb(0:NT,kbm1))  ; rn_of_turb = ZERO
   ALLOCATE(turb_area(0:NT,kbm1))   ; turb_area  = ZERO
   ALLOCATE(thrust_coef(0:NT,kbm1))   ; thrust_coef  = ZERO
   ALLOCATE(cut_in(0:NT,kbm1))   ; cut_in  = ZERO
   ALLOCATE(cut_out(0:NT,kbm1))  ; cut_out  = ZERO
   ALLOCATE(EMS_X(0:NT,kbm1))  ; EMS_X  = ZERO
   ALLOCATE(EMS_Y(0:NT,kbm1))  ; EMS_Y  = ZERO

   ALLOCATE(MHK_POWER_MS(0:NT,kbm1))  ; MHK_POWER_MS = ZERO
   ALLOCATE(MHK_POWER_MS_DENS(0:NT,kbm1))  ; MHK_POWER_MS_DENS = ZERO
   ALLOCATE(KINETIC_DENS(0:NT,kbm1))  ; KINETIC_DENS = ZERO

   RETURN
   END SUBROUTINE ALLOC_VARS_MHK

!!----------------------------------------------------------------------------!!
   SUBROUTINE IOFILES_MHK
   USE ALL_VARS
   USE MOD_UTILS
   USE CONTROL
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif

   IMPLICIT NONE

   INTEGER :: eid, kid, I,J,K, MHK_MAP
   DOUBLE PRECISION, ALLOCATABLE :: rn_of_turb_gl(:,:), turb_area_gl(:,:)
   DOUBLE PRECISION, ALLOCATABLE :: thrust_coef_gl(:,:), cut_in_gl(:,:), cut_out_gl(:,:)

   MHK_MAP = 135

!-----------------OPEN MHK_MAP FILE (SPATIAL VARYING MHK DEVICE MAP)-----------!

   OPEN(MHK_MAP, file = TRIM(INPUT_DIR)//TRIM(CASENAME)//'_mhk_map.dat')
   
   ALLOCATE(rn_of_turb_gl(0:NGL,kbm1),turb_area_gl(0:NGL,kbm1),thrust_coef_gl(0:NGL,kbm1))
   ALLOCATE(cut_in_gl(0:NGL,kbm1),cut_out_gl(0:NGL,kbm1))

   rn_of_turb_gl      = 0.0
   turb_area_gl       = 0.0
   thrust_coef_gl     = 0.0
   cut_in_gl          = 0.0
   cut_out_gl         = 0.0

   rewind(MHK_MAP)
   READ(MHK_MAP,*)  !skip the title line
   DO I=1, NGL      !read in element-based inputs
      do k=1, kbm1
         READ(MHK_MAP,*) eid, kid, rn_of_turb_gl(i,k), turb_area_gl(i,k),&
                & thrust_coef_gl(i,k), cut_in_gl(i,k), cut_out_gl(i,k)
      end do
   END DO
   CLOSE(MHK_MAP)

!!--Transfer to local domains from global 
!!FOR SERIAL MODE FIRST...
   IF(SERIAL) THEN
     rn_of_turb(1:N,:)  = rn_of_turb_gl(1:N,:)
     turb_area(1:N,:)   = turb_area_gl(1:N,:)

     thrust_coef(1:N,:) = thrust_coef_gl(1:N,:)
     cut_in(1:N,:)  = cut_in_gl(1:N,:)
     cut_out(1:N,:)  = cut_out_gl(1:N,:)
   END IF

!!FOR PARALLEL MODE
#  if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       DO I=1,N
          rn_of_turb(i,:)  = rn_of_turb_gl(EGID(i),:)
          turb_area(i,:)   = turb_area_gl(EGID(i),:)

          thrust_coef(i,:) = thrust_coef_gl(EGID(i),:)
          cut_in(i,:) = cut_in_gl(EGID(i),:)
          cut_out(i,:) = cut_out_gl(EGID(i),:)
       END DO
       
       DO I=1,NHE
          rn_of_turb(i+N,:)  = rn_of_turb_gl(HE_LST(I),:)
          turb_area(i+N,:)   = turb_area_gl(HE_LST(I),:)

          thrust_coef(i+N,:) = thrust_coef_gl(HE_LST(i),:)
          cut_in(i+N,:) = cut_in_gl(HE_LST(i),:)
          cut_out(i+N,:) = cut_out_gl(HE_LST(i),:)
       END DO
     END IF
#  endif

   DEALLOCATE(rn_of_turb_gl,turb_area_gl,thrust_coef_gl,cut_in_gl,cut_out_gl)

   RETURN
   END SUBROUTINE IOFILES_MHK
!==============================================================================!


!==============================================================================|
!   CALCULATES:  momentum sink by MHK structures for internal mode             !
!==============================================================================|

   SUBROUTINE MS_MHK
!==============================================================================!
   USE ALL_VARS
   USE CONTROL
   IMPLICIT NONE
   INTEGER :: i, k
   DOUBLE PRECISION :: rtemp

!!---momentum_sink = 0.5 * Ce * A * u^2-----!
   EMS_X = 0.0
   EMS_Y = 0.0

   DO i=1,N
      do k=1,kbm1
         rtemp = sqrt(u(i,k)**2+v(i,k)**2)

         IF(rtemp < cut_in(i,k) .or. rtemp > cut_out(i,k)) rtemp = 0.0

         EMS_X(i,k) = 0.5*rn_of_turb(i,k)*thrust_coef(i,k)*turb_area(i,k)*rtemp*u(i,k)
         EMS_Y(i,k) = 0.5*rn_of_turb(i,k)*thrust_coef(i,k)*turb_area(i,k)*rtemp*v(i,k)

         MHK_POWER_MS(i,k) = 0.5*rn_of_turb(i,k)*thrust_coef(i,k)*turb_area(i,k)*rtemp**3*(rho(i,k)+1.0)*1000.0

         MHK_POWER_MS_DENS(i,k) = 0.5*thrust_coef(i,k)*rtemp**3*(rho(i,k)+1.0)*1000.0

         KINETIC_DENS(i,k) = 0.5*rtemp**3*(rho(i,k)+1.0)*1000.0
      end do
   END DO

   RETURN
   END SUBROUTINE MS_MHK
!==============================================================================|

# endif

END MODULE MOD_MHK


