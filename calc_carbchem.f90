! **********************************************************************************************************************************
! CALC_carbchem.f90
! CARBONATE CHEMISTRY TEST HARNESS
! **********************************************************************************************************************************


PROGRAM calc_carbchem


  USE gem_carbchem
  USE gem_cmn
  USE gem_util
  IMPLICIT NONE
  SAVE


  ! define variables
  real::loc_var
  REAL::loc_D
  REAL,DIMENSION(n_ocn)::loc_ocn                                       ! ocean tracers 
  REAL,DIMENSION(n_carb)::loc_carb                                     ! carbonate chemistry 
  REAL,DIMENSION(n_carbconst)::loc_carbconst                           ! carbonate chemistry constants
  REAL,DIMENSION(n_carbalk)::loc_carbalk
  REAL,DIMENSION(n_carbisor)::loc_carbisor
  real::loc_DIC_13C,loc_CO2_13C,loc_HCO3_13C,loc_CO3_13C
  real::loc_DIC_14C,loc_CO2_14C,loc_HCO3_14C,loc_CO3_14C
  ! initialize variables
  loc_ocn(:) = 0.0
  loc_carb(:) = 0.0
  loc_carbconst(:) = 0.0
  loc_carbalk(:) = 0.0
  loc_carbisor(:) = 0.0
  ! 
  par_carbchem_pH_iterationmax = 100
  par_carbchem_pH_tolerance = 0.001

  ! *** loop ***
  DO
     ! user inputs
     ! NOTE: ensure units are converted before passing ot dissolution function;
     !       ocean depth: m
     !       dCO3 concentration: mol kg-1
     !       calcite content: normalized (max = 1.0)
     !       organic carbon rain flux: mol cm-2 yr-1
     !       dissolved O2 concentration: mol kg-1
     PRINT*,"Input salinity (o/oo)"
     READ*,loc_var
     loc_ocn(io_S) = loc_var
     PRINT*,"Input [PO4] (umol kg-1)"
     READ*,loc_var
     loc_ocn(io_PO4) = 1.0E-06*loc_var
     PRINT*,"Input [SiO2] (umol kg-1)"
     READ*,loc_var
     loc_ocn(io_SiO2) = 1.0E-06*loc_var
     PRINT*,"Input temperature (degrees C)"
     READ*,loc_var
     loc_ocn(io_T) = loc_var + const_zeroC
     PRINT*,"Input water column depth (m)"
     READ*,loc_var
     loc_D = loc_var
     PRINT*,"Input [ALK] (umol kg-1)"
     READ*,loc_var
     loc_ocn(io_ALK) = 1.0E-06*loc_var
     PRINT*,"Input [DIC] (umol kg-1)"
     READ*,loc_var
     loc_ocn(io_DIC) = 1.0E-06*loc_var
     ! seed carbonate chemistry calculation with a 'reasonable' H+ concentration (pH)
     loc_carb(ic_H) = 10**(-7.8)          
     ! lets have a little H2S and NH4!
     loc_ocn(io_H2S) = 1.0E-9
     loc_ocn(io_NH4) = 1.0E-9
     ! calculate carbonate chemistry constants
     call sub_calc_carbconst( &
          & loc_D,            &
          & loc_ocn(io_T),    &
          & loc_ocn(io_S),    &
          & loc_carbconst(:)  &
          & )
!!$     ! estimate DIC ('tco2') concentration from ALK and dCO3 (both given)
!!$     tco2 = fun_find_DIC_from_dCO3(                              &
!!$          & dco3,loc_DIC_low,loc_DIC_high,                       &
!!$          & loc_PO4,loc_SiO2,alk,loc_B,loc_Ca,loc_SO4,loc_F,     &
!!$          & loc_carbconst(:),loc_carb(:)                         &
!!$          & )
     ! estimate Ca and borate etc concentrations from salinity
     loc_ocn(io_Ca)  = fun_calc_Ca(loc_ocn(io_S))
     loc_ocn(io_B)   = fun_calc_Btot(loc_ocn(io_S))
     loc_ocn(io_SO4) = fun_calc_SO4tot(loc_ocn(io_S))
     loc_ocn(io_F)   = fun_calc_Ftot(loc_ocn(io_S))
     ! calculate carbonate chemistry system
     CALL sub_calc_carb(      &
          & loc_ocn(io_DIC),  &
          & loc_ocn(io_ALK),  &
          & loc_ocn(io_Ca),   &
          & loc_ocn(io_PO4),  &
          & loc_ocn(io_SiO2), &
          & loc_ocn(io_B),    &
          & loc_ocn(io_SO4),  &
          & loc_ocn(io_F),    &
          & loc_ocn(io_H2S),  &
          & loc_ocn(io_NH4),  &
          & loc_carbconst(:), & 
          & loc_carb(:),      &
          & loc_carbalk(:)    & 
          & )
     loc_DIC_13C = fun_calc_isotope_fraction(0.0,const_standards(11))*loc_ocn(io_DIC)
     call sub_calc_carb_r13C( &
          & loc_ocn(io_T), &
          & loc_ocn(io_DIC), &
          & loc_DIC_13C, &
          & loc_carb(:), &
          & loc_carbisor(:) &
          & )
     loc_DIC_14C = fun_calc_isotope_fraction(0.0,const_standards(12))*loc_ocn(io_DIC)
     call sub_calc_carb_r14C( &
          & loc_ocn(io_T), &
          & loc_ocn(io_DIC), &
          & loc_DIC_14C, &
          & loc_carb(:), &
          & loc_carbisor(:) &
          & )
     loc_DIC_13C  = loc_carbisor(ici_DIC_r13C)*loc_ocn(io_DIC)
     loc_CO2_13C  = loc_carbisor(ici_CO2_r13C)*loc_carb(ic_conc_CO2)
     loc_HCO3_13C = loc_carbisor(ici_HCO3_r13C)*loc_carb(ic_conc_HCO3)
     loc_CO3_13C  = loc_carbisor(ici_CO3_r13C)*loc_carb(ic_conc_CO3)
     loc_DIC_14C  = loc_carbisor(ici_DIC_r14C)*loc_ocn(io_DIC)
     loc_CO2_14C  = loc_carbisor(ici_CO2_r14C)*loc_carb(ic_conc_CO2)
     loc_HCO3_14C = loc_carbisor(ici_HCO3_r14C)*loc_carb(ic_conc_HCO3)
     loc_CO3_14C  = loc_carbisor(ici_CO3_r14C)*loc_carb(ic_conc_CO3)
     ! print some info
     PRINT*,'--- INPUTS ---'
     PRINT*,'T          = ',loc_ocn(io_T) - const_zeroC,' degrees C'
     PRINT*,'S          = ',loc_ocn(io_S),' o/oo'
     PRINT*,'D          = ',loc_D,' m'
     PRINT*,'ALK        = ',loc_ocn(io_ALK)*1.0E+06,' umol kg-1'
     PRINT*,'DIC        = ',loc_ocn(io_DIC)*1.0E+06,' umol kg-1'
     PRINT*,'PO4        = ',loc_ocn(io_PO4)*1.0E+06,' umol kg-1'
     PRINT*,'SiO2       = ',loc_ocn(io_SiO2)*1.0E+06,' umol kg-1'
     PRINT*,'--- DERIVED INPUTS ---'
     PRINT*,'Ca         = ',loc_ocn(io_Ca)*1.0E+06,' umol kg-1'
     PRINT*,'B          = ',loc_ocn(io_B)*1.0E+06,' umol kg-1'
     PRINT*,'SO4        = ',loc_ocn(io_SO4)*1.0E+06,' umol kg-1'
     PRINT*,'F          = ',loc_ocn(io_F)*1.0E+06,' umol kg-1'
     PRINT*,'--- BASIC CARBONATE PROPERTIES ---'
     PRINT*,'pH(SWS)    = ',-log10(loc_carb(ic_H)),' '
     PRINT*,'fCO2       = ',loc_carb(ic_fug_CO2)*1.0E+06,' uatm'
     PRINT*,'[CO2]      = ',loc_carb(ic_conc_CO2)*1.0E+06,' umol kg-1'
     PRINT*,'[HCO3-]    = ',loc_carb(ic_conc_HCO3)*1.0E+06,' umol kg-1'
     PRINT*,'[CO3-]     = ',loc_carb(ic_conc_CO3)*1.0E+06,' umol kg-1'
     PRINT*,'ohm_cal    = ',loc_carb(ic_ohm_cal),' '
     PRINT*,'ohm_arg    = ',loc_carb(ic_ohm_arg),' '
     PRINT*,'--- CARBONATE ALKALINITY CONTRIBUTIONS ---'
     PRINT*,'Boron alk  = ',loc_carbalk(ica_H4BO4)*1.0E+06,' umol kg-1'
     PRINT*,'OH alk     = ',loc_carbalk(ica_OH)*1.0E+06,' umol kg-1'
     PRINT*,'pTOT alk   = ',(loc_carbalk(ica_HPO4) + loc_carbalk(ica_PO4) + loc_carbalk(ica_H3PO4))*1.0E+06,' umol kg-1'
     PRINT*,'H3SiO4 alk = ',loc_carbalk(ica_H3SiO4)*1.0E+06,' umol kg-1'
     PRINT*,'HPO4 alk   = ',loc_carbalk(ica_HPO4)*1.0E+06,' umol kg-1'
     PRINT*,'PO4 alk    = ',loc_carbalk(ica_PO4)*1.0E+06,' umol kg-1'
     PRINT*,'NH3 alk    = ',loc_carbalk(ica_NH3)*1.0E+06,' umol kg-1'
     PRINT*,'HS alk     = ',loc_carbalk(ica_HS)*1.0E+06,' umol kg-1'
     PRINT*,'HSO4 alk   = ',loc_carbalk(ica_HSO4)*1.0E+06,' umol kg-1'
     PRINT*,'HF alk     = ',loc_carbalk(ica_HF)*1.0E+06,' umol kg-1'
     PRINT*,'H3PO4 alk  = ',loc_carbalk(ica_H3PO4)*1.0E+06,' umol kg-1'
     PRINT*,'--- DISSOCIATION CONSTANTS ---'
     PRINT*,'pk1        = ',-log10(loc_carbconst(icc_k1)),' '
     PRINT*,'pk2        = ',-log10(loc_carbconst(icc_k2)),' '
     PRINT*,'pkW        = ',-log10(loc_carbconst(icc_kW)),' '
     PRINT*,'pkB        = ',-log10(loc_carbconst(icc_kB)),' '
     PRINT*,'pkSi       = ',-log10(loc_carbconst(icc_kSi)),' '
     PRINT*,'pkHF       = ',-log10(loc_carbconst(icc_kHF)),' '
     PRINT*,'pkHSO4     = ',-log10(loc_carbconst(icc_kHSO4)),' '
     PRINT*,'pkP1       = ',-log10(loc_carbconst(icc_kP1)),' '
     PRINT*,'pkP2       = ',-log10(loc_carbconst(icc_kP2)),' '
     PRINT*,'pkP3       = ',-log10(loc_carbconst(icc_kP3)),' '
     PRINT*,'pkH2S      = ',-log10(loc_carbconst(icc_kH2S)),' '
     PRINT*,'pkNH4      = ',-log10(loc_carbconst(icc_kNH4)),' '
     PRINT*,'pkcal      = ',-log10(loc_carbconst(icc_kcal)),' '
     PRINT*,'pkarg      = ',-log10(loc_carbconst(icc_karg)),' '
     PRINT*,'--- ISOTOPIE PROPERTIES ---'
     PRINT*,'DIC_d13C  = ',fun_calc_isotope_delta(loc_ocn(io_DIC),loc_DIC_13C,const_standards(11),.true.),' o/oo'
     PRINT*,'CO2_d13C  = ',fun_calc_isotope_delta(loc_carb(ic_conc_CO2),loc_CO2_13C,const_standards(11),.true.),' o/oo'
     PRINT*,'HCO3_d13C = ',fun_calc_isotope_delta(loc_carb(ic_conc_HCO3),loc_HCO3_13C,const_standards(11),.true.),' o/oo'
     PRINT*,'CO3_d13C  = ',fun_calc_isotope_delta(loc_carb(ic_conc_CO3),loc_CO3_13C,const_standards(11),.true.),' o/oo'
     PRINT*,'DIC_d14C  = ',fun_calc_isotope_delta(loc_ocn(io_DIC),loc_DIC_14C,const_standards(12),.true.),' o/oo'
     PRINT*,'CO2_d14C  = ',fun_calc_isotope_delta(loc_carb(ic_conc_CO2),loc_CO2_14C,const_standards(12),.true.),' o/oo'
     PRINT*,'HCO3_d14C = ',fun_calc_isotope_delta(loc_carb(ic_conc_HCO3),loc_HCO3_14C,const_standards(12),.true.),' o/oo'
     PRINT*,'CO3_d14C  = ',fun_calc_isotope_delta(loc_carb(ic_conc_CO3),loc_CO3_14C,const_standards(12),.true.),' o/oo'
     PRINT*,'--- END ---'
     PRINT*,' '
  END DO


END PROGRAM calc_carbchem
