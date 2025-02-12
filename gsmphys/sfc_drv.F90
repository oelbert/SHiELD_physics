       module module_sfc_drv
       contains
!>  \file sfc_drv.f
!!  This file contains the NOAH land surface scheme.
!> \defgroup NOAH NOAH Land Surface
!! @{
!!
!!  The Noah LSM (Chen et al., 1996; Koren et al., 1999; Ek et al., 2003) is targeted for moderate complexity and good computational efficiency for numerical weather prediction and climate models. Thus, it omits subgrid surface tiling and uses a single-layer snowpack. The surface energy balance is solved via a Penman-based approximation for latent heat flux. The Noah model includes packages to simulate soil moisture, soil ice, soil temperature, skin temperature, snow depth, snow water equivalent, energy fluxes such as latent heat, sensible heat and ground heat, and water fluxes such as evaporation and total runoff. The Noah surface infiltration scheme follows that of Schaake et al. (1996) for its treatment of the subgrid variability of precipitation and soil moisture.
!!
!!  On 31 May and 14 June 2005, NCEP extensively upgraded the land-surface component of its Global Forecast System (GFS), including its Global Data Assimilation System (GDAS). The Noah LSM upgrade includes an increase from two (10, 190 cm thick) to four soil layers (10, 30, 60, 100 cm thick), addition of frozen soil physics, new formulations for infiltration and runoff (giving more runoff for unsaturated soils), revised physics of the snowpack and its influence on surface heat fluxes and albedo, tuning and adding canopy resistance parameters, allowing spatially varying root depth, revised treatment of ground heat flux and soil thermal conductivity, reformulation for dependence of direct surface evaporation on first layer soil moisture, and improved seasonality of green vegetation cover. The frozen soil physics includes soil heat sinks/sources from freezing/thawing and influences vertical transport of soil moisture, soil thermal conductivity and heat capacity, and surface infiltration. The prognostic states of snowpack depth and liquid soil moisture were added to the already present prognostic states of snowpack water-equivalent (SWE), total soil moisture (liquid plus frozen), soil temperature, canopy water, and skin temperature. SWE divided by the snowpack depth gives the snowpack density. Total soil moisture minus liquid soil moisture gives the frozen soil moisture (Mitchell et al. 2005)
!!
!!  The addition of Noah LSM greatly reduced the two prominent biases in land-surface processes: 1) an early depletion of snowpack; and 2) a high bias in both surface evaporation and precipitation in the warm season in non-arid mid-latitudes. However, a lower tropospheric warm bias as well as increased surface sensible heat flux emerged, particularly over the arid areas during the daytime. Extensive tests attributed this bias mainly to improper treatment of the thermal roughness length. In May 2011, a new thermal roughness length formulation, which assigned a smaller value for the thermal roughness length compared to the momentum roughness length, was implemented. This greatly reduced the warm surface air temperature bias and the cold skin temperature bias over the arid areas during the daytime (Wei et al. 2009; Zheng et al. 2012).
!!
!!  In January 2015, CFS/GLDAS soil moisture climatology at T574 was used for soil moisture nudge to replace the out-of-date coarse resolution bucket soil moisture climatology; a dependence of the ratio of the thermal and momentum roughness on vegetation type was added to address the land-atmosphere coupling strength; a look-up table based on vegetation type was used to replace 1.0 degree momentum roughness length climatology. After this implementation summer warm/dry biases were found over cropland/grassland areas. Some evaporation-related parameters were refined to increase the evaporation to address this issue. The refinement was implemented in May 2016.
!!
!!  In July 2017, new high-resolution MODIS-based snow-free albedo, maximum snow albedo, soil type and vegetation type were used to address the cold biases over the snow area and the blockiness of surface fields due to the coarse resolution data of soil type and vegetation type. The surface layer parameterization scheme was upgraded to modify the roughness-length formulation and introduce a stability parameter constraint in the Monin-Obukhov similarity theory to prevent the land-atmosphere system from decoupling which causes the rapid temperature drop during the sunset (Zheng et al. 2017).
!!
!!  \section diagram Calling Hierarchy Diagram
!!  \section intraphysics Intraphysics Communication
!!
!> \brief Brief description of the subroutine
!!
!!
!! \section arg_table_Noah_run Arguments
!! | local var name | longname                                           | description                        | units   | rank | type    |    kind   | intent | optional |
!! |----------------|----------------------------------------------------|------------------------------------|---------|------|---------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                             | horizontal loop extent, start at 1 | index   |    0 | integer |           | in     | F        |
!!
!!  \section general General Algorithm
!!  \section detailed Detailed Algorithm
!!  @{
! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  usage:                                                               !
!                                                                       !
!      call sfc_drv                                                     !
!  ---  inputs:                                                         !
!          ( im, km, ps, t1, q1, soiltyp, vegtype, sigmaf,              !
!            sfcemis, dlwflx, dswsfc, snet, delt, tg3, cm, ch,          !
!            prsl1, prslki, zf, land, wind,  slopetyp,                  !
!            shdmin, shdmax, snoalb, sfalb, flag_iter, flag_guess,      !
!            lheatstrg, isot, ivegsrc,                                  !
!  ---  in/outs:                                                        !
!            weasd, snwdph, tskin, tprcp, srflag, smc, stc, slc,        !
!            canopy, trans, tsurf, zorl,                                !
!  ---  outputs:                                                        !
!            sncovr1, qsurf, gflux, drain, evap, hflx, ep, runoff,      !
!            cmm, chh, evbs, evcw, sbsno, snowc, stm, snohf,            !
!            smcwlt2, smcref2, wet1 )                                   !
!                                                                       !
!                                                                       !
!  subprogram called:  sflx                                             !
!                                                                       !
!  program history log:                                                 !
!         xxxx  --             created                                  !
!         200x  -- sarah lu    modified                                 !
!    oct  2006  -- h. wei      modified                                 !
!    apr  2009  -- y.-t. hou   modified to include surface emissivity   !
!                     effect on lw radiation. replaced the comfussing   !
!                     slrad (net sw + dlw) with sfc net sw snet=dsw-usw !
!    sep  2009  -- s. moorthi modification to remove rcl and unit change!
!    nov  2011  -- sarah lu    corrected wet1 calculation
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     im       - integer, horiz dimention and num of used pts      1    !
!     km       - integer, vertical soil layer dimension            1    !
!     ps       - real, surface pressure (pa)                       im   !
!     t1       - real, surface layer mean temperature (k)          im   !
!     q1       - real, surface layer mean specific humidity        im   !
!     soiltyp  - integer, soil type (integer index)                im   !
!     vegtype  - integer, vegetation type (integer index)          im   !
!     sigmaf   - real, areal fractional cover of green vegetation  im   !
!     sfcemis  - real, sfc lw emissivity ( fraction )              im   !
!     dlwflx   - real, total sky sfc downward lw flux ( w/m**2 )   im   !
!     dswflx   - real, total sky sfc downward sw flux ( w/m**2 )   im   !
!     snet     - real, total sky sfc netsw flx into ground(w/m**2) im   !
!     delt     - real, time interval (second)                      1    !
!     tg3      - real, deep soil temperature (k)                   im   !
!     cm       - real, surface exchange coeff for momentum (m/s)   im   !
!     ch       - real, surface exchange coeff heat & moisture(m/s) im   !
!     prsl1    - real, sfc layer 1 mean pressure (pa)              im   !
!     prslki   - real,                                             im   !
!     zf       - real, height of bottom layer (m)                  im   !
!     land     - logical, = T if a point with any land             im   !
!     wind     - real, wind speed (m/s)                            im   !
!     slopetyp - integer, class of sfc slope (integer index)       im   !
!     shdmin   - real, min fractional coverage of green veg        im   !
!     shdmax   - real, max fractnl cover of green veg (not used)   im   !
!     snoalb   - real, upper bound on max albedo over deep snow    im   !
!     sfalb    - real, mean sfc diffused sw albedo (fractional)    im   !
!     flag_iter- logical,                                          im   !
!     flag_guess-logical,                                          im   !
!     lheatstrg- logical, flag for canopy heat storage             1    !
!                         parameterization                              !
!     isot     - integer, sfc soil type data source zobler or statsgo   !
!     ivegsrc  - integer, sfc veg type data source umd or igbp          !
!                                                                       !
!  input/outputs:                                                       !
!     weasd    - real, water equivalent accumulated snow depth (mm) im  !
!     snwdph   - real, snow depth (water equiv) over land          im   !
!     tskin    - real, ground surface skin temperature ( k )       im   !
!     tprcp    - real, total precipitation                         im   !
!     srflag   - real, snow/rain flag for precipitation            im   !
!     smc      - real, total soil moisture content (fractional)   im,km !
!     stc      - real, soil temp (k)                              im,km !
!     slc      - real, liquid soil moisture                       im,km !
!     canopy   - real, canopy moisture content (m)                 im   !
!     trans    - real, total plant transpiration (m/s)             im   !
!     tsurf    - real, surface skin temperature (after iteration)  im   !
!     zorl     - real, surface roughness                           im   !
!                                                                       !
!  outputs:                                                             !
!     sncovr1  - real, snow cover over land (fractional)           im   !
!     qsurf    - real, specific humidity at sfc                    im   !
!     gflux    - real, soil heat flux (w/m**2)                     im   !
!     drain    - real, subsurface runoff (mm/s)                    im   !
!     evap     - real, evaperation from latent heat flux           im   !
!     hflx     - real, sensible heat flux                          im   !
!     ep       - real, potential evaporation                       im   !
!     runoff   - real, surface runoff (m/s)                        im   !
!     cmm      - real,                                             im   !
!     chh      - real,                                             im   !
!     evbs     - real, direct soil evaporation (m/s)               im   !
!     evcw     - real, canopy water evaporation (m/s)              im   !
!     sbsno    - real, sublimation/deposit from snopack (m/s)      im   !
!     snowc    - real, fractional snow cover                       im   !
!     stm      - real, total soil column moisture content (m)      im   !
!     snohf    - real, snow/freezing-rain latent heat flux (w/m**2)im   !
!     smcwlt2  - real, dry soil moisture threshold                 im   !
!     smcref2  - real, soil moisture threshold                     im   !
!     wet1     - real, normalized soil wetness                     im   !
!                                                                       !
!  ====================    end of description    =====================  !

!-----------------------------------
      subroutine sfc_drv                                                &
!...................................
!  ---  inputs:
           ( im, km, ps, t1, q1, soiltyp, vegtype, sigmaf,              &
             sfcemis, dlwflx, dswsfc, snet, delt, tg3, cm, ch,          &
             prsl1, prslki, zf, land, wind, slopetyp,                   &
             shdmin, shdmax, snoalb, sfalb, flag_iter, flag_guess,      &
             lheatstrg, isot, ivegsrc,                                  &
             bexppert, xlaipert, vegfpert,pertvegf,                     &  ! sfc perts, mgehne
             !$ser verbatim sfc_iter,&
!  ---  in/outs:
             weasd, snwdph, tskin, tprcp, srflag, smc, stc, slc,        &
             canopy, trans, tsurf, zorl,                                &
!  ---  outputs:
             sncovr1, qsurf, gflux, drain, evap, hflx, ep, runoff,      &
             cmm, chh, evbs, evcw, sbsno, snowc, stm, snohf,            &
             smcwlt2, smcref2, wet1                                     &
           )
!
      !$ser verbatim use mpi
      !$ser verbatim USE m_serialize, ONLY: fs_is_serialization_on
      use machine , only : kind_phys
      use funcphys, only : fpvs
      use physcons, only : grav   => con_g,    cp   => con_cp,          &
                           hvap   => con_hvap, rd   => con_rd,          &
                           eps    => con_eps, epsm1 => con_epsm1,       &
                           rvrdm1 => con_fvirt

      use surface_perturbation, only : ppfbet

      implicit none

!  ---  constant parameters:
      real(kind=kind_phys), parameter :: cpinv   = 1.0/cp
      real(kind=kind_phys), parameter :: hvapi   = 1.0/hvap
      real(kind=kind_phys), parameter :: elocp   = hvap/cp
      real(kind=kind_phys), parameter :: rhoh2o  = 1000.0
      real(kind=kind_phys), parameter :: a2      = 17.2693882
      real(kind=kind_phys), parameter :: a3      = 273.16
      real(kind=kind_phys), parameter :: a4      = 35.86
      real(kind=kind_phys), parameter :: a23m4   = a2*(a3-a4)

      real(kind=kind_phys), save         :: zsoil_noah(4)
      data zsoil_noah / -0.1, -0.4, -1.0, -2.0 /

!  ---  input:
      integer, intent(in) :: im, km, isot, ivegsrc
      !$ser verbatim integer, intent(in) :: sfc_iter
      real (kind=kind_phys), dimension(5), intent(in) :: pertvegf

      integer, dimension(im), intent(in) :: soiltyp, vegtype, slopetyp

      real (kind=kind_phys), dimension(im), intent(in) :: ps,           &
             t1, q1, sigmaf, sfcemis, dlwflx, dswsfc, snet, tg3, cm,    &
             ch, prsl1, prslki, wind, shdmin, shdmax,                   &
             snoalb, sfalb, zf, bexppert, xlaipert, vegfpert

      real (kind=kind_phys),  intent(in) :: delt

      logical, dimension(im), intent(in) :: flag_iter, flag_guess, land

      logical, intent(in) :: lheatstrg

!  ---  in/out:
      real (kind=kind_phys), dimension(im), intent(inout) :: weasd,     &
             snwdph, tskin, tprcp, srflag, canopy, trans, tsurf, zorl

      real (kind=kind_phys), dimension(im,km), intent(inout) ::         &
             smc, stc, slc

!  ---  output:
      real (kind=kind_phys), dimension(im), intent(out) :: sncovr1,     &
             qsurf, gflux, drain, evap, hflx, ep, runoff, cmm, chh,     &
             evbs, evcw, sbsno, snowc, stm, snohf, smcwlt2, smcref2,    &
             wet1

!  ---  locals:
      real (kind=kind_phys), dimension(im) :: rch, rho,                 &
             q0, qs1, theta1,       weasd_old, snwdph_old,              &
             tprcp_old, srflag_old, tskin_old, canopy_old

      !$ser verbatim real (kind=kind_phys), dimension(im) :: can_swdn, can_ch,&
             !$ser verbatim can_q2, can_q2sat, can_dqsdt2, can_sfctmp, can_sfcprs,&
             !$ser verbatim can_sfcems, can_smcwlt, can_smcref, can_rsmin,&
             !$ser verbatim can_rsmax, can_topt, can_rgl, can_hs, can_xlai, can_rc,&
             !$ser verbatim can_pc, can_rcs, can_rct, can_rcq, can_rcsoil, zerobuff_2d,&
             !$ser verbatim nop_etp, nop_prcp, nop_smcmax, nop_smcwlt, nop_smcref,&
             !$ser verbatim nop_smcdry, nop_cmcmax, nop_shdfac, nop_sbeta,&
             !$ser verbatim nop_sfctmp, nop_sfcems, nop_t24, nop_th2, nop_fdown,&
             !$ser verbatim nop_epsca, nop_bexp, nop_pc, nop_rch, nop_rr, nop_cfactr,&
             !$ser verbatim nop_slope, nop_kdt, nop_frzx, nop_psisat, nop_dksat,&
             !$ser verbatim nop_dwsat, nop_zbot, nop_quartz, nop_fxexp, nop_csoil,&
             !$ser verbatim nop_cmc_in, nop_t1_in, nop_tbot_in, nop_beta_in, nop_eta,&
             !$ser verbatim nop_ssoil, nop_runoff1, nop_runoff2, nop_runoff3, nop_edir,&
             !$ser verbatim nop_ec, nop_ett, nop_drip, nop_dew, nop_flx1, nop_flx3, sop_df1,&
             !$ser verbatim sop_flx2, sop_prcp1_in, sop_sncovr_in, sop_sneqv_in,&
             !$ser verbatim sop_sndens_in, sop_snowh_in, sop_snomlt, sop_esnow,&
             !$ser verbatim nop_cmc_out, nop_t1_out, nop_tbot_out, nop_beta_out,&
             !$ser verbatim sop_cmc_out, sop_t1_out, sop_prcp1_out,&
             !$ser verbatim sop_sncovr_out, sop_sneqv_out, sop_sndens_out, sop_snowh_out,&
             !$ser verbatim sop_tbot_out, sop_beta_out, sop_eta, sop_ssoil, sop_runoff1,&
             !$ser verbatim sop_runoff2, sop_runoff3, sop_edir, sop_ec, sop_ett, sop_drip,&
             !$ser verbatim sop_dew, sop_flx1, sop_flx3, sop_ffrozp

      real (kind=kind_phys), dimension(km) :: et, sldpth, stsoil,       &
             smsoil, slsoil

      !$ser verbatim real (kind=kind_phys), dimension(km) :: cn_zsoil, cn_sh2o,&
             !$ser verbatim np_zsoil, np_rtdis, np_stc_in, np_sh2o_in, np_stc_out, np_sh2o_out,&
             !$ser verbatim np_smc, np_et, sp_stc_out, sp_sh2o_out, sp_smc, sp_et

      real (kind=kind_phys), dimension(im,km) :: zsoil, smc_old,        &
             stc_old, slc_old
      !$ser verbatim real (kind=kind_phys), dimension(im,km) :: can_zsoil, can_sh2o, zerobuff_3d,&
             !$ser verbatim nop_zsoil, nop_rtdis, nop_stc_in, nop_sh2o_in, nop_stc_out, nop_sh2o_out,&
             !$ser verbatim nop_smc, nop_et, sop_stc_out, sop_sh2o_out, sop_smc, sop_et

      real (kind=kind_phys) :: alb, albedo, beta, chx, cmx, cmc,        &
             dew, drip, dqsdt2, ec, edir, ett, eta, esnow, etp,         &
             flx1, flx2, flx3, ffrozp, lwdn, pc, prcp, ptu, q2,         &
             q2sat, solnet, rc, rcs, rct, rcq, rcsoil, rsmin,           &
             runoff1, runoff2, runoff3, sfcspd, sfcprs, sfctmp,         &
             sfcems, sheat, shdfac, shdmin1d, shdmax1d, smcwlt,         &
             smcdry, smcref, smcmax, sneqv, snoalb1d, snowh,            &
             snomlt, sncovr, soilw, soilm, ssoil, tsea, th2, tbot,      &
             xlai, zlvl, swdn, tem, z0, bexpp, xlaip, vegfp,            &
             mv,sv,alphav,betav,vegftmp

      !$ser verbatim real (kind=kind_phys) :: cn_ch, cn_rcsoil, cn_smcwlt, cn_smcref,&
             !$ser verbatim cn_rsmin, cn_rsmax, cn_topt, cn_rgl, cn_hs, cn_xlai, cn_rc,&
             !$ser verbatim cn_pc, cn_rcs, cn_rct, cn_rcq,&
             !$ser verbatim np_etp, np_prcp, np_smcmax, np_smcwlt, np_smcref,&
             !$ser verbatim np_smcdry, np_cmcmax, np_shdfac, np_sbeta,&
             !$ser verbatim np_sfctmp, np_sfcems, np_t24, np_th2, np_fdown,&
             !$ser verbatim np_epsca, np_bexp, np_pc, np_rch, np_rr, np_cfactr,&
             !$ser verbatim np_slope, np_kdt, np_frzx, np_psisat, np_dksat,&
             !$ser verbatim np_dwsat, np_zbot, np_quartz, np_fxexp, np_csoil,&
             !$ser verbatim np_cmc_in, np_t1_in, np_tbot_in, np_beta_in, np_eta,&
             !$ser verbatim np_ssoil, np_runoff1, np_runoff2, np_runoff3, np_edir,&
             !$ser verbatim np_ec, np_ett, np_drip, np_dew, np_flx1, np_flx3, sp_df1,&
             !$ser verbatim sp_flx2, sp_prcp1_in, sp_sncovr_in, sp_sneqv_in,&
             !$ser verbatim sp_sndens_in, sp_snowh_in, sp_snomlt, sp_esnow,&
             !$ser verbatim np_cmc_out, np_t1_out, np_tbot_out, np_beta_out,&
             !$ser verbatim sp_cmc_out, sp_t1_out, sp_prcp1_out,&
             !$ser verbatim sp_sncovr_out, sp_sneqv_out, sp_sndens_out, sp_snowh_out,&
             !$ser verbatim sp_tbot_out, sp_beta_out, sp_eta, sp_ssoil, sp_runoff1,&
             !$ser verbatim sp_runoff2, sp_runoff3, sp_edir, sp_ec, sp_ett, sp_drip,&
             !$ser verbatim sp_dew, sp_flx1, sp_flx3, sp_ffrozp

      !$ser verbatim integer :: cn_nroot, np_nroot, np_ice
      integer :: couple, ice, nsoil, nroot, slope, stype, vtype
      integer :: i, k, iflag
      !$ser verbatim integer, dimension(im) :: can_nroot, nop_nroot, nop_ice

      !$ser verbatim logical :: np_mask, np_lheatstrg, sp_snowng, sp_mask
      !$ser verbatim logical, dimension(im) :: nop_mask, nop_lheatstrg, sop_snowng, sop_mask
!
!===> ...  begin here
!
!  --- ...  save land-related prognostic fields for guess run
      !$ser verbatim print *, 'INFO: inside LSM, serialization is ', ser_on

      !$ser verbatim cn_ch = 0.
      !$ser verbatim cn_rcsoil = 0.
      !$ser verbatim cn_smcwlt = 0.
      !$ser verbatim cn_smcref = 0.
      !$ser verbatim cn_rsmin = 0.
      !$ser verbatim cn_rsmax = 0.
      !$ser verbatim cn_topt = 0.
      !$ser verbatim cn_rgl = 0.
      !$ser verbatim cn_hs = 0.
      !$ser verbatim cn_xlai = 0.
      !$ser verbatim cn_rc = 0.
      !$ser verbatim cn_pc = 0.
      !$ser verbatim cn_rcs = 0.
      !$ser verbatim cn_rct = 0.
      !$ser verbatim cn_rcq = 0.
      !$ser verbatim cn_nroot = 0
      !$ser verbatim do k = 1, km
        !$ser verbatim cn_zsoil(k) = 0.
        !$ser verbatim cn_sh2o(k) = 0.
      !$ser verbatim enddo
      do i = 1, im
        !$ser verbatim can_nroot(i) = 0
        !$ser verbatim can_swdn(i) = 0.
        !$ser verbatim can_ch(i) = 0.
        !$ser verbatim can_q2(i) = 0.
        !$ser verbatim can_q2sat(i) = 0.
        !$ser verbatim can_dqsdt2(i) = 0.
        !$ser verbatim can_sfctmp(i) = 0.
        !$ser verbatim can_sfcprs(i) = 0.
        !$ser verbatim can_sfcems(i) = 0.
        !$ser verbatim can_smcwlt(i) = 0.
        !$ser verbatim can_smcref(i) = 0.
        !$ser verbatim can_rsmin(i) = 0.
        !$ser verbatim can_rsmax(i) = 0.
        !$ser verbatim can_topt(i) = 0.
        !$ser verbatim can_rgl(i) = 0.
        !$ser verbatim can_hs(i) = 0.
        !$ser verbatim can_xlai(i) = 0.
        !$ser verbatim can_rc(i) = 0.
        !$ser verbatim can_pc(i) = 0.
        !$ser verbatim can_rcs(i) = 0.
        !$ser verbatim can_rct(i) = 0.
        !$ser verbatim can_rcq(i) = 0.
        !$ser verbatim can_rcsoil(i) = 0.
        !$ser verbatim zerobuff_2d(i) = 0.
        !$ser verbatim do k = 1, km
          !$ser verbatim can_zsoil(i, k) = 0.
          !$ser verbatim can_sh2o(i, k) = 0.
          !$ser verbatim zerobuff_3d(i, k) = 0.
        !$ser verbatim enddo
        !$ser verbatim nop_mask = .false.
        !$ser verbatim sop_mask = .false.

        if (land(i) .and. flag_guess(i)) then
          weasd_old(i)  = weasd(i)
          snwdph_old(i) = snwdph(i)
          tskin_old(i)  = tskin(i)
          canopy_old(i) = canopy(i)
          tprcp_old(i)  = tprcp(i)
          srflag_old(i) = srflag(i)

          do k = 1, km
            smc_old(i,k) = smc(i,k)
            stc_old(i,k) = stc(i,k)
            slc_old(i,k) = slc(i,k)
          enddo
        endif   ! land & flag_guess
      enddo

!  --- ...  initialization block

      do i = 1, im
        if (flag_iter(i) .and. land(i)) then
          ep(i)     = 0.0
          evap (i)  = 0.0
          hflx (i)  = 0.0
          gflux(i)  = 0.0
          drain(i)  = 0.0
          canopy(i) = max(canopy(i), 0.0)

          evbs (i)  = 0.0
          evcw (i)  = 0.0
          trans(i)  = 0.0
          sbsno(i)  = 0.0
          snowc(i)  = 0.0
          snohf(i)  = 0.0
        endif   ! flag_iter & land
      enddo

!  --- ...  initialize variables

      do i = 1, im
        if (flag_iter(i) .and. land(i)) then
          q0(i)   = max(q1(i), 1.e-8)   !* q1=specific humidity at level 1 (kg/kg)
          theta1(i) = t1(i) * prslki(i) !* adiabatic temp at level 1 (k)

          rho(i) = prsl1(i) / (rd*t1(i)*(1.0+rvrdm1*q0(i)))
          qs1(i) = fpvs( t1(i) )        !* qs1=sat. humidity at level 1 (kg/kg)
          qs1(i) = max(eps*qs1(i) / (prsl1(i)+epsm1*qs1(i)), 1.e-8)
          q0 (i) = min(qs1(i), q0(i))
        endif   ! flag_iter & land
      enddo

      do i = 1, im
        if (flag_iter(i) .and. land(i)) then
          do k = 1, km
            zsoil(i,k) = zsoil_noah(k)
          enddo
        endif   ! flag_iter & land
      enddo

      do i = 1, im
        if (flag_iter(i) .and. land(i)) then

!  --- ...  noah: prepare variables to run noah lsm
!   1. configuration information (c):
!      ------------------------------
!    couple  - couple-uncouple flag (=1: coupled, =0: uncoupled)
!    ffrozp  - flag for snow-rain detection (1.=all snow, 0.=all rain, 0-1 mixed)
!    ice     - sea-ice flag (=1: sea-ice, =0: land)
!    dt      - timestep (sec) (dt should not exceed 3600 secs) = delt
!    zlvl    - height (m) above ground of atmospheric forcing variables
!    nsoil   - number of soil layers (at least 2)
!    sldpth  - the thickness of each soil layer (m)

          couple = 1                      ! run noah lsm in 'couple' mode
! use srflag directly to allow fractional rain/snow
!          if     (srflag(i) == 1.0) then  ! snow phase
!            ffrozp = 1.0
!          elseif (srflag(i) == 0.0) then  ! rain phase
!            ffrozp = 0.0
!          endif
          ffrozp = srflag(i)
          ice = 0

          zlvl = zf(i)

          nsoil = km
          sldpth(1) = - zsoil(i,1)
          do k = 2, km
            sldpth(k) = zsoil(i,k-1) - zsoil(i,k)
          enddo

!   2. forcing data (f):
!      -----------------
!    lwdn    - lw dw radiation flux (w/m2)
!    solnet  - net sw radiation flux (dn-up) (w/m2)
!    sfcprs  - pressure at height zlvl above ground (pascals)
!    prcp    - precip rate (kg m-2 s-1)
!    sfctmp  - air temperature (k) at height zlvl above ground
!    th2     - air potential temperature (k) at height zlvl above ground
!    q2      - mixing ratio at height zlvl above ground (kg kg-1)

          lwdn   = dlwflx(i)         !..downward lw flux at sfc in w/m2
          swdn   = dswsfc(i)         !..downward sw flux at sfc in w/m2
          solnet = snet(i)           !..net sw rad flx (dn-up) at sfc in w/m2
          sfcems = sfcemis(i)

          sfcprs = prsl1(i) 
          prcp   = rhoh2o * tprcp(i) / delt
          sfctmp = t1(i)  
          th2    = theta1(i)
          q2     = q0(i)

!   3. other forcing (input) data (i):
!      ------------------------------
!    sfcspd  - wind speed (m s-1) at height zlvl above ground
!    q2sat   - sat mixing ratio at height zlvl above ground (kg kg-1)
!    dqsdt2  - slope of sat specific humidity curve at t=sfctmp (kg kg-1 k-1)

          sfcspd = wind(i)
          q2sat  =  qs1(i)
          dqsdt2 = q2sat * a23m4/(sfctmp-a4)**2

!   4. canopy/soil characteristics (s):
!      --------------------------------
!    vegtyp  - vegetation type (integer index)                       -> vtype
!    soiltyp - soil type (integer index)                             -> stype
!    slopetyp- class of sfc slope (integer index)                    -> slope
!    shdfac  - areal fractional coverage of green vegetation (0.0-1.0)
!    shdmin  - minimum areal fractional coverage of green vegetation -> shdmin1d
!    ptu     - photo thermal unit (plant phenology for annuals/crops)
!    alb     - backround snow-free surface albedo (fraction)
!    snoalb  - upper bound on maximum albedo over deep snow          -> snoalb1d
!    tbot    - bottom soil temperature (local yearly-mean sfc air temp)

          vtype  = vegtype(i)
          stype  = soiltyp(i)
          slope  = slopetyp(i)
          shdfac = sigmaf(i)

!  perturb vegetation fraction that goes into sflx, use the same
!  perturbation strategy as for albedo (percentile matching)
          vegfp  = vegfpert(i)                    ! sfc-perts, mgehne
          if (pertvegf(1) > 0.0) then
                ! compute beta distribution parameters for vegetation fraction
                mv = shdfac
                sv = pertvegf(1)*mv*(1.-mv)
                alphav = mv*mv*(1.0-mv)/(sv*sv)-mv
                betav  = alphav*(1.0-mv)/mv
! compute beta distribution value corresponding
! to the given percentile albPpert to use as new albedo
                call ppfbet(vegfp,alphav,betav,iflag,vegftmp)
                shdfac = vegftmp
          endif
! *** sfc-perts, mgehne

          shdmin1d = shdmin(i)
          shdmax1d = shdmax(i)
          snoalb1d = snoalb(i)

          ptu  = 0.0
          alb  = sfalb(i)
          tbot = tg3(i)

!   5. history (state) variables (h):
!      ------------------------------
!    cmc     - canopy moisture content (m)
!    t1      - ground/canopy/snowpack) effective skin temperature (k)   -> tsea
!    stc(nsoil) - soil temp (k)                                         -> stsoil
!    smc(nsoil) - total soil moisture content (volumetric fraction)     -> smsoil
!    sh2o(nsoil)- unfrozen soil moisture content (volumetric fraction)  -> slsoil
!    snowh   - actual snow depth (m)
!    sneqv   - liquid water-equivalent snow depth (m)
!    albedo  - surface albedo including snow effect (unitless fraction)
!    ch      - surface exchange coefficient for heat and moisture (m s-1) -> chx
!    cm      - surface exchange coefficient for momentum (m s-1)          -> cmx

          cmc  = canopy(i) * 0.001           ! convert from mm to m
          tsea = tsurf(i)                    ! clu_q2m_iter

          do k = 1, km
            stsoil(k) = stc(i,k)
            smsoil(k) = smc(i,k)
            slsoil(k) = slc(i,k)
          enddo

          snowh = snwdph(i) * 0.001         ! convert from mm to m
          sneqv = weasd(i)  * 0.001         ! convert from mm to m
          if (sneqv /= 0.0 .and. snowh == 0.0) then
            snowh = 10.0 * sneqv
          endif

          chx    = ch(i)  * wind(i)              ! compute conductance
          cmx    = cm(i)  * wind(i)
          chh(i) = chx * rho(i)
          cmm(i) = cmx

!  ---- ... outside sflx, roughness uses cm as unit
          z0 = zorl(i)/100.
!  ---- mgehne, sfc-perts
          bexpp  = bexppert(i)                   ! sfc perts, mgehne
          xlaip  = xlaipert(i)                   ! sfc perts, mgehne

!  --- ...  call noah lsm
          call sflx                                                     &
!  ---  inputs:
           ( nsoil, couple, ice, ffrozp, delt, zlvl, sldpth,            &
             swdn, solnet, lwdn, sfcems, sfcprs, sfctmp,                &
             sfcspd, prcp, q2, q2sat, dqsdt2, th2, ivegsrc,             &
             vtype, stype, slope, shdmin1d, alb, snoalb1d,              &
             bexpp, xlaip,                                              & ! sfc-perts, mgehne
             lheatstrg,                                                 &
!  ---  input/outputs:
             tbot, cmc, tsea, stsoil, smsoil, slsoil, sneqv, chx, cmx,  &
             z0,                                                        &
!  ---  outputs:
             nroot, shdfac, snowh, albedo, eta, sheat, ec,              &
             edir, et, ett, esnow, drip, dew, beta, etp, ssoil,         &
             flx1, flx2, flx3, runoff1, runoff2, runoff3,               &
             snomlt, sncovr, rc, pc, rsmin, xlai, rcs, rct, rcq,        &
             !$ser verbatim cn_nroot, cn_ch, cn_zsoil, cn_rcsoil,&
             !$ser verbatim cn_sh2o, cn_smcwlt, cn_smcref, cn_rsmin,&
             !$ser verbatim cn_rsmax, cn_topt, cn_rgl, cn_hs, cn_xlai,&
             !$ser verbatim cn_rc, cn_pc, cn_rcs, cn_rct, cn_rcq,& 
             !$ser verbatim np_mask, sp_mask, np_lheatstrg, sp_snowng, np_nroot,&
             !$ser verbatim np_ice, np_etp, np_prcp, np_smcmax, np_smcwlt, np_smcref,&
             !$ser verbatim np_smcdry, np_cmcmax, np_shdfac, np_sbeta,&
             !$ser verbatim np_sfctmp, np_sfcems, np_t24, np_th2, np_fdown,&
             !$ser verbatim np_epsca, np_bexp, np_pc, np_rch, np_rr, np_cfactr,&
             !$ser verbatim np_slope, np_kdt, np_frzx, np_psisat, np_dksat,&
             !$ser verbatim np_dwsat, np_zbot, np_quartz, np_fxexp, np_csoil,&
             !$ser verbatim np_cmc_in, np_t1_in, np_tbot_in, np_beta_in, np_eta,&
             !$ser verbatim np_ssoil, np_runoff1, np_runoff2, np_runoff3, np_edir,&
             !$ser verbatim np_ec, np_ett, np_drip, np_dew, np_flx1, np_flx3, sp_df1,&
             !$ser verbatim sp_flx2, sp_prcp1_in, sp_sncovr_in, sp_sneqv_in,&
             !$ser verbatim sp_sndens_in, sp_snowh_in, sp_snomlt, sp_esnow,&
             !$ser verbatim np_cmc_out, np_t1_out, np_tbot_out, np_beta_out,&
             !$ser verbatim sp_cmc_out, sp_t1_out, sp_prcp1_out,&
             !$ser verbatim sp_sncovr_out, sp_sneqv_out, sp_sndens_out, sp_snowh_out,&
             !$ser verbatim sp_tbot_out, sp_beta_out, sp_eta, sp_ssoil, sp_runoff1,&
             !$ser verbatim sp_runoff2, sp_runoff3, sp_edir, sp_ec, sp_ett, sp_drip,&
             !$ser verbatim sp_dew, sp_flx1, sp_flx3,&
             !$ser verbatim np_zsoil, np_rtdis, np_stc_in, np_sh2o_in, np_stc_out, np_sh2o_out,&
             !$ser verbatim np_smc, np_et, sp_stc_out, sp_sh2o_out, sp_smc, sp_et, sp_ffrozp,&
             rcsoil, soilw, soilm, smcwlt, smcdry, smcref, smcmax)

             !$ser verbatim do k = 1, km
                !$ser verbatim can_zsoil(i, k) = cn_zsoil(k)
                !$ser verbatim can_sh2o(i, k) = cn_sh2o(k)
                !$ser verbatim nop_zsoil(i, k) = np_zsoil(k)
                !$ser verbatim nop_rtdis(i, k) = np_rtdis(k)
                !$ser verbatim nop_stc_in(i, k) = np_stc_in(k)
                !$ser verbatim nop_sh2o_in(i, k) = np_sh2o_in(k)
                !$ser verbatim nop_stc_out(i, k) = np_stc_out(k)
                !$ser verbatim nop_sh2o_out(i, k) = np_sh2o_out(k)
                !$ser verbatim nop_smc(i, k) = np_smc(k)
                !$ser verbatim nop_et(i, k) = np_et(k)
                !$ser verbatim sop_stc_out(i, k) = sp_stc_out(k)
                !$ser verbatim sop_sh2o_out(i, k) = sp_sh2o_out(k)
                !$ser verbatim sop_smc(i, k) = sp_smc(k)
                !$ser verbatim sop_et(i, k) = sp_et(k)
             !$ser verbatim enddo
             !$ser verbatim can_nroot(i) = cn_nroot
             !$ser verbatim can_rcsoil(i) = cn_rcsoil
             !$ser verbatim can_rc(i) = cn_rc
             !$ser verbatim can_pc(i) = cn_pc
             !$ser verbatim can_rcs(i) = cn_rcs
             !$ser verbatim can_rct(i) = cn_rct
             !$ser verbatim can_rcq(i) = cn_rcq
             !$ser verbatim can_ch(i) = cn_ch
             !$ser verbatim can_smcwlt(i) = cn_smcwlt
             !$ser verbatim can_smcref(i) = cn_smcref
             !$ser verbatim can_rsmin(i) = cn_rsmin
             !$ser verbatim can_rsmax(i) = cn_rsmax
             !$ser verbatim can_topt(i) = cn_topt
             !$ser verbatim can_rgl(i) = cn_rgl
             !$ser verbatim can_hs(i) = cn_hs
             !$ser verbatim can_xlai(i) = cn_xlai

             !$ser verbatim can_swdn = swdn
             !$ser verbatim can_q2 = q2
             !$ser verbatim can_q2sat = q2sat
             !$ser verbatim can_dqsdt2 = dqsdt2
             !$ser verbatim can_sfctmp = sfctmp
             !$ser verbatim can_sfcems = sfcems
             !$ser verbatim can_sfcprs = sfcprs

             !$ser verbatim nop_mask(i) = np_mask
             !$ser verbatim sop_mask(i) = sp_mask
             !$ser verbatim nop_lheatstrg(i) = np_lheatstrg
             !$ser verbatim sop_snowng(i) = sp_snowng
             !$ser verbatim nop_nroot(i) = np_nroot
             !$ser verbatim nop_ice(i) = np_ice
             !$ser verbatim nop_etp(i) = np_etp
             !$ser verbatim nop_prcp(i) = np_prcp
             !$ser verbatim nop_smcmax(i) = np_smcmax
             !$ser verbatim nop_smcwlt(i) = np_smcwlt
             !$ser verbatim nop_smcref(i) = np_smcref
             !$ser verbatim nop_smcdry(i) = np_smcdry
             !$ser verbatim nop_cmcmax(i) = np_cmcmax
             !$ser verbatim nop_shdfac(i) = np_shdfac
             !$ser verbatim nop_sbeta(i) = np_sbeta
             !$ser verbatim nop_sfctmp(i) = np_sfctmp
             !$ser verbatim nop_sfcems(i) = np_sfcems
             !$ser verbatim nop_t24(i) = np_t24
             !$ser verbatim nop_th2(i) = np_th2
             !$ser verbatim nop_fdown(i) = np_fdown
             !$ser verbatim nop_epsca(i) = np_epsca
             !$ser verbatim nop_bexp(i) = np_bexp
             !$ser verbatim nop_pc(i) = np_pc
             !$ser verbatim nop_rch(i) = np_rch
             !$ser verbatim nop_rr(i) = np_rr
             !$ser verbatim nop_cfactr(i) = np_cfactr
             !$ser verbatim nop_slope(i) = np_slope
             !$ser verbatim nop_kdt(i) = np_kdt
             !$ser verbatim nop_frzx(i) = np_frzx
             !$ser verbatim nop_psisat(i) = np_psisat
             !$ser verbatim nop_dksat(i) = np_dksat
             !$ser verbatim nop_dwsat(i) = np_dwsat
             !$ser verbatim nop_zbot(i) = np_zbot
             !$ser verbatim nop_quartz(i) = np_quartz
             !$ser verbatim nop_fxexp(i) = np_fxexp
             !$ser verbatim nop_csoil(i) = np_csoil
             !$ser verbatim nop_cmc_in(i) = np_cmc_in
             !$ser verbatim nop_t1_in(i) = np_t1_in
             !$ser verbatim nop_tbot_in(i) = np_tbot_in
             !$ser verbatim nop_beta_in(i) = np_beta_in
             !$ser verbatim nop_eta(i) = np_eta
             !$ser verbatim nop_ssoil(i) = np_ssoil
             !$ser verbatim nop_runoff1(i) = np_runoff1
             !$ser verbatim nop_runoff2(i) = np_runoff2
             !$ser verbatim nop_runoff3(i) = np_runoff3
             !$ser verbatim nop_edir(i) = np_edir
             !$ser verbatim nop_ec(i) = np_ec
             !$ser verbatim nop_ett(i) = np_ett
             !$ser verbatim nop_drip(i) = np_drip
             !$ser verbatim nop_dew(i) = np_dew
             !$ser verbatim nop_flx1(i) = np_flx1
             !$ser verbatim nop_flx3(i) = np_flx3
             !$ser verbatim sop_df1(i) = sp_df1
             !$ser verbatim sop_flx2(i) = sp_flx2
             !$ser verbatim sop_prcp1_in(i) = sp_prcp1_in
             !$ser verbatim sop_sncovr_in(i) = sp_sncovr_in
             !$ser verbatim sop_sneqv_in(i) = sp_sneqv_in
             !$ser verbatim sop_sndens_in(i) = sp_sndens_in
             !$ser verbatim sop_snowh_in(i) = sp_snowh_in
             !$ser verbatim nop_cmc_out(i) = np_cmc_out
             !$ser verbatim nop_t1_out(i) = np_t1_out
             !$ser verbatim nop_tbot_out(i) = np_tbot_out
             !$ser verbatim nop_beta_out(i) = np_beta_out
             
             !$ser verbatim sop_cmc_out(i) = sp_cmc_out
             !$ser verbatim sop_t1_out(i) = sp_t1_out
             !$ser verbatim sop_prcp1_out(i) = sp_prcp1_out
             !$ser verbatim sop_sncovr_out(i) = sp_sncovr_out
             !$ser verbatim sop_sneqv_out(i) = sp_sneqv_out
             !$ser verbatim sop_sndens_out(i) = sp_sndens_out
             !$ser verbatim sop_snowh_out(i) = sp_snowh_out
             !$ser verbatim sop_tbot_out(i) = sp_tbot_out
             !$ser verbatim sop_beta_out(i) = sp_beta_out
             !$ser verbatim sop_eta(i) = sp_eta
             !$ser verbatim sop_ssoil(i) = sp_ssoil
             !$ser verbatim sop_runoff1(i) = sp_runoff1
             !$ser verbatim sop_runoff2(i) = sp_runoff2
             !$ser verbatim sop_runoff3(i) = sp_runoff3
             !$ser verbatim sop_edir(i) = sp_edir
             !$ser verbatim sop_ec(i) = sp_ec
             !$ser verbatim sop_ett(i) = sp_ett
             !$ser verbatim sop_drip(i) = sp_drip
             !$ser verbatim sop_dew(i) = sp_dew
             !$ser verbatim sop_flx1(i) = sp_flx1
             !$ser verbatim sop_flx3(i) = sp_flx3
             !$ser verbatim sop_snomlt(i) = sp_snomlt
             !$ser verbatim sop_esnow(i) = sp_esnow
             !$ser verbatim sop_ffrozp(i) = sp_ffrozp

!  --- ...  noah: prepare variables for return to parent mode
!   6. output (o):
!      -----------
!    eta     - actual latent heat flux (w m-2: positive, if upward from sfc)
!    sheat   - sensible heat flux (w m-2: positive, if upward from sfc)
!    beta    - ratio of actual/potential evap (dimensionless)
!    etp     - potential evaporation (w m-2)
!    ssoil   - soil heat flux (w m-2: negative if downward from surface)
!    runoff1 - surface runoff (m s-1), not infiltrating the surface
!    runoff2 - subsurface runoff (m s-1), drainage out bottom

          evap(i)  = eta
          hflx(i)  = sheat
          gflux(i) = ssoil

          evbs(i)  = edir
          evcw(i)  = ec
          trans(i) = ett
          sbsno(i) = esnow
          snowc(i) = sncovr
          stm(i)   = soilm * 1000.0 ! unit conversion (from m to kg m-2)
          snohf(i) = flx1 + flx2 + flx3

          smcwlt2(i) = smcwlt
          smcref2(i) = smcref

          ep(i)      = etp
          tsurf(i)   = tsea

          do k = 1, km
            stc(i,k) = stsoil(k) 
            smc(i,k) = smsoil(k)
            slc(i,k) = slsoil(k)
          enddo
          wet1(i) = smsoil(1) / smcmax !Sarah Lu added 09/09/2010 (for GOCART)

!  --- ...  unit conversion (from m s-1 to mm s-1 and kg m-2 s-1)
          runoff(i)  = runoff1 * 1000.0
          drain (i)  = runoff2 * 1000.0

!  --- ...  unit conversion (from m to mm)
          canopy(i)  = cmc   * 1000.0
          snwdph(i)  = snowh * 1000.0
          weasd(i)   = sneqv * 1000.0
          sncovr1(i) = sncovr
!  ---- ... outside sflx, roughness uses cm as unit (update after snow's
!  effect)
          zorl(i) = z0*100.

!  --- ...  do not return the following output fields to parent model
!    ec      - canopy water evaporation (m s-1)
!    edir    - direct soil evaporation (m s-1)
!    et(nsoil)-plant transpiration from a particular root layer (m s-1)
!    ett     - total plant transpiration (m s-1)
!    esnow   - sublimation from (or deposition to if <0) snowpack (m s-1)
!    drip    - through-fall of precip and/or dew in excess of canopy
!              water-holding capacity (m)
!    dew     - dewfall (or frostfall for t<273.15) (m)
!    beta    - ratio of actual/potential evap (dimensionless)
!    flx1    - precip-snow sfc (w m-2)
!    flx2    - freezing rain latent heat flux (w m-2)
!    flx3    - phase-change heat flux from snowmelt (w m-2)
!    snomlt  - snow melt (m) (water equivalent)
!    sncovr  - fractional snow cover (unitless fraction, 0-1)
!    runoff3 - numerical trunctation in excess of porosity (smcmax)
!              for a given soil layer at the end of a time step
!    rc      - canopy resistance (s m-1)
!    pc      - plant coefficient (unitless fraction, 0-1) where pc*etp
!              = actual transp
!    xlai    - leaf area index (dimensionless)
!    rsmin   - minimum canopy resistance (s m-1)
!    rcs     - incoming solar rc factor (dimensionless)
!    rct     - air temperature rc factor (dimensionless)
!    rcq     - atmos vapor pressure deficit rc factor (dimensionless)
!    rcsoil  - soil moisture rc factor (dimensionless)
!    soilw   - available soil moisture in root zone (unitless fraction
!              between smcwlt and smcmax)
!    soilm   - total soil column moisture content (frozen+unfrozen) (m)
!    smcwlt  - wilting point (volumetric)
!    smcdry  - dry soil moisture threshold where direct evap frm top
!              layer ends (volumetric)
!    smcref  - soil moisture threshold where transpiration begins to
!              stress (volumetric)
!    smcmax  - porosity, i.e. saturated value of soil moisture
!              (volumetric)
!    nroot   - number of root layers, a function of veg type, determined
!              in subroutine redprm.

        endif   ! flag_iter and flag
      enddo   ! end do_i_loop

      !$ser verbatim if (sfc_iter == 1) then
          !$ser savepoint Canres1-In
        !$ser verbatim else
          !$ser savepoint Canres2-In
        !$ser verbatim end if
      !$ser data nsoil=nsoil nroot=can_nroot swdn=can_swdn ch=can_ch q2=can_q2
      !$ser data q2sat=can_q2sat dqsdt2=can_dqsdt2 sfctmp=can_sfctmp sfcprs=can_sfcprs
      !$ser data sfcems=can_sfcems sh2o=can_sh2o smcwlt=can_smcwlt smcref=can_smcref
      !$ser data zsoil=can_zsoil rsmin=can_rsmin rsmax=can_rsmax topt=can_topt
      !$ser data rgl=can_rgl hs=can_hs xlai=can_xlai rc=zerobuff_2d pc=zerobuff_2d
      !$ser data rcs=zerobuff_2d rct=zerobuff_2d rcq=zerobuff_2d rcsoil=zerobuff_2d
      !$ser data lsm_mask=land can_shdfac=nop_shdfac

      !$ser verbatim if (sfc_iter == 1) then
        !$ser savepoint Nopack1-In
      !$ser verbatim else
        !$ser savepoint Nopack2-In
      !$ser verbatim end if
      !$ser data nsoil=nsoil nopac_mask=nop_mask lheatstrg=nop_lheatstrg nroot=nop_nroot
      !$ser data ice=nop_ice etp=nop_etp nop_prcp=nop_prcp smcmax=nop_smcmax smcwlt=nop_smcwlt smcref=nop_smcref
      !$ser data smcdry=nop_smcdry dt=delt cmcmax=nop_cmcmax nop_shdfac=nop_shdfac sbeta=nop_sbeta
      !$ser data sfctmp=nop_sfctmp sfcems=nop_sfcems t24=nop_t24 th2=nop_th2 fdown=nop_fdown
      !$ser data epsca=nop_epsca bexp=nop_bexp pc=nop_pc rch=nop_rch rr=nop_rr cfactr=nop_cfactr
      !$ser data slope=nop_slope kdt=nop_kdt frzx=nop_frzx psisat=nop_psisat dksat=nop_dksat
      !$ser data dwsat=nop_dwsat zbot=nop_zbot quartz=nop_quartz fxexp=nop_fxexp csoil=nop_csoil
      !$ser data cmc=nop_cmc_in nop_t1=nop_t1_in tbot=nop_tbot_in beta=nop_beta_in
      !$ser data ssoil=zerobuff_2d runoff1=zerobuff_2d runoff2=zerobuff_2d runoff3=zerobuff_2d edir=zerobuff_2d
      !$ser data ec=zerobuff_2d ett=zerobuff_2d drip=zerobuff_2d dew=zerobuff_2d flx1=zerobuff_2d flx3=zerobuff_2d
      !$ser data eta=zerobuff_2d
      !$ser data zsoil=nop_zsoil rtdis=nop_rtdis stc=nop_stc_in sh2o=nop_sh2o_in
      !$ser data smc=zerobuff_3d et=zerobuff_3d vegtype=vegtype

      !$ser verbatim if (sfc_iter == 1) then
        !$ser savepoint Snopack1-In
      !$ser verbatim else
        !$ser savepoint Snopack2-In
      !$ser verbatim end if
      !$ser data nsoil=nsoil snopac_mask=sop_mask lheatstrg=nop_lheatstrg snowng=sop_snowng nroot=nop_nroot
      !$ser data ice=nop_ice etp=nop_etp nop_prcp=nop_prcp smcmax=nop_smcmax smcwlt=nop_smcwlt smcref=nop_smcref
      !$ser data smcdry=nop_smcdry cmcmax=nop_cmcmax dt=delt df1=sop_df1 sop_shdfac=nop_shdfac vegtype=vegtype
      !$ser data sfctmp=nop_sfctmp sfcems=nop_sfcems t24=nop_t24 th2=nop_th2 fdown=nop_fdown
      !$ser data epsca=nop_epsca bexp=nop_bexp pc=nop_pc rch=nop_rch rr=nop_rr cfactr=nop_cfactr
      !$ser data slope=nop_slope kdt=nop_kdt frzx=nop_frzx psisat=nop_psisat dksat=nop_dksat
      !$ser data dwsat=nop_dwsat zbot=nop_zbot quartz=nop_quartz fxexp=nop_fxexp csoil=nop_csoil
      !$ser data cmc=nop_cmc_in nop_t1=nop_t1_in tbot=nop_tbot_in beta=nop_beta_in flx2=sop_flx2
      !$ser data ssoil=zerobuff_2d runoff1=zerobuff_2d runoff2=zerobuff_2d runoff3=zerobuff_2d edir=zerobuff_2d
      !$ser data ec=zerobuff_2d ett=zerobuff_2d drip=zerobuff_2d dew=zerobuff_2d flx1=zerobuff_2d flx3=zerobuff_2d
      !$ser data eta=zerobuff_2d prcp1=sop_prcp1_in sncovr=sop_sncovr_in sneqv=sop_sneqv_in sndens=sop_sndens_in
      !$ser data snowh=sop_snowh_in zsoil=nop_zsoil rtdis=nop_rtdis stc=nop_stc_in sh2o=nop_sh2o_in
      !$ser data smc=zerobuff_3d et=zerobuff_3d snomlt=zerobuff_2d esnow=zerobuff_2d ffrozp=sop_ffrozp

      !$ser verbatim if (sfc_iter == 1) then
        !$ser savepoint Canres1-Out
      !$ser verbatim else
        !$ser savepoint Canres2-Out
      !$ser verbatim end if
      !$ser data rc=can_rc pc=can_pc rcs=can_rcs rct=can_rct rcq=can_rcq rcsoil=can_rcsoil

      !$ser verbatim if (sfc_iter == 1) then
        !$ser savepoint Nopack1-Out
      !$ser verbatim else
        !$ser savepoint Nopack2-Out
      !$ser verbatim end if
      !$ser data cmc=nop_cmc_out nop_t1=nop_t1_out tbot=nop_tbot_out beta=nop_beta_out et=nop_et
      !$ser data stc=nop_stc_out sh2o=nop_sh2o_out eta=nop_eta smc=nop_smc ssoil=nop_ssoil
      !$ser data runoff1=nop_runoff1 runoff2=nop_runoff2 runoff3=nop_runoff3 edir=nop_edir
      !$ser data ec=nop_ec ett=nop_ett drip=nop_drip dew=nop_dew flx1=nop_flx1 flx3=nop_flx3

      !$ser verbatim if (sfc_iter == 1) then
        !$ser savepoint Snopack1-Out
      !$ser verbatim else
        !$ser savepoint Snopack2-Out
      !$ser verbatim end if
      !$ser data prcp1=sop_prcp1_out cmc=sop_cmc_out sop_t1=sop_t1_out stc=sop_stc_out
      !$ser data sncovr=sop_sncovr_out sneqv=sop_sneqv_out sndens=sop_sndens_out
      !$ser data snowh=sop_snowh_out sh2o=sop_sh2o_out tbot=sop_tbot_out beta=sop_beta_out
      !$ser data smc=sop_smc ssoil=sop_ssoil runoff1=sop_runoff1 runoff2=sop_runoff2
      !$ser data runoff3=sop_runoff3 edir=sop_edir ec=sop_ec et=sop_et ett=sop_ett
      !$ser data snomlt=sop_snomlt drip=sop_drip dew=sop_dew flx1=sop_flx1 flx3=sop_flx3
      !$ser data esnow=sop_esnow

      !  ---  inputs:
!          ( nsoil, nroot, swdn, ch, q2, q2sat, dqsdt2, sfctmp,         &
!            sfcprs, sfcems, sh2o, smcwlt, smcref, zsoil, rsmin,        &
!            rsmax, topt, rgl, hs, xlai,                                &
!  ---  outputs:
!            rc, pc, rcs, rct, rcq, rcsoil                              &
!          )
!   --- ...  compute qsurf (specific humidity at sfc)

      do i = 1, im
        if (flag_iter(i) .and. land(i)) then
          rch(i)   = rho(i) * cp * ch(i) * wind(i)
          qsurf(i) = q1(i)  + evap(i) / (elocp * rch(i))
        endif   ! flag_iter & flag
      enddo

      do i = 1, im
        if (flag_iter(i) .and. land(i)) then
          tem     = 1.0 / rho(i)
          hflx(i) = hflx(i) * tem * cpinv
          evap(i) = evap(i) * tem * hvapi
        endif   ! flag_iter & flag
      enddo

!  --- ...  restore land-related prognostic fields for guess run

      do i = 1, im
        if (land(i)) then
          if (flag_guess(i)) then
            weasd(i)  = weasd_old(i)
            snwdph(i) = snwdph_old(i)
            tskin(i)  = tskin_old(i)
            canopy(i) = canopy_old(i)
            tprcp(i)  = tprcp_old(i)
            srflag(i) = srflag_old(i)

            do k = 1, km
              smc(i,k) = smc_old(i,k)
              stc(i,k) = stc_old(i,k)
              slc(i,k) = slc_old(i,k)
            enddo
          else    ! flag_guess = F
            tskin(i) = tsurf(i)
          endif   ! flag_guess
        endif     ! flag

      enddo
!
      return
!...................................
      end subroutine sfc_drv
!-----------------------------------
!> @}
!> @}
       end module module_sfc_drv
