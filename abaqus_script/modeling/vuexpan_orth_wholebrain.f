C =============================================================================
C Copyright (c) 2026 Jixin Hou et al.
C All rights reserved.
C
C This code is provided as part of the research software accompanying:
C [Grow2Fold: mapping heterogeneous developmental growth to human brain folding]
C
C Use, modification, and redistribution are permitted under the terms of the
C license provided in the LICENSE file of this repository.
C
C Repository: [Ghttps://github.com/BioDMX-UGA/Grow2Fold]
C =============================================================================

      subroutine vuexpan (
C Read only variables -
     * nblock, nDir, nShr, nExpanType,
     * nElem, nIntPt, nLayer, nSectPt,
     * stepTime, totalTime, dt, cmname,
     * nstatev, nfieldv, nprops, props,
     * tempOld, tempNew, fieldOld, fieldNew,
     * stateOld,
C Write only variables -
     * stateNew, strainThInc, dStrainTherDT )
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * nElem(nblock), props(nprops),
     * tempOld(nblock), tempNew(nblock),
     * fieldOld(nblock,nfieldv),
     * fieldNew(nblock,nfieldv),
     * stateOld(nblock,nstatev),
     * stateNew(nblock,nstatev)
C
      character*80 cmname
C
      if (cmname(1:14) .eq. 'GRAY_REGION_00') then
         call vuexpan_region_00(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_01') then
         call vuexpan_region_01(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_02') then
         call vuexpan_region_02(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_03') then
         call vuexpan_region_03(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_04') then
         call vuexpan_region_04(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_05') then
         call vuexpan_region_05(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_06') then
         call vuexpan_region_06(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_07') then
         call vuexpan_region_07(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_08') then
         call vuexpan_region_08(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_09') then
         call vuexpan_region_09(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_10') then
         call vuexpan_region_10(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_11') then
         call vuexpan_region_11(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_12') then
         call vuexpan_region_12(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_13') then
         call vuexpan_region_13(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_14') then
         call vuexpan_region_14(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_15') then
         call vuexpan_region_15(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_16') then
         call vuexpan_region_16(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_17') then
         call vuexpan_region_17(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_18') then
         call vuexpan_region_18(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_19') then
         call vuexpan_region_19(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_20') then
         call vuexpan_region_20(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_21') then
         call vuexpan_region_21(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_22') then
         call vuexpan_region_22(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_23') then
         call vuexpan_region_23(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_24') then
         call vuexpan_region_24(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_25') then
         call vuexpan_region_25(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_26') then
         call vuexpan_region_26(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_27') then
         call vuexpan_region_27(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_28') then
         call vuexpan_region_28(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_29') then
         call vuexpan_region_29(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_30') then
         call vuexpan_region_30(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_31') then
         call vuexpan_region_31(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_32') then
         call vuexpan_region_32(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_33') then
         call vuexpan_region_33(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_34') then
         call vuexpan_region_34(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_35') then
         call vuexpan_region_35(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_36') then
         call vuexpan_region_36(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_37') then
         call vuexpan_region_37(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_38') then
         call vuexpan_region_38(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_39') then
         call vuexpan_region_39(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_40') then
         call vuexpan_region_40(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:14) .eq. 'GRAY_REGION_41') then
         call vuexpan_region_41(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else if (cmname(1:5) .eq. 'WHITE') then
         call vuexpan_region_white(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
      else
         call xplb_abqerr(-2,'User subroutine VUEXPAN missing!',
     *       intv,zero,' ')
         call xplb_exit
      end if
C
      return
      end

      subroutine vuexpan_region_00(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus00 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus00) then
          print *, 'User subroutine vuexpan_region_00 is being used.'
          PrintStatus00 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (0.1d0)*timeAmp
      alpha2 = (0.1d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_01(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus01 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus01) then
          print *, 'User subroutine vuexpan_region_01 is being used.'
          PrintStatus01 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (1.462d0 - 0.816d0*(timeAmp*totalTime)**1.0d0)*timeAmp
      alpha2 = (1.462d0 - 0.816d0*(timeAmp*totalTime)**1.0d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_02(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus02 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus02) then
          print *, 'User subroutine vuexpan_region_02 is being used.'
          PrintStatus02 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (1.164d0*(timeAmp*totalTime)**1.0d0 + 0.862d0)*timeAmp
      alpha2 = (1.164d0*(timeAmp*totalTime)**1.0d0 + 0.862d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_03(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus03 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus03) then
          print *, 'User subroutine vuexpan_region_03 is being used.'
          PrintStatus03 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (1.122d0*(timeAmp*totalTime)**1.0d0 + 1)*timeAmp
      alpha2 = (1.122d0*(timeAmp*totalTime)**1.0d0 + 1)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_04(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus04 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus04) then
          print *, 'User subroutine vuexpan_region_04 is being used.'
          PrintStatus04 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (0.56d0*(timeAmp*totalTime)**1.0d0 + 0.863d0)*timeAmp
      alpha2 = (0.56d0*(timeAmp*totalTime)**1.0d0 + 0.863d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_05(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus05 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus05) then
          print *, 'User subroutine vuexpan_region_05 is being used.'
          PrintStatus05 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (1.532d0*(timeAmp*totalTime)**1.0d0 + 1.207d0)*timeAmp
      alpha2 = (1.532d0*(timeAmp*totalTime)**1.0d0 + 1.207d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_06(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus06 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus06) then
          print *, 'User subroutine vuexpan_region_06 is being used.'
          PrintStatus06 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (2.0d0*(timeAmp*totalTime)**1.0d0 + 0.717d0)*timeAmp
      alpha2 = (2.0d0*(timeAmp*totalTime)**1.0d0 + 0.717d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_07(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus07 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus07) then
          print *, 'User subroutine vuexpan_region_07 is being used.'
          PrintStatus07 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (0.892d0*(timeAmp*totalTime)**1.0d0 + 1.329d0)*timeAmp
      alpha2 = (0.892d0*(timeAmp*totalTime)**1.0d0 + 1.329d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_08(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus08 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus08) then
          print *, 'User subroutine vuexpan_region_08 is being used.'
          PrintStatus08 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (1.266d0*(timeAmp*totalTime)**1.0d0 + 1.45d0)*timeAmp
      alpha2 = (1.266d0*(timeAmp*totalTime)**1.0d0 + 1.45d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_09(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus09 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus09) then
          print *, 'User subroutine vuexpan_region_09 is being used.'
          PrintStatus09 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (2.702d0*(timeAmp*totalTime)**1.0d0 + 1.207d0)*timeAmp
      alpha2 = (2.702d0*(timeAmp*totalTime)**1.0d0 + 1.207d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_10(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus10 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus10) then
          print *, 'User subroutine vuexpan_region_10 is being used.'
          PrintStatus10 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (2.804d0*(timeAmp*totalTime)**1.0d0 + 1.397d0)*timeAmp
      alpha2 = (2.804d0*(timeAmp*totalTime)**1.0d0 + 1.397d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_11(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus11 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus11) then
          print *, 'User subroutine vuexpan_region_11 is being used.'
          PrintStatus11 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (2.0d0*(timeAmp*totalTime)**1.0d0 + 2.397d0)*timeAmp
      alpha2 = (2.0d0*(timeAmp*totalTime)**1.0d0 + 2.397d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_12(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus12 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus12) then
          print *, 'User subroutine vuexpan_region_12 is being used.'
          PrintStatus12 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (2.73d0*(timeAmp*totalTime)**1.0d0 + 1.185d0)*timeAmp
      alpha2 = (2.73d0*(timeAmp*totalTime)**1.0d0 + 1.185d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_13(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus13 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus13) then
          print *, 'User subroutine vuexpan_region_13 is being used.'
          PrintStatus13 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (2.688d0*(timeAmp*totalTime)**1.0d0 + 1)*timeAmp
      alpha2 = (2.688d0*(timeAmp*totalTime)**1.0d0 + 1)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_14(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus14 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus14) then
          print *, 'User subroutine vuexpan_region_14 is being used.'
          PrintStatus14 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (2.0d0*(timeAmp*totalTime)**1.0d0 + 1.477d0)*timeAmp
      alpha2 = (2.0d0*(timeAmp*totalTime)**1.0d0 + 1.477d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_15(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus15 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus15) then
          print *, 'User subroutine vuexpan_region_15 is being used.'
          PrintStatus15 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (2.0d0*(timeAmp*totalTime)**1.0d0 + 1.759d0)*timeAmp
      alpha2 = (2.0d0*(timeAmp*totalTime)**1.0d0 + 1.759d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_16(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus16 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus16) then
          print *, 'User subroutine vuexpan_region_16 is being used.'
          PrintStatus16 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (2.0d0*(timeAmp*totalTime)**1.0d0 + 1.315d0)*timeAmp
      alpha2 = (2.0d0*(timeAmp*totalTime)**1.0d0 + 1.315d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_17(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus17 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus17) then
          print *, 'User subroutine vuexpan_region_17 is being used.'
          PrintStatus17 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (2.0d0*(timeAmp*totalTime)**1.0d0 + 1.197d0)*timeAmp
      alpha2 = (2.0d0*(timeAmp*totalTime)**1.0d0 + 1.197d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_18(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus18 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus18) then
          print *, 'User subroutine vuexpan_region_18 is being used.'
          PrintStatus18 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (2.67600000000000d0)*timeAmp
      alpha2 = (2.67600000000000d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_19(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus19 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus19) then
          print *, 'User subroutine vuexpan_region_19 is being used.'
          PrintStatus19 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (0.978d0*(timeAmp*totalTime)**1.0d0 + 2.312d0)*timeAmp
      alpha2 = (0.978d0*(timeAmp*totalTime)**1.0d0 + 2.312d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_20(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus20 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus20) then
          print *, 'User subroutine vuexpan_region_20 is being used.'
          PrintStatus20 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (2.656d0*(timeAmp*totalTime)**1.0d0 + 0.701d0)*timeAmp
      alpha2 = (2.656d0*(timeAmp*totalTime)**1.0d0 + 0.701d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_21(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus21 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus21) then
          print *, 'User subroutine vuexpan_region_21 is being used.'
          PrintStatus21 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (0.1d0)*timeAmp
      alpha2 = (0.1d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_22(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus22 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus22) then
          print *, 'User subroutine vuexpan_region_22 is being used.'
          PrintStatus22 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (1.178d0 - 0.756d0*(timeAmp*totalTime)**1.0d0)*timeAmp
      alpha2 = (1.178d0 - 0.756d0*(timeAmp*totalTime)**1.0d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_23(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus23 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus23) then
          print *, 'User subroutine vuexpan_region_23 is being used.'
          PrintStatus23 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (0.968d0*(timeAmp*totalTime)**1.0d0 + 1)*timeAmp
      alpha2 = (0.968d0*(timeAmp*totalTime)**1.0d0 + 1)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_24(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus24 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus24) then
          print *, 'User subroutine vuexpan_region_24 is being used.'
          PrintStatus24 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (1.106d0*(timeAmp*totalTime)**1.0d0 + 1)*timeAmp
      alpha2 = (1.106d0*(timeAmp*totalTime)**1.0d0 + 1)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_25(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus25 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus25) then
          print *, 'User subroutine vuexpan_region_25 is being used.'
          PrintStatus25 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (0.936d0*(timeAmp*totalTime)**1.0d0 + 0.861d0)*timeAmp
      alpha2 = (0.936d0*(timeAmp*totalTime)**1.0d0 + 0.861d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_26(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus26 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus26) then
          print *, 'User subroutine vuexpan_region_26 is being used.'
          PrintStatus26 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (1.586d0*(timeAmp*totalTime)**1.0d0 + 1)*timeAmp
      alpha2 = (1.586d0*(timeAmp*totalTime)**1.0d0 + 1)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_27(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus27 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus27) then
          print *, 'User subroutine vuexpan_region_27 is being used.'
          PrintStatus27 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (2.0d0*(timeAmp*totalTime)**1.0d0 + 0.796d0)*timeAmp
      alpha2 = (2.0d0*(timeAmp*totalTime)**1.0d0 + 0.796d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_28(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus28 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus28) then
          print *, 'User subroutine vuexpan_region_28 is being used.'
          PrintStatus28 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (0.806d0*(timeAmp*totalTime)**1.0d0 + 1.338d0)*timeAmp
      alpha2 = (0.806d0*(timeAmp*totalTime)**1.0d0 + 1.338d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_29(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus29 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus29) then
          print *, 'User subroutine vuexpan_region_29 is being used.'
          PrintStatus29 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (1.796d0*(timeAmp*totalTime)**1.0d0 + 1)*timeAmp
      alpha2 = (1.796d0*(timeAmp*totalTime)**1.0d0 + 1)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_30(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus30 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus30) then
          print *, 'User subroutine vuexpan_region_30 is being used.'
          PrintStatus30 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (3.072d0*(timeAmp*totalTime)**1.0d0 + 1)*timeAmp
      alpha2 = (3.072d0*(timeAmp*totalTime)**1.0d0 + 1)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_31(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus31 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus31) then
          print *, 'User subroutine vuexpan_region_31 is being used.'
          PrintStatus31 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (3.156d0*(timeAmp*totalTime)**1.0d0 + 1)*timeAmp
      alpha2 = (3.156d0*(timeAmp*totalTime)**1.0d0 + 1)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_32(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus32 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus32) then
          print *, 'User subroutine vuexpan_region_32 is being used.'
          PrintStatus32 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (2.0d0*(timeAmp*totalTime)**1.0d0 + 2.071d0)*timeAmp
      alpha2 = (2.0d0*(timeAmp*totalTime)**1.0d0 + 2.071d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_33(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus33 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus33) then
          print *, 'User subroutine vuexpan_region_33 is being used.'
          PrintStatus33 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (3.056d0*(timeAmp*totalTime)**1.0d0 + 1)*timeAmp
      alpha2 = (3.056d0*(timeAmp*totalTime)**1.0d0 + 1)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_34(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus34 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus34) then
          print *, 'User subroutine vuexpan_region_34 is being used.'
          PrintStatus34 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (2.77d0*(timeAmp*totalTime)**1.0d0 + 1)*timeAmp
      alpha2 = (2.77d0*(timeAmp*totalTime)**1.0d0 + 1)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_35(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus35 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus35) then
          print *, 'User subroutine vuexpan_region_35 is being used.'
          PrintStatus35 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (2.0d0*(timeAmp*totalTime)**1.0d0 + 1.442d0)*timeAmp
      alpha2 = (2.0d0*(timeAmp*totalTime)**1.0d0 + 1.442d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_36(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus36 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus36) then
          print *, 'User subroutine vuexpan_region_36 is being used.'
          PrintStatus36 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (2.0d0*(timeAmp*totalTime)**1.0d0 + 1.947d0)*timeAmp
      alpha2 = (2.0d0*(timeAmp*totalTime)**1.0d0 + 1.947d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_37(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus37 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus37) then
          print *, 'User subroutine vuexpan_region_37 is being used.'
          PrintStatus37 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (2.0d0*(timeAmp*totalTime)**1.0d0 + 1.052d0)*timeAmp
      alpha2 = (2.0d0*(timeAmp*totalTime)**1.0d0 + 1.052d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_38(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus38 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus38) then
          print *, 'User subroutine vuexpan_region_38 is being used.'
          PrintStatus38 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (2.0d0*(timeAmp*totalTime)**1.0d0 + 1.29d0)*timeAmp
      alpha2 = (2.0d0*(timeAmp*totalTime)**1.0d0 + 1.29d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_39(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus39 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus39) then
          print *, 'User subroutine vuexpan_region_39 is being used.'
          PrintStatus39 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (2.739d0 - 0.532d0*(timeAmp*totalTime)**1.0d0)*timeAmp
      alpha2 = (2.739d0 - 0.532d0*(timeAmp*totalTime)**1.0d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_40(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus40 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus40) then
          print *, 'User subroutine vuexpan_region_40 is being used.'
          PrintStatus40 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (2.90000000000000d0)*timeAmp
      alpha2 = (2.90000000000000d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_41(nblock, nDir, nShr, nExpanType,
     *       stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus41 = .false.
      logical isStepStart
      real*8 alpha1, alpha2, alpha3
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus41) then
          print *, 'User subroutine vuexpan_region_41 is being used.'
          PrintStatus41 = .true.
      endif
C
      nExpanType = 2
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Orthotropic cortical growth coefficients
      timeAmp = 4.0d0
C
      alpha1 = (2.736d0*(timeAmp*totalTime)**1.0d0 + 0.794d0)*timeAmp
      alpha2 = (2.736d0*(timeAmp*totalTime)**1.0d0 + 0.794d0)*timeAmp
      alpha3 = 0.0d0
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha1*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha2*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha3*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha1*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha2*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha3*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

      subroutine vuexpan_region_white(nblock, nDir, nShr,
     *       nExpanType, stepTime, totalTime, dt,
     *       tempOld, tempNew, strainThInc, dStrainTherDT)
C
      include 'vaba_param.inc'
C
      dimension strainThInc(nblock,nDir+nShr),
     * dStrainTherDT(nblock,nDir+nShr),
     * tempOld(nblock), tempNew(nblock)
C
      logical, save :: PrintStatus_white = .false.
      logical isStepStart
      real*8 alpha
      real*8 epsTime, dTClamp, timeStep1, timeStep2
      integer km
C
      if (.not. PrintStatus_white) then
          print *, 'User subroutine vuexpan_region_white is being used.'
          PrintStatus_white = .true.
      endif
C
      nExpanType = 1
C
C     Step end times
      timeStep1 = 0.08d0
      timeStep2 = 0.25d0
      epsTime = 0.0001d0
      dTClamp = 0.001d0
C
C     Detect start of each internal step
      isStepStart = .false.
      if ((totalTime .gt. timeStep1) .and. (totalTime .lt. timeStep1+epsTime)) then
          isStepStart = .true.
      endif
C
C     Isotropic white matter growth coefficient

      timeAmp = 4.0d0
C      alpha = (0.4d0 / 0.091d0)*timeAmp

      alpha = (3.396d0*(timeAmp*totalTime)**5.0d0 - (10.300d0*(timeAmp*totalTime)**4.0d0 
     * + 7.744d0*(timeAmp*totalTime)**3.0d0 + 1.698d0*(timeAmp*totalTime)**2.0d0 
     * - 4.120d0*(timeAmp*totalTime) + 1.936d0) * 0.4d0 / 0.091d0*timeAmp
C
      do 100 km = 1, nblock
         if (isStepStart) then
      strainThInc(km,1) = alpha*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,2) = alpha*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
      strainThInc(km,3) = alpha*max(-dTClamp, min(dTClamp, tempNew(km)-tempOld(km)))
         else
            strainThInc(km,1) = alpha*(tempNew(km)-tempOld(km))
            strainThInc(km,2) = alpha*(tempNew(km)-tempOld(km))
            strainThInc(km,3) = alpha*(tempNew(km)-tempOld(km))
         endif
 100  continue
C
      return
      end

