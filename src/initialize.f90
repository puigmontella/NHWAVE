










!------------------------------------------------------------------------------------------
!
!
!    initialize.F
!
!-----------------------------------------------------------------------------------------
!
!
!    This file is part of NHWAVE.
!
!    This file contains the subroutines:
!
!        (1) wall_time_secs
!        (2) read_input
!        index
!        allocate_variables
!        generate_grid
!        read_bathymetry
!        initial
!
!------------------------------------------------------------------------------------------
!
!   BSD 2-Clause License
!
!   Copyright (c) 2019, NHWAVE Development Group
!   All rights reserved.
!
!   Redistribution and use in source and binary forms, with or without
!   modification, are permitted provided that the following conditions are met:
!
!   * Redistributions of source code must retain the above copyright notice, this
!     list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above copyright notice,
!     this list of conditions and the following disclaimer in the documentation
!     and/or other materials provided with the distribution.
!
!   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
!   FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!   DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!   CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!   OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!----------------------------------------------------------------------------------------
!
!----------------------------------------------------------------------------------------
!
!   (1) wall_time_secs
!
!   Calculate current wall time
!
!   Called by: main
!
!   Gangfeng Ma, 09/12/2011
!-----------------------------------------------------------------------------------------
!
    subroutine wall_time_secs(tcurrent)
!
    use global, only: SP
    implicit none
    integer, dimension(8) :: walltime
    real(SP), intent(out) :: tcurrent
    real(SP) :: msecs,secs,mins,hrs,days,months,mscale,years

    call date_and_time(VALUES=walltime)

    msecs = real(walltime(8))
    secs = real(walltime(7))
    mins = real(walltime(6))
    hrs = real(walltime(5))
    days = real(walltime(3))
    months = real(walltime(2))
    years = real(walltime(1))

    if((months.eq.1).or.(months.eq.3).or.(months.eq.5).or.  &
          (months.eq.7).or.(months.eq.8).or.(months.eq.10).or.  &                                                                                   
          (months.eq.12)) then
      mscale = 31.0
    elseif((months.eq.4).or.(months.eq.6).or.  &
          (months.eq.9).or.(months.eq.11)) then
      mscale = 30.0
    elseif(years.eq.4*int(years/4)) then
      mscale = 29.0
    else
      mscale = 28.0
    endif

    tcurrent = months*mscale*24.0*60.0*60.0+days*24.0*60.0*60.0+  &
         hrs*60.0*60.0+60.0*mins+secs+msecs/1000.0

    return
!
    end subroutine wall_time_secs
!
!---------------------------------------------------------------------------------------
!
     subroutine read_input
!
!----------------------------------------------------------------------------------------
!
!    (2) read_input
!
!    This subroutine is used to read input.txt
!
!    Called by: main
!
!    Gangfeng Ma, 20/12/2010
!
!----------------------------------------------------------------------------------------
!
     use global
     use input_util


     implicit none
     character(len=80) :: FILE_NAME
     real(SP) :: Segma,Celerity,Wave_Length,Wave_Number,  &
                 Fk,Fkdif,Theta_Calc,Wnumy,tmp,tmp1,DFreq, &
                 Freq_Peak,gam,sa,sb,SumInt,A_Jon
     integer :: line,ierr,Iter,i,j,nw,n,nn
     !Added by M.Derakhti for FOCUSED wavemaker
     real(SP) :: Xltmp, Ctmp, xxk,segmasqr
     real(SP), dimension(12000) :: period_i 

     ! log and error file
     open(3,file='log.txt')

     ! read from input.txt
     FILE_NAME='input.txt'
!
!    title
!
     CALL GET_STRING_VAL(TITLE,FILE_NAME,'TITLE',line,ierr)

     IF(ierr==1)THEN

     if(myid.eq.0) write(3,*) 'No TITLE in ', FILE_NAME, 'use default'

     TITLE='---TEST RUN---'
     ENDIF

     if(myid.eq.0) WRITE(3,*)'---- LOG FILE ---'
     if(myid.eq.0) WRITE(3,*)TITLE
     if(myid.eq.0) WRITE(3,*)'--------------input start --------------'
!
!    dimension     
!                                        
     CALL GET_INTEGER_VAL(Mglob,FILE_NAME,'Mglob',line)
     CALL GET_INTEGER_VAL(Nglob,FILE_NAME,'Nglob',line)
     CALL GET_INTEGER_VAL(Kglob,FILE_NAME,'Kglob',line)
     if(myid.eq.0) WRITE(3,'(A7,I5)')'Mglob= ',Mglob
     if(myid.eq.0) WRITE(3,'(A7,I5)')'Nglob= ',Nglob
     if(myid.eq.0) WRITE(3,'(A7,I5)')'Kglob= ',Kglob

     ! processor number
     CALL GET_INTEGER_VAL(PX,FILE_NAME,'PX',line)
     CALL GET_INTEGER_VAL(PY,FILE_NAME,'PY',line)
     if(myid.eq.0) WRITE(3,'(A4,I5)')'PX= ',PX
     if(myid.eq.0) WRITE(3,'(A4,I5)')'PY= ',PY
     if(PX*PY.ne.NumP) then
       if(myid.eq.0) WRITE(3,'(A6,I5)') 'NumP= ',NumP
       stop
     endif

     ! grid sizes
     CALL GET_Float_VAL(dx,FILE_NAME,'DX',line)
     CALL GET_Float_VAL(dy,FILE_NAME,'DY',line)
     if(myid.eq.0) WRITE(3,'(A4,F8.4)')'DX= ',dx
     if(myid.eq.0) WRITE(3,'(A4,F8.4)')'DY= ',dy

     ! vertical grid option
     call GET_INTEGER_VAL(Ivgrd,FILE_NAME,'IVGRD',line)
     CALL GET_Float_VAL(Grd_R,FILE_NAME,'GRD_R',line)
     if(myid.eq.0) WRITE(3,'(A7,I3)')'Ivgrd= ',Ivgrd
     if(myid.eq.0) WRITE(3,'(A7,f5.2)')'Grd_R= ',Grd_R

     ! time step
     CALL GET_Float_VAL(dt_ini,FILE_NAME,'DT_INI',line)
     CALL GET_Float_VAL(dt_min,FILE_NAME,'DT_MIN',line)
     CALL GET_Float_VAL(dt_max,FILE_NAME,'DT_MAX',line)
     if(myid.eq.0) WRITE(3,'(A8,F8.4)')'DT_INI= ',dt_ini
     if(myid.eq.0) WRITE(3,'(A8,F8.4)')'DT_MIN= ',dt_min
     if(myid.eq.0) WRITE(3,'(A8,F8.4)')'DT_MAX= ',dt_max

! result folder                                     
     CALL GET_STRING_VAL(RESULT_FOLDER,FILE_NAME,'RESULT_FOLDER',line,ierr)
     if(myid.eq.0) WRITE(3,'(A15,A50)')'RESULT_FOLDER= ', RESULT_FOLDER

     ! simulation steps and time
     call GET_INTEGER_VAL(SIM_STEPS,FILE_NAME,'SIM_STEPS',line)
     CALL GET_Float_VAL(TOTAL_TIME,FILE_NAME,'TOTAL_TIME',line)
     CALL GET_Float_VAL(Plot_Start,FILE_NAME,'PLOT_START',line)
     CALL GET_Float_VAL(Plot_Intv,FILE_NAME,'PLOT_INTV',line)
     CALL GET_Float_VAL(Screen_Intv,FILE_NAME,'SCREEN_INTV',line)
	 CALL GET_LOGICAL_VAL(HOTSTART,FILE_NAME,'HOTSTART',line) ! added by Cheng for hot start
     if(myid.eq.0) WRITE(3,'(A11,I12)')'SIM_STEPS= ', SIM_STEPS
     if(myid.eq.0) WRITE(3,'(A12,F8.2)')'TOTAL_TIME= ', TOTAL_TIME
     if(myid.eq.0) WRITE(3,'(A11,F8.2)')'PLOT_START= ', Plot_Start
     if(myid.eq.0) WRITE(3,'(A11,F8.2)')'PLOT_INTV= ', Plot_Intv
     if(myid.eq.0) WRITE(3,'(A13,F8.2)')'SCREEN_INTV= ', Screen_Intv
	 if(myid.eq.0) WRITE(3,'(A11,L4)')'HOTSTART= ',HOTSTART ! added by Cheng for hot start

! added by Cheng for fluid slide

! added by Cheng for deformable slide
     CALL GET_STRING_VAL(Slide_File,FILE_NAME,'SLIDE_FILE',line,ierr)
	 CALL GET_STRING_VAL(RHEO_OPT,FILE_NAME,'RHEO_OPT',line,ierr)
	 CALL GET_LOGICAL_VAL(NON_HYDRO_SLD,FILE_NAME,'NON_HYDRO_SLD',line)
	 CALL GET_LOGICAL_VAL(DISP_CORR_SLD,FILE_NAME,'DISP_CORR_SLD',line)
	 CALL GET_LOGICAL_VAL(REDU_GRAV_SLD,FILE_NAME,'REDU_GRAV_SLD',line)
	 CALL GET_LOGICAL_VAL(NON_HYDRO_UP,FILE_NAME,'NON_HYDRO_UP',line)
     CALL GET_Float_VAL(SLIDE_DENSITY ,FILE_NAME,'SLIDE_DENSITY',line)
     CALL GET_Float_VAL(SLIDE_VISCOSITY ,FILE_NAME,'SLIDE_VISCOSITY',line)
     CALL GET_Float_VAL(SLIDE_MINTHICK ,FILE_NAME,'SLIDE_MINTHICK',line)
     CALL GET_Float_VAL(SLIDE_GAMMA ,FILE_NAME,'SLIDE_GAMMA',line)
     CALL GET_Float_VAL(SLIDE_CONC ,FILE_NAME,'SLIDE_CONC',line)
     CALL GET_Float_VAL(GRAIN_DENSITY ,FILE_NAME,'GRAIN_DENSITY',line)
     CALL GET_Float_VAL(PhiInt_A ,FILE_NAME,'PhiInt_A',line)
     CALL GET_Float_VAL(PhiBed_A ,FILE_NAME,'PhiBed_A',line)
     CALL GET_Float_VAL(PhiInt_F ,FILE_NAME,'PhiInt_F',line)
     CALL GET_Float_VAL(PhiBed_F ,FILE_NAME,'PhiBed_F',line)
	 CALL GET_Float_VAL(SLIDE_LAMBDA ,FILE_NAME,'SLIDE_LAMBDA',line)
     CALL GET_Float_VAL(SLIDE_INIU ,FILE_NAME,'SLIDE_INIU',line)
	 CALL GET_Float_VAL(SLIDE_INIV ,FILE_NAME,'SLIDE_INIV',line)
	 CALL GET_Float_VAL(SLIDE_INIW ,FILE_NAME,'SLIDE_INIW',line)
       if(myid.eq.0) WRITE(3,'(A12,A15)')'Slide File= ', TRIM(Slide_File)
	   if(myid.eq.0) WRITE(3,'(A10,A10)')'RHEO_OPT= ',TRIM(RHEO_OPT)
	   if(myid.eq.0) WRITE(3,'(A15,L4)')'NON_HYDRO_SLD= ',NON_HYDRO_SLD
	   if(myid.eq.0) WRITE(3,'(A15,L4)')'DISP_CORR_SLD= ',DISP_CORR_SLD
	   if(myid.eq.0) WRITE(3,'(A15,L4)')'REDU_GRAV_SLD= ',REDU_GRAV_SLD
	   if(myid.eq.0) WRITE(3,'(A14,L4)')'NON_HYDRO_UP= ',NON_HYDRO_UP
       if(myid.eq.0) WRITE(3,'(A15,F8.2)')'SLIDE_DENSITY= ', SLIDE_DENSITY
       if(myid.eq.0) WRITE(3,'(A17,F8.2)')'SLIDE_VISCOSITY= ', SLIDE_VISCOSITY
       if(myid.eq.0) WRITE(3,'(A16,F8.2)')'SLIDE_MINTHICK= ', SLIDE_MINTHICK
       if(myid.eq.0) WRITE(3,'(A13,F8.2)')'SLIDE_GAMMA= ', SLIDE_GAMMA
       if(myid.eq.0) WRITE(3,'(A12,F8.2)')'SLIDE_CONC= ', SLIDE_CONC
       if(myid.eq.0) WRITE(3,'(A15,F8.2)')'GRAIN_DENSITY= ', GRAIN_DENSITY
       if(myid.eq.0) WRITE(3,'(A8,F10.2)')'PhiInt_A= ', PhiInt_A
       if(myid.eq.0) WRITE(3,'(A8,F10.2)')'PhiBed_A= ', PhiBed_A
       if(myid.eq.0) WRITE(3,'(A8,F10.2)')'PhiInt_F= ', PhiInt_F
       if(myid.eq.0) WRITE(3,'(A8,F10.2)')'PhiBed_F= ', PhiBed_F
	   if(myid.eq.0) WRITE(3,'(A14,F8.2)')'SLIDE_LAMBDA= ', SLIDE_LAMBDA
       if(myid.eq.0) WRITE(3,'(A12,F8.2)')'SLIDE_INIU= ', SLIDE_INIU
       if(myid.eq.0) WRITE(3,'(A12,F8.2)')'SLIDE_INIV= ', SLIDE_INIV
       if(myid.eq.0) WRITE(3,'(A12,F8.2)')'SLIDE_INIW= ', SLIDE_INIW
     PhiInt_A = PhiInt_A*pi/180.
     PhiBed_A = PhiBed_A*pi/180.
     PhiInt_F = PhiInt_F*pi/180.
     PhiBed_F = PhiBed_F*pi/180.

     ! courant number
     CALL GET_Float_VAL(CFL,FILE_NAME,'CFL',line)
     if(myid.eq.0) WRITE(3,'(A5,F8.3)')'CFL= ',CFL

     ! viscous number
     CALL GET_Float_VAL(VISCOUS_NUMBER,FILE_NAME,'VISCOUS_NUMBER',line)
     if(myid.eq.0) WRITE(3,'(A16,F8.3)')'VISCOUS_NUMBER= ',VISCOUS_NUMBER

! added by Cheng for limiting the maximum Froude number
     CALL GET_Float_VAL(FROUDECAP,FILE_NAME,'FROUDE_CAP',line)
     if(myid.eq.0) WRITE(3,'(A13,F8.3)')'FROUDE_CAP= ',FROUDECAP

! end froude cap

     ! minimum depth
     CALL GET_Float_VAL(MinDep,FILE_NAME,'MinDep',line)
     if(myid.eq.0) WRITE(3,'(A8,F8.3)')'MinDep= ',MinDep

     ! laminar viscosity
     CALL GET_LOGICAL_VAL(VISCOUS_FLOW,FILE_NAME,'VISCOUS_FLOW',line)
     CALL GET_INTEGER_VAL(IVturb,FILE_NAME,'IVTURB',line)
     CALL GET_INTEGER_VAL(IHturb,FILE_NAME,'IHTURB',line)
     CALL GET_INTEGER_VAL(ProdType,FILE_NAME,'PRODTYPE',line)
     CALL GET_Float_VAL(Visc,FILE_NAME,'VISCOSITY',line)
     CALL GET_Float_VAL(Schmidt,FILE_NAME,'Schmidt',line)
     CALL GET_Float_VAL(Cvs,FILE_NAME,'Cvs',line)
     CALL GET_Float_VAL(Chs,FILE_NAME,'Chs',line)
     CALL GET_LOGICAL_VAL(RNG,FILE_NAME,'RNG',line)
     if(myid.eq.0) WRITE(3,'(A14,L4)')'VISCOUS_FLOW= ',VISCOUS_FLOW
     if(myid.eq.0) WRITE(3,'(A8,I2)')'IVTURB= ',IVturb
     if(myid.eq.0) WRITE(3,'(A8,I2)')'IHTURB= ',IHturb
     if(myid.eq.0) WRITE(3,'(A11,F8.3)')'VISCOSITY= ',Visc
     if(myid.eq.0) WRITE(3,'(A9,F8.3)')'Schmidt= ',Schmidt
     if(myid.eq.0) WRITE(3,'(A6,F8.3)')'Cvs= ',Cvs
     if(myid.eq.0) WRITE(3,'(A6,F8.3)')'Chs= ',Chs
     if(myid.eq.0) WRITE(3,'(A5,L4)')'RNG= ',RNG

! added by Cheng for hot start
     IF(HOTSTART)THEN
       CALL GET_STRING_VAL(Eta_HotStart_File,FILE_NAME,'Eta_HotStart_File',line,ierr)
       CALL GET_STRING_VAL(U_HotStart_File,FILE_NAME,'U_HotStart_File',line,ierr)
       CALL GET_STRING_VAL(V_HotStart_File,FILE_NAME,'V_HotStart_File',line,ierr)
       CALL GET_STRING_VAL(W_HotStart_File,FILE_NAME,'W_HotStart_File',line,ierr)
       CALL GET_STRING_VAL(P_HotStart_File,FILE_NAME,'P_HotStart_File',line,ierr)
       CALL GET_STRING_VAL(Depth_HotStart_File,FILE_NAME,'Depth_HotStart_File',line,ierr)
       CALL GET_STRING_VAL(Us_HotStart_File,FILE_NAME,'Us_HotStart_File',line,ierr)
	   CALL GET_STRING_VAL(Vs_HotStart_File,FILE_NAME,'Vs_HotStart_File',line,ierr)
       CALL GET_STRING_VAL(Ws_HotStart_File,FILE_NAME,'Ws_HotStart_File',line,ierr)
      IF(VISCOUS_FLOW)THEN
       ! CALL GET_STRING_VAL(Rho_HotStart_File,FILE_NAME,'Rho_HotStart_File',line,ierr)
       CALL GET_STRING_VAL(TKE_HotStart_File,FILE_NAME,'TKE_HotStart_File',line,ierr)
       CALL GET_STRING_VAL(EPS_HotStart_File,FILE_NAME,'EPS_HotStart_File',line,ierr)
      ENDIF
     ENDIF 
! end hotstart

     ! bathymetry     
     CALL GET_STRING_VAL(DEPTH_TYPE,FILE_NAME,'DEPTH_TYPE',line,ierr)
     CALL GET_LOGICAL_VAL(ANA_BATHY,FILE_NAME,'ANA_BATHY',line)
     CALL GET_Float_VAL(DepConst,FILE_NAME,'DepConst',line)
     if(myid.eq.0) WRITE(3,'(A12,A50)')'DEPTH_TYPE= ',DEPTH_TYPE
     if(myid.eq.0) WRITE(3,'(A11,L4)')'ANA_BATHY= ',ANA_BATHY
     if(myid.eq.0) WRITE(3,'(A10,F8.5)')'DepConst= ',DepConst

     ! initial conditions
     CALL GET_LOGICAL_VAL(INITIAL_EUVW,FILE_NAME,'INITIAL_EUVW',line)
     if(myid.eq.0) WRITE(3,'(A14,L4)')'INITIAL_EUVW= ',INITIAL_EUVW

     ! bottom roughness
     CALL GET_INTEGER_VAL(Ibot,FILE_NAME,'Ibot',line)
     CALL GET_Float_VAL(Cd0,FILE_NAME,'Cd0',line)
     CALL GET_Float_VAL(Zob,FILE_NAME,'Zob',line)
     if(myid.eq.0) WRITE(3,'(A6,I2)')'Ibot= ',Ibot
     if(myid.eq.0) WRITE(3,'(A5,F8.5)')'Cd0= ',Cd0
     if(myid.eq.0) WRITE(3,'(A5,F8.5)')'Zob= ',Zob


     ! wind speed/stress
     CALL GET_INTEGER_VAL(Iws,FILE_NAME,'Iws',line)
     if(Iws==1) then
       CALL GET_Float_VAL(WindU,FILE_NAME,'WindU',line)
       CALL GET_Float_VAL(WindV,FILE_NAME,'WindV',line)
     endif

     ! Coriolis
     CALL GET_Float_VAL(slat,FILE_NAME,'slat',line)

     slat = slat*pi/180.0
     fcor = 2.0*7.29e-5*sin(slat)

     ! barotropic or baroclinic
     CALL GET_LOGICAL_VAL(BAROTROPIC,FILE_NAME,'BAROTROPIC',line)
     if(myid.eq.0) WRITE(3,'(A12,L4)')'BAROTROPIC= ',BAROTROPIC

     ! numerical scheme
     CALL GET_STRING_VAL(HIGH_ORDER,FILE_NAME,'HIGH_ORDER',line,ierr)
     CALL GET_STRING_VAL(TIME_ORDER,FILE_NAME,'TIME_ORDER',line,ierr)
     CALL GET_STRING_VAL(CONVECTION,FILE_NAME,'CONVECTION',line,ierr)
     CALL GET_LOGICAL_VAL(ADV_HLLC,FILE_NAME,'HLLC',line)
     IF(ierr==1)THEN
       if(myid.eq.0) WRITE(3,'(A12,A50)')'HIGH_ORDER', 'NOT DEFINED, USE DEFAULT'
       HIGH_ORDER='SECOND'
     ENDIF
     if(myid.eq.0) WRITE(3,'(A12,A50)')'HIGH_ORDER= ', HIGH_ORDER
     IF(ierr==1)THEN
       if(myid.eq.0) WRITE(3,'(A12,A50)')'TIME_ORDER', 'NOT DEFINED, USE DEFAULT'
       TIME_ORDER='THIRD'
     ENDIF
     if(myid.eq.0) WRITE(3,'(A12,A50)')'TIME_ORDER= ', TIME_ORDER

     ! ramp up the simulation
     CALL GET_Float_VAL(TRamp,FILE_NAME,'TRAMP',line)
     if(myid.eq.0) WRITE(3,'(A7,E12.3)')'TRAMP= ',TRamp

     ! if non-hydrostatic simulation
     CALL GET_LOGICAL_VAL(NON_HYDRO,FILE_NAME,'NON_HYDRO',line)
     if(myid.eq.0) WRITE(3,'(A11,L4)')'NON_HYDRO= ',NON_HYDRO

     ! poisson solver
     CALL GET_INTEGER_VAL(isolver,FILE_NAME,'ISOLVER',line)
     CALL GET_INTEGER_VAL(itmax,FILE_NAME,'ITMAX',line)
     CALL GET_Float_VAL(tol,FILE_NAME,'TOL',line)
     if(myid.eq.0) WRITE(3,'(A9,I2)')'ISOLVER= ',isolver
     if(myid.eq.0) WRITE(3,'(A7,I5)')'ITMAX= ',itmax
     if(myid.eq.0) WRITE(3,'(A5,E12.3)')'TOL= ',tol

     ! periodic bc
     CALL GET_LOGICAL_VAL(PERIODIC_X,FILE_NAME,'PERIODIC_X',line)
     CALL GET_LOGICAL_VAL(PERIODIC_Y,FILE_NAME,'PERIODIC_Y',line)
     if(myid.eq.0) WRITE(3,'(A12,L4)')'PERIODIC_X= ',PERIODIC_X
     if(myid.eq.0) WRITE(3,'(A12,L4)')'PERIODIC_Y= ',PERIODIC_Y

     ! boundary type
     CALL GET_INTEGER_VAL(Bc_X0,FILE_NAME,'BC_X0',line)
     CALL GET_INTEGER_VAL(Bc_Xn,FILE_NAME,'BC_Xn',line)
     CALL GET_INTEGER_VAL(Bc_Y0,FILE_NAME,'BC_Y0',line)
     CALL GET_INTEGER_VAL(Bc_Yn,FILE_NAME,'BC_Yn',line)
     CALL GET_INTEGER_VAL(Bc_Z0,FILE_NAME,'BC_Z0',line)
     CALL GET_INTEGER_VAL(Bc_Zn,FILE_NAME,'BC_Zn',line)
     if(myid.eq.0) WRITE(3,'(A7,I2)')'BC_X0= ',Bc_X0
     if(myid.eq.0) WRITE(3,'(A7,I2)')'BC_Xn= ',Bc_Xn
     if(myid.eq.0) WRITE(3,'(A7,I2)')'BC_Y0= ',Bc_Y0
     if(myid.eq.0) WRITE(3,'(A7,I2)')'BC_Yn= ',Bc_Yn
     if(myid.eq.0) WRITE(3,'(A7,I2)')'BC_Z0= ',Bc_Z0
     if(myid.eq.0) WRITE(3,'(A7,I2)')'BC_Zn= ',Bc_Zn

! wavemaker 
     CALL GET_STRING_VAL(WaveMaker,FILE_NAME,'WAVEMAKER',line,ierr)
     if(myid.eq.0) WRITE(3,'(A11,A50)')'WAVEMAKER= ', WAVEMAKER
     IF(WaveMaker(1:3)=='LEF'.or.WaveMaker(1:3)=='RIG'  &
         .or.WaveMaker(1:3)=='INT'.or.WaveMaker(1:3)=='FLU'.or.WaveMaker(1:3)=='WAV')THEN
       CALL GET_Float_VAL(Amp_Wave,FILE_NAME,'AMP',line)
       CALL GET_Float_VAL(Per_Wave,FILE_NAME,'PER',line)
       CALL GET_Float_VAL(Dep_Wave,FILE_NAME,'DEP',line)
       CALL GET_Float_VAL(Theta_Wave,FILE_NAME,'THETA',line)
       CALL GET_Float_VAL(Cur_Wave,FILE_NAME,'CUR',line)
       CALL GET_Float_VAL(sd_return,FILE_NAME,'sd_return',line)
       if(myid.eq.0) WRITE(3,'(A9,F6.3)')'AMP_WAVE= ', Amp_Wave
       if(myid.eq.0) WRITE(3,'(A9,F6.3)')'PER_WAVE= ', Per_Wave
       if(myid.eq.0) WRITE(3,'(A9,F6.3)')'DEP_WAVE= ', Dep_Wave
       if(myid.eq.0) WRITE(3,'(A12,F6.3)')'THETA_WAVE= ', Theta_Wave
       if(myid.eq.0) WRITE(3,'(A11,F6.3)')'sd_return= ',sd_return

       IF(WaveMaker(1:3)=='INT') then
         CALL GET_Float_VAL(Xsource_West,FILE_NAME,'Xsource_West',line)
         CALL GET_Float_VAL(Xsource_East,FILE_NAME,'Xsource_East',line)
         CALL GET_Float_VAL(Ysource_Suth,FILE_NAME,'Ysource_Suth',line)
         CALL GET_Float_VAL(Ysource_Nrth,FILE_NAME,'Ysource_Nrth',line)
         if(myid.eq.0) WRITE(3,'(A14,F6.3)')'Xsource_West= ',Xsource_West
         if(myid.eq.0) WRITE(3,'(A14,F6.3)')'Xsource_East= ',Xsource_East
         if(myid.eq.0) WRITE(3,'(A14,F6.3)')'Ysource_Suth= ',Ysource_Suth
         if(myid.eq.0) WRITE(3,'(A14,F6.3)')'Ysource_Nrth= ',Ysource_Nrth
       ENDIF

! test periodicity
       IF(PERIODIC_Y.and.Theta_Wave.ne.Zero) then
         ! find wave number
         Segma = 2.0*pi/Per_Wave
         Celerity = sqrt(Grav*Dep_Wave)
         Wave_Length = Celerity*Per_Wave
         Wave_Number = 2.0*pi/Wave_Length

         Iter = 0
 75      Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
         if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 85
         Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  &                
            Grav*tanh(Wave_Number*Dep_Wave)
         Wave_Number = Wave_Number-Fk/Fkdif
         Iter = Iter+1
         goto 75
 85      continue

         if(Theta_Wave>Zero) then
           ! find right angle for periodic bc
           tmp = Large       
           do Iter = 1,10000
             Wnumy = Iter*2.0*pi/(Nglob*dy)
             if(WnumY<Wave_Number) then
               ! theta based on Ky = K*sin(theta)
               tmp1 = asin(Wnumy/Wave_Number)*180./pi
               if(abs(tmp1-Theta_Wave)<tmp) then
                 tmp = abs(tmp1-Theta_Wave)
                 Theta_Calc = tmp1
               endif
             endif
           enddo
         elseif(Theta_Wave<Zero) then
           ! find right angle for periodic bc 
           tmp = Large
           do Iter = 1,10000
             Wnumy = Iter*2.0*pi/(Nglob*dy)
             if(WnumY<Wave_Number) then
               ! theta based on Ky = K*sin(theta)
               tmp1 = -asin(Wnumy/Wave_Number)*180./pi
               if(abs(tmp1-Theta_Wave)<tmp) then
                 tmp = abs(tmp1-Theta_Wave)
                 Theta_Calc = tmp1
               endif
             endif
           enddo
         endif

         if(myid.eq.0) then
           write(3,'(A20,F6.3)') 'Wave angle you set= ',Theta_Wave         
           write(3,'(A28,F6.3)') 'Wave angle for periodic bc= ',Theta_Calc
         endif
         Theta_Wave = Theta_Calc
       ENDIF
     ENDIF

 
!added by M.Derakhti for FOCUSED wavemaker
     IF(WaveMaker(1:3)=='FOC') then
         CALL GET_INTEGER_VAL(nwave,FILE_NAME,'nwave',line)
         CALL GET_INTEGER_VAL(Component_Amp_Type,FILE_NAME,'Component_Amp_Type',line)
         CALL GET_Float_VAL(k_center,FILE_NAME,'k_center',line)
         CALL GET_Float_VAL(f_center,FILE_NAME,'f_center',line)
         CALL GET_Float_VAL(x_breaking,FILE_NAME,'x_breaking',line)
         CALL GET_Float_VAL(t_breaking,FILE_NAME,'t_breaking',line)
         CALL GET_Float_VAL(Slope_group,FILE_NAME,'Slope_group',line)
         CALL GET_Float_VAL(normalized_delta_f,FILE_NAME,'normalized_delta_f',line)
         CALL GET_Float_VAL(depth_comp,FILE_NAME,'depth_comp',line)

         if(myid.eq.0) WRITE(3,'(A7,I2)')'nwave= ',nwave
         if(myid.eq.0) WRITE(3,'(A20,I2)')'Component_Amp_Type= ',Component_Amp_Type
         if(myid.eq.0) WRITE(3,'(A10,F6.3)')'k_center= ',k_center
         if(myid.eq.0) WRITE(3,'(A10,F6.3)')'f_center= ',f_center
         if(myid.eq.0) WRITE(3,'(A12,F6.3)')'x_breaking= ',x_breaking
         if(myid.eq.0) WRITE(3,'(A12,F6.3)')'t_breaking= ',t_breaking
         if(myid.eq.0) WRITE(3,'(A13,F6.3)')'Slope_group= ',Slope_group
         if(myid.eq.0) WRITE(3,'(A20,F6.3)')'normalized_delta_f= ',normalized_delta_f
         if(myid.eq.0) WRITE(3,'(A12,F6.3)')'depth_comp= ',depth_comp

         !calculate f_i based on packet center frequency and band width
         period_i = 0.0
         do nw = 1 , nwave
            f_i (nw) = f_center - 0.5 * normalized_delta_f * f_center + &
              (nw - 1.0) * normalized_delta_f * f_center / (nwave - 1.0)
            period_i (nw)  = 1.0_SP / f_i (nw)
         end do
!calculate the wave number based on dispersion relation using
!Newton Ralphson method,xxk means k
         k_i = 0.0
         do nw = 1 , nwave
            segmasqr = ( 2.0d0 * pi / period_i (nw) ) ** 2
            ctmp = dsqrt ( Grav * depth_comp)
            xltmp = ctmp * period_i ( nw )
            xxk = 2.0d0 * pi / xltmp
            n = 0
 355        Fk = Grav * xxk * tanh ( xxk * depth_comp ) - segmasqr
            if ( abs ( Fk ) .le. 1.0e-6 .or. n .gt. 1000 ) goto 365
            Fkdif = Grav * xxk * depth_comp * ( 1.0d0 - &
            tanh ( xxk * depth_comp) ** 2 ) &
            + Grav * tanh( xxk * depth_comp )
            xxk = xxk - Fk / Fkdif
            n = n+1
            goto 355
 365        continue
            k_i (nw) = xxk
         end do
         !calculate componets amplitude
         Amp_i = 0.0
         If (Component_Amp_Type == 1) then
            do nw = 1 , nwave
               Amp_i (nw) = Slope_group / nwave / k_center
            end do
         else
            do nw = 1 , nwave
               Amp_i (nw) = Slope_group / nwave / k_i(nw)
            end do
         endif
     ENDIF

     ! random wave, read in 2d spectrum
     IF(WaveMaker(5:7)=='SPC') then
       open(14,file='spc2d.txt')
       read(14,*) NumFreq,NumDir
       if(NumFreq>MaxNumFreq) then
         if(myid.eq.0) then
           write(3,'(A)') 'Please set a larger MaxNumFreq in mod_glob.F'
           stop
         endif
       endif
       if(NumDir>MaxNumDir) then
         if(myid.eq.0) then
           write(3,'(A)') 'Please set a larger MaxNumDir in mod_glob.F'
           stop
         endif
       endif
       do i = 1,NumFreq
         read(14,*) Freq(i)
       enddo
       do i = 1,NumDir
         read(14,*) Dire(i)
       enddo
       do j = 1,NumFreq
       do i = 1,NumDir
         read(14,*) Wave_Spc2d(i,j)
       enddo
       enddo
       close(14)

       ! random phase for each component
       do j = 1,NumFreq
       do i = 1,NumDir
         Random_Phs(i,j) = rand(0)*2.0*pi
       enddo
       enddo
     ENDIF
	 
     ! Added by Cheng Zhang,irregular wave internal wavemaker
     IF(WaveMaker(5:7)=='IRR') then
       open(14,file='irr2d.txt')
       read(14,*) NumFreq,NumDir
       if(NumFreq>MaxNumFreq) then
         if(myid.eq.0) then
           write(3,'(A)') 'Please set a larger MaxNumFreq in mod_glob.F'
           stop
         endif
       endif
       if(NumDir>MaxNumDir) then
         if(myid.eq.0) then
           write(3,'(A)') 'Please set a larger MaxNumDir in mod_glob.F'
           stop
         endif
       endif
       do i = 1,NumFreq
         read(14,*) Freq(i)
       enddo
       do i = 1,NumDir
         read(14,*) Dire(i)
       enddo
       do j = 1,NumFreq
       do i = 1,NumDir
         read(14,*) Wave_Spc2d(i,j)
       enddo
       enddo
	   do j = 1,NumFreq
       do i = 1,NumDir
         read(14,*) Random_Phs(i,j)
       enddo
       enddo
       close(14)
     ENDIF

! JONSWAP spectrum
     if((WaveMaker(5:7)=='JON').or.(WaveMaker(5:7)=='TMA')) then
       CALL GET_Float_VAL(Hm0,FILE_NAME,'Hm0',line)
       CALL GET_Float_VAL(Tp,FILE_NAME,'Tp',line)
       CALL GET_Float_VAL(Freq_Min,FILE_NAME,'Freq_Min',line)
       CALL GET_Float_VAL(Freq_Max,FILE_NAME,'Freq_Max',line)
       CALL GET_INTEGER_VAL(NumFreq,FILE_NAME,'NumFreq',line) 
       if(myid.eq.0) WRITE(3,'(A5,f6.2)')'Hm0= ', Hm0
       if(myid.eq.0) WRITE(3,'(A4,f6.2)')'Tp= ', Tp
       if(myid.eq.0) WRITE(3,'(A10,f6.2)')'Freq_Min= ', Freq_Min
       if(myid.eq.0) WRITE(3,'(A10,f6.2)')'Freq_Max= ', Freq_Max
       if(myid.eq.0) WRITE(3,'(A9,I5)')'NumFreq= ', NumFreq

! jonswap spectrum
       gam = 3.3; sa = 0.07; sb = 0.09
       Freq_Peak = 1.0/Tp

       DFreq = (Freq_Max-Freq_Min)/NumFreq
       do i = 1,NumFreq
         Freq(i) = Freq_Min+0.5*DFreq+(i-1)*DFreq

         Per_Wave = 1.0/Freq(i)
         Segma = 2.0*pi/Per_Wave
         Celerity = sqrt(Grav*Dep_Wave)
         Wave_Length = Celerity*Per_Wave
         Wave_Number = 2.0*pi/Wave_Length
       
         Iter = 0
 76      Fk = Grav*Wave_Number*tanh(Wave_Number*Dep_Wave)-Segma**2
         if(abs(Fk)<=1.0e-8.or.Iter>1000) goto 86
         Fkdif = Grav*Wave_Number*Dep_Wave*(1.0-tanh(Wave_Number*Dep_Wave)**2)+  & 
                 Grav*tanh(Wave_Number*Dep_Wave)
         Wave_Number = Wave_Number-Fk/Fkdif
         Iter = Iter+1
         goto 76
 86      continue

         if(Freq(i)<Freq_Peak) then
           Jon_Spc(i) = Grav**2/Freq(i)**5*exp(-1.25*(Freq_Peak/Freq(i))**4)*  &
               gam**exp(-0.5*(Freq(i)/Freq_Peak-1.0)**2/sa**2)
         else
           Jon_Spc(i) = Grav**2/Freq(i)**5*exp(-1.25*(Freq_Peak/Freq(i))**4)*  &
               gam**exp(-0.5*(Freq(i)/Freq_Peak-1.0)**2/sb**2)
         endif

         if(WaveMaker(5:7)=='TMA') then
            Jon_Spc(i) = Jon_Spc(i)*tanh(Wave_Number*Dep_Wave)**2/  &
                 (1.0+2*Wave_Number*Dep_Wave/sinh(2.*Wave_Number*Dep_Wave))
         endif
       enddo
         
       ! make sure m0=Hm0**2/16=int S(f)df
       SumInt = Zero
       do i = 1,NumFreq
         SumInt = SumInt+Jon_Spc(i)*DFreq
       enddo
       A_Jon = Hm0**2/16.0/SumInt

       do i = 1,NumFreq
         Jon_Spc(i) = Jon_Spc(i)*A_Jon
         RanPhs(i) = rand()*2.0*pi
       enddo
     endif

! sponge layer
     CALL GET_LOGICAL_VAL(SPONGE_ON,FILE_NAME,'SPONGE_ON',line)
     if(myid.eq.0) WRITE(3,'(A11,L4)')'SPONGE_ON= ', SPONGE_ON
     IF(SPONGE_ON)THEN
       CALL GET_Float_VAL(Sponge_West_Width,FILE_NAME,'Sponge_West_Width',line)
       CALL GET_Float_VAL(Sponge_East_Width,FILE_NAME,'Sponge_East_Width',line)
       CALL GET_Float_VAL(Sponge_South_Width,FILE_NAME,'Sponge_South_Width',line)
       CALL GET_Float_VAL(Sponge_North_Width,FILE_NAME,'Sponge_North_Width',line)
       if(myid.eq.0) WRITE(3,'(A19,F6.3)')'Sponge_West_Width= ', Sponge_West_Width
       if(myid.eq.0) WRITE(3,'(A19,F6.3)')'Sponge_East_Width= ', Sponge_East_Width
       if(myid.eq.0) WRITE(3,'(A20,F6.3)')'Sponge_South_Width= ', Sponge_South_Width
       if(myid.eq.0) WRITE(3,'(A20,F6.3)')'Sponge_North_Width= ', Sponge_North_Width
     ENDIF

     ! wave average control
     CALL GET_LOGICAL_VAL(WAVE_AVERAGE_ON,FILE_NAME,'WAVE_AVERAGE_ON',line)
     CALL GET_Float_VAL(Wave_Ave_Start,FILE_NAME,'WAVE_AVERAGE_START',line)
     CALL GET_Float_VAL(Wave_Ave_End,FILE_NAME,'WAVE_AVERAGE_END',line)
     CALL GET_INTEGER_VAL(WaveheightID,FILE_NAME,'WaveheightID',line)






! end landslide comprehensive

     
     ! whether to consider rheology
     CALL GET_LOGICAL_VAL(RHEOLOGY_ON,FILE_NAME,'RHEOLOGY_ON',line)
     CALL GET_Float_VAL(Yield_Stress,FILE_NAME,'Yield_Stress',line)
     CALL GET_Float_VAL(Plastic_Visc,FILE_NAME,'Plastic_Visc',line)


     ! if there is external forcing
     CALL GET_LOGICAL_VAL(EXTERNAL_FORCING,FILE_NAME,'EXTERNAL_FORCING',line)
     if(EXTERNAL_FORCING) then
       CALL GET_Float_VAL(Pgrad0,FILE_NAME,'Pgrad0',line)
     endif

     ! probe output
     CALL GET_INTEGER_VAL(NSTAT,FILE_NAME,'NSTAT',line)
     CALL GET_Float_VAL(Plot_Intv_Stat,FILE_NAME,'PLOT_INTV_STAT',line)
       if(myid.eq.0) WRITE(3,'(A7,I3)')'NSTAT= ', NSTAT
       if(myid.eq.0) WRITE(3,'(A19,F6.3)')'Plot_Intv_Stat= ', Plot_Intv_Stat

     if(NSTAT>0) then
       open(15,file='stat.txt',status='old')
       do i = 1,NSTAT
         read(15,*) xstat(i),ystat(i),zstat(i)
       enddo
       close(15)
     endif

     ! output parameters  
     CALL GET_LOGICAL_VAL(OUT_H,FILE_NAME,'OUT_H',line)
     CALL GET_LOGICAL_VAL(OUT_E,FILE_NAME,'OUT_E',line)
     CALL GET_LOGICAL_VAL(OUT_U,FILE_NAME,'OUT_U',line)
     CALL GET_LOGICAL_VAL(OUT_V,FILE_NAME,'OUT_V',line)
     CALL GET_LOGICAL_VAL(OUT_W,FILE_NAME,'OUT_W',line)
     CALL GET_LOGICAL_VAL(OUT_P,FILE_NAME,'OUT_P',line)
     CALL GET_LOGICAL_VAL(OUT_K,FILE_NAME,'OUT_K',line)
     CALL GET_LOGICAL_VAL(OUT_D,FILE_NAME,'OUT_D',line)
     CALL GET_LOGICAL_VAL(OUT_S,FILE_NAME,'OUT_S',line)
     CALL GET_LOGICAL_VAL(OUT_C,FILE_NAME,'OUT_C',line)
     CALL GET_LOGICAL_VAL(OUT_B,FILE_NAME,'OUT_B',line)
     CALL GET_LOGICAL_VAL(OUT_A,FILE_NAME,'OUT_A',line)
     CALL GET_LOGICAL_VAL(OUT_F,FILE_NAME,'OUT_F',line)
     CALL GET_LOGICAL_VAL(OUT_T,FILE_NAME,'OUT_T',line)
     CALL GET_LOGICAL_VAL(OUT_G,FILE_NAME,'OUT_G',line)
     CALL GET_LOGICAL_VAL(OUT_I,FILE_NAME,'OUT_I',line)
	 CALL GET_LOGICAL_VAL(OUT_Z,FILE_NAME,'OUT_Z',line) !added by Cheng for varying depth
	 CALL GET_LOGICAL_VAL(OUT_M,FILE_NAME,'OUT_M',line) !added by Cheng for recording Hmax
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_H= ',OUT_H
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_E= ',OUT_E
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_U= ',OUT_U
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_V= ',OUT_V
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_W= ',OUT_W
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_P= ',OUT_P
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_K= ',OUT_K
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_D= ',OUT_D
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_S= ',OUT_S
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_C= ',OUT_C
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_B= ',OUT_B
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_A= ',OUT_A
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_F= ',OUT_F
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_T= ',OUT_T
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_G= ',OUT_G
     if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_I= ',OUT_I
	 if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_Z= ',OUT_Z !added by Cheng for varying depth
	 if(myid.eq.0) WRITE(3,'(A7,L4)') 'OUT_M= ',OUT_M !added by Cheng for recording Hmax

! added by Cheng for nesting

     if(myid.eq.0) WRITE(3,*)'--------------input end --------------'
!
     end subroutine read_input
!
!---------------------------------------------------------------------------------------
!
     subroutine index
!
!----------------------------------------------------------------------------------------
!
!    (3) index
!
!    This subroutine is used to create work index (what's that?)
!
!    Called by: main  
!                                                              
!    Last update: 20/12/2010, Gangfeng Ma       
!                              
!---------------------------------------------------------------------------------------
!
     use global
     implicit none

     dims(1)=PX
     dims(2)=PY
     periods(1)=.false.
     periods(2)=.false.
     if(PERIODIC_X) periods(1)=.true.
     if(PERIODIC_Y) periods(2)=.true.
     coords(1)=0
     coords(2)=0

     call MPI_CART_CREATE(MPI_COMM_WORLD,ndims,dims, &
         periods,reorder,comm2d,ier)
     call MPI_CART_COORDS(comm2d,myid,2,coords,ier)

     npx=coords(1)
     npy=coords(2)
 
     call MPI_CART_SHIFT(comm2d,0,1,n_west,n_east,ier)
     call MPI_CART_SHIFT(comm2d,1,1,n_suth,n_nrth,ier)

     ! local index
     Mloc = Mglob/PX+2*Nghost
     Nloc = Nglob/PY+2*Nghost
     Kloc = Kglob+2*Nghost
     Mloc1 = Mloc+1
     Nloc1 = Nloc+1
     Kloc1 = Kloc+1

     Ibeg = Nghost+1
     Iend = Mloc-Nghost
     Iend1 = Mloc1-Nghost
     Jbeg = Nghost+1
     Jend = Nloc-Nghost
     Jend1 = Nloc1-Nghost
     Kbeg = Nghost+1
     Kend = Kloc-Nghost
     Kend1 = Kloc1-Nghost

	 if (npx<(PX-1).and.npy<(PY-1)) then
	   IendC=Iend;JendC=Jend
	   MlocC=Mloc;NlocC=Nloc
	 elseif (npx==(PX-1).and.npy<(PY-1)) then
	   IendC=Iend1;JendC=Jend
       MlocC=Mloc1;NlocC=Nloc
	 elseif (npx<(PX-1).and.npy==(PY-1)) then
	   IendC=Iend;JendC=Jend1
       MlocC=Mloc;NlocC=Nloc1
	 elseif (npx==(PX-1).and.npy==(PY-1)) then
	   IendC=Iend1;JendC=Jend1
       MlocC=Mloc1;NlocC=Nloc1
	 endif
!
     end subroutine index
!
!-----------------------------------------------------------------------------------------
!
     subroutine allocate_variables
!
!----------------------------------------------------------------------------------------
!
!    (4) allocate_variables
! 
!    This subroutine is used to allocate variables
!
!    Called by: main  
!                                                                                
!    Last update: 23/12/2010, Gangfeng Ma    
!                                               
!----------------------------------------------------------------------------------------
!
     use global
     implicit none

     ! one-dimensional vars
     ALLOCATE(x(Mloc1),xc(Mloc),y(Nloc1),yc(Nloc),sig(Kloc1),dsig(Kloc),sigc(Kloc),  &
              Ein_X0(Nloc),Din_X0(Nloc),Ein_Xn(Nloc),Din_Xn(Nloc))

     ! two-dimensional vars
     ALLOCATE(HCG(Mglob,Nglob),Ho(Mloc,Nloc),H(Mloc1,Nloc1),Hc(Mloc,Nloc),Hc0(Mloc,Nloc), &
              Hfx(Mloc1,Nloc),Hfy(Mloc,Nloc1),Hfx0(Mloc1,Nloc),Hfy0(Mloc,Nloc1),D(Mloc,Nloc),Eta0(Mloc,Nloc),Eta00(Mloc,Nloc),  &
              D0(Mloc,Nloc),DeltH(Mloc,Nloc),DelxH(Mloc,Nloc),DelyH(Mloc,Nloc),Eta(Mloc,Nloc),Mask(Mloc,Nloc),  &
              Mask_Struct(Mloc,Nloc),Mask9(Mloc,Nloc),SourceC(Mloc,Nloc),SourceX(Mloc,Nloc), &
              SourceY(Mloc,Nloc),DeltHo(Mloc,Nloc),Uin_X0(Nloc,Kloc),Vin_X0(Nloc,Kloc),  &
              Win_X0(Nloc,Kloc),Uin_Xn(Nloc,Kloc),Uin_Xni(Nloc,Kloc),Uin_Xni0(Nloc,Kloc),Delt2H(Mloc,Nloc),  &
              Vin_Xn(Nloc,Kloc),Win_Xn(Nloc,Kloc),Bc_Prs(Mloc,Nloc),Brks(Mloc,Nloc))
     ALLOCATE(DxL(Mloc1,Nloc),DxR(Mloc1,Nloc),DyL(Mloc,Nloc1),DyR(Mloc,Nloc1), &
              EtaxL(Mloc1,Nloc),EtaxR(Mloc1,Nloc),EtayL(Mloc,Nloc1),EtayR(Mloc,Nloc1), &
              DelxEta(Mloc,Nloc),DelyEta(Mloc,Nloc),DelxD(Mloc,Nloc),DelyD(Mloc,Nloc),Sponge(Mloc,Nloc), &
              Setup(Mloc,Nloc),WaveHeight(Mloc,Nloc),Umean(Mloc,Nloc),Vmean(Mloc,Nloc),Num_Zero_Up(Mloc,Nloc), &
              Emax(Mloc,Nloc),Emin(Mloc,Nloc),WdU(Mloc,Nloc),WdV(Mloc,Nloc),Wsx(Mloc,Nloc),Wsy(Mloc,Nloc), &
			  HeightMax(Mloc,Nloc))

     ! three-dimensional vars


     ALLOCATE(U(Mloc,Nloc,Kloc),V(Mloc,Nloc,Kloc),W(Mloc,Nloc,Kloc),Omega(Mloc,Nloc,Kloc1), &
              P(Mloc,Nloc,Kloc1),DU(Mloc,Nloc,Kloc),DV(Mloc,Nloc,Kloc),DW(Mloc,Nloc,Kloc),  &
              U0(Mloc,Nloc,Kloc),V0(Mloc,Nloc,Kloc),W0(Mloc,Nloc,Kloc),  &
              U00(Mloc,Nloc,Kloc),V00(Mloc,Nloc,Kloc),W00(Mloc,Nloc,Kloc),  &
              DU0(Mloc,Nloc,Kloc),DV0(Mloc,Nloc,Kloc),DW0(Mloc,Nloc,Kloc),Uf(Mloc,Nloc,Kloc1), &
              Vf(Mloc,Nloc,Kloc1),Wf(Mloc,Nloc,Kloc1),Cmu(Mloc,Nloc,Kloc),CmuR(Mloc,Nloc,Kloc), &
              Diffxx(Mloc,Nloc,Kloc),Diffxy(Mloc,Nloc,Kloc),Diffxz(Mloc,Nloc,Kloc), &
              Diffyx(Mloc,Nloc,Kloc),Diffyy(Mloc,Nloc,Kloc),Diffyz(Mloc,Nloc,Kloc),Diffzx(Mloc,Nloc,Kloc), &
              Diffzy(Mloc,Nloc,Kloc),Diffzz(Mloc,Nloc,Kloc),DelxSc(Mloc,Nloc,Kloc),DelySc(Mloc,Nloc,Kloc), &
              CmuHt(Mloc,Nloc,Kloc),CmuVt(Mloc,Nloc,Kloc),Rho(Mloc,Nloc,Kloc),Rmean(Mloc,Nloc,Kloc),Tke(Mloc,Nloc,Kloc), &
              Eps(Mloc,Nloc,Kloc),Skl(Mloc,Nloc,Kloc),DTke(Mloc,Nloc,Kloc),DEps(Mloc,Nloc,Kloc),DTke0(Mloc,Nloc,Kloc), &
              DEps0(Mloc,Nloc,Kloc),Prod_s(Mloc,Nloc,Kloc),Prod_b(Mloc,Nloc,Kloc),Lag_Umean(Mloc,Nloc,Kloc), &
              Lag_Vmean(Mloc,Nloc,Kloc),Lag_Wmean(Mloc,Nloc,Kloc),Euler_Umean(Mloc,Nloc,Kloc),Euler_Vmean(Mloc,Nloc,Kloc), &
              Euler_Wmean(Mloc,Nloc,Kloc),DRhoX(Mloc,Nloc,Kloc),DRhoY(Mloc,Nloc,Kloc),ExtForceX(Mloc,Nloc,Kloc), &
              ExtForceY(Mloc,Nloc,Kloc),UpWp(Mloc,Nloc,Kloc),IsMove(Mloc,Nloc,Kloc),Richf(Mloc,Nloc,Kloc), &
              DelxSl(Mloc,Nloc,Kloc1),DelySl(Mloc,Nloc,Kloc1))

     ! fluxes for construction at cell faces    
     ALLOCATE(UxL(Mloc1,Nloc,Kloc),UxR(Mloc1,Nloc,Kloc),VxL(Mloc1,Nloc,Kloc),VxR(Mloc1,Nloc,Kloc), &
              WxL(Mloc1,Nloc,Kloc),WxR(Mloc1,Nloc,Kloc),DUxL(Mloc1,Nloc,Kloc),DUxR(Mloc1,Nloc,Kloc), &
              DVxL(Mloc1,Nloc,Kloc),DVxR(Mloc1,Nloc,Kloc),DWxL(Mloc1,Nloc,Kloc),DWxR(Mloc1,Nloc,Kloc), &
              UyL(Mloc,Nloc1,Kloc),UyR(Mloc,Nloc1,Kloc),VyL(Mloc,Nloc1,Kloc),VyR(Mloc,Nloc1,Kloc), &
              WyL(Mloc,Nloc1,Kloc),WyR(Mloc,Nloc1,Kloc),DUyL(Mloc,Nloc1,Kloc),DUyR(Mloc,Nloc1,Kloc), &
              DVyL(Mloc,Nloc1,Kloc),DVyR(Mloc,Nloc1,Kloc),DWyL(Mloc,Nloc1,Kloc),DWyR(Mloc,Nloc1,Kloc), &
              UzL(Mloc,Nloc,Kloc1),UzR(Mloc,Nloc,Kloc1),VzL(Mloc,Nloc,Kloc1),VzR(Mloc,Nloc,Kloc1), &
              WzL(Mloc,Nloc,Kloc1),WzR(Mloc,Nloc,Kloc1),OzL(Mloc,Nloc,Kloc1),OzR(Mloc,Nloc,Kloc1), &
              SxL(Mloc1,Nloc,Kloc),SxR(Mloc1,Nloc,Kloc),SxS(Mloc1,Nloc,Kloc), &
              SyL(Mloc,Nloc1,Kloc),SyR(Mloc,Nloc1,Kloc),SyS(Mloc,Nloc1,Kloc), &
              ExL(Mloc1,Nloc,Kloc),ExR(Mloc1,Nloc,Kloc),FxL(Mloc1,Nloc,Kloc),FxR(Mloc1,Nloc,Kloc), &
              GxL(Mloc1,Nloc,Kloc),GxR(Mloc1,Nloc,Kloc),HxL(Mloc1,Nloc,Kloc),HxR(Mloc1,Nloc,Kloc), &
              EyL(Mloc,Nloc1,Kloc),EyR(Mloc,Nloc1,Kloc),FyL(Mloc,Nloc1,Kloc),FyR(Mloc,Nloc1,Kloc), &
              GyL(Mloc,Nloc1,Kloc),GyR(Mloc,Nloc1,Kloc),HyL(Mloc,Nloc1,Kloc),HyR(Mloc,Nloc1,Kloc), &
              Ex(Mloc1,Nloc,Kloc),Ey(Mloc,Nloc1,Kloc),Fx(Mloc1,Nloc,Kloc),Fy(Mloc,Nloc1,Kloc), &
              Gx(Mloc1,Nloc,Kloc),Gy(Mloc,Nloc1,Kloc),Hx(Mloc1,Nloc,Kloc),Hy(Mloc,Nloc1,Kloc), &
              Fz(Mloc,Nloc,Kloc1),Gz(Mloc,Nloc,Kloc1),Hz(Mloc,Nloc,Kloc1),DelxU(Mloc,Nloc,Kloc), &
              DelyU(Mloc,Nloc,Kloc),DelzU(Mloc,Nloc,Kloc),DelxV(Mloc,Nloc,Kloc),DelyV(Mloc,Nloc,Kloc), &
              DelzV(Mloc,Nloc,Kloc),DelxW(Mloc,Nloc,Kloc),DelyW(Mloc,Nloc,Kloc),DelzW(Mloc,Nloc,Kloc), &
              DelxDU(Mloc,Nloc,Kloc),DelyDU(Mloc,Nloc,Kloc),DelxDV(Mloc,Nloc,Kloc),DelyDV(Mloc,Nloc,Kloc), &
              DelxDW(Mloc,Nloc,Kloc),DelyDW(Mloc,Nloc,Kloc),DelzO(Mloc,Nloc,Kloc)) 










! added by Cheng for fluid slide

! added by Cheng for deformable slide
     ALLOCATE(Ugs(Mloc,Nloc),Vgs(Mloc,Nloc),Wgs(Mloc,Nloc),DUgs(Mloc,Nloc),DVgs(Mloc,Nloc),DWgs(Mloc,Nloc), &
	          Dgs(Mloc,Nloc),Rhogs(Mloc,Nloc),Hgs(Mloc,Nloc),Wtgs(Mloc,Nloc),Qtgs(Mloc,Nloc), &
	          Ugs0(Mloc,Nloc),Vgs0(Mloc,Nloc),Wgs0(Mloc,Nloc),DUgs0(Mloc,Nloc),DVgs0(Mloc,Nloc),DWgs0(Mloc,Nloc), &
			  Dgs0(Mloc,Nloc),Rhogs0(Mloc,Nloc),Wgs00(Mloc,Nloc),Wbgs(Mloc,Nloc),Qbgs(Mloc,Nloc), &
              DelxUgs(Mloc,Nloc),DelxVgs(Mloc,Nloc),DelxWgs(Mloc,Nloc), &
			  DelyUgs(Mloc,Nloc),DelyVgs(Mloc,Nloc),DelyWgs(Mloc,Nloc), &
			  DelxDUgs(Mloc,Nloc),DelxDVgs(Mloc,Nloc),DelxDWgs(Mloc,Nloc),DelxDgs(Mloc,Nloc), &
			  DelyDUgs(Mloc,Nloc),DelyDVgs(Mloc,Nloc),DelyDWgs(Mloc,Nloc),DelyDgs(Mloc,Nloc), &
			  DelxHgs(Mloc,Nloc),DelxH0(Mloc,Nloc),Delx2H0(Mloc,Nloc),Kx(Mloc,Nloc),Cxgs(Mloc,Nloc), &
			  DelyHgs(Mloc,Nloc),DelyH0(Mloc,Nloc),Dely2H0(Mloc,Nloc),Ky(Mloc,Nloc),Cygs(Mloc,Nloc), &
			  SrcmgsX(Mloc,Nloc),SrcpgsX(Mloc,Nloc),SrctgsX(Mloc,Nloc), &
			  SrcmgsY(Mloc,Nloc),SrcpgsY(Mloc,Nloc),SrctgsY(Mloc,Nloc), &
			  SrctgsZ(Mloc,Nloc),UgsP(Mloc,Nloc),VgsP(Mloc,Nloc),WgsP(Mloc,Nloc),DwDt(Mloc,Nloc), &
			  Cxx(Mloc,Nloc),Cxy(Mloc,Nloc),Cxz(Mloc,Nloc),Cyx(Mloc,Nloc),Cyy(Mloc,Nloc), &
			  Cyz(Mloc,Nloc),Czx(Mloc,Nloc),Czy(Mloc,Nloc),Czz(Mloc,Nloc),Czz0(Mloc,Nloc), &
			  Taxx(Mloc,Nloc),Tayx(Mloc,Nloc),Taxy(Mloc,Nloc),Tayy(Mloc,Nloc),Taxz(Mloc,Nloc),Tayz(Mloc,Nloc), &
			  Tezz(Mloc,Nloc),Pss(Mloc,Nloc),Qsgs(Mloc,Nloc), &
			  Tbxx(Mloc,Nloc),Tbxy(Mloc,Nloc),Tbxz(Mloc,Nloc),Tbyx(Mloc,Nloc),Tbyy(Mloc,Nloc), &
			  Tbyz(Mloc,Nloc),Tbzx(Mloc,Nloc),Tbzy(Mloc,Nloc),Tbzz(Mloc,Nloc),Maskgs(Mloc,Nloc), &
			  PhiInt(Mloc,Nloc),PhiBed(Mloc,Nloc))
     ALLOCATE(DgsxL(Mloc1,Nloc),DgsxR(Mloc1,Nloc),HSgsxL(Mloc1,Nloc),HSgsxR(Mloc1,Nloc),H0fx(Mloc1,Nloc), &
	          UgsxL(Mloc1,Nloc),UgsxR(Mloc1,Nloc),VgsxL(Mloc1,Nloc),VgsxR(Mloc1,Nloc),WgsxL(Mloc1,Nloc),WgsxR(Mloc1,Nloc), &
              DUgsxL(Mloc1,Nloc),DUgsxR(Mloc1,Nloc),DVgsxL(Mloc1,Nloc),DVgsxR(Mloc1,Nloc),DWgsxL(Mloc1,Nloc),DWgsxR(Mloc1,Nloc), &
			  Cxgsx(Mloc1,Nloc),SgsxL(Mloc1,Nloc),SgsxR(Mloc1,Nloc), &
			  Egsx(Mloc1,Nloc),Fgsx(Mloc1,Nloc),Ggsx(Mloc1,Nloc),Hgsx(Mloc1,Nloc), &
			  EgsxL(Mloc1,Nloc),EgsxR(Mloc1,Nloc),FgsxL(Mloc1,Nloc),FgsxR(Mloc1,Nloc), &
			  GgsxL(Mloc1,Nloc),GgsxR(Mloc1,Nloc),HgsxL(Mloc1,Nloc),HgsxR(Mloc1,Nloc), &
			  SrcmgsxL(Mloc1,Nloc),SrcmgsxR(Mloc1,Nloc),SrcpgsxL(Mloc1,Nloc),SrcpgsxR(Mloc1,Nloc))
     ALLOCATE(DgsyL(Mloc,Nloc1),DgsyR(Mloc,Nloc1),HSgsyL(Mloc,Nloc1),HSgsyR(Mloc,Nloc1),H0fy(Mloc,Nloc1), &
	          UgsyL(Mloc,Nloc1),UgsyR(Mloc,Nloc1),VgsyL(Mloc,Nloc1),VgsyR(Mloc,Nloc1),WgsyL(Mloc,Nloc1),WgsyR(Mloc,Nloc1), &
              DUgsyL(Mloc,Nloc1),DUgsyR(Mloc,Nloc1),DVgsyL(Mloc,Nloc1),DVgsyR(Mloc,Nloc1),DWgsyL(Mloc,Nloc1),DWgsyR(Mloc,Nloc1), &
			  Cygsy(Mloc,Nloc1),SgsyL(Mloc,Nloc1),SgsyR(Mloc,Nloc1), &
			  Egsy(Mloc,Nloc1),Fgsy(Mloc,Nloc1),Ggsy(Mloc,Nloc1),Hgsy(Mloc,Nloc1), &
			  EgsyL(Mloc,Nloc1),EgsyR(Mloc,Nloc1),FgsyL(Mloc,Nloc1),FgsyR(Mloc,Nloc1), &
			  GgsyL(Mloc,Nloc1),GgsyR(Mloc,Nloc1),HgsyL(Mloc,Nloc1),HgsyR(Mloc,Nloc1), &
			  SrcmgsyL(Mloc,Nloc1),SrcmgsyR(Mloc,Nloc1),SrcpgsyL(Mloc,Nloc1),SrcpgsyR(Mloc,Nloc1))
     ALLOCATE(QbgsC(MlocC,NlocC),WgsC(MlocC,NlocC))
! poisson solver (for NSPCG use)
     neqnsgs = (IendC-Ibeg+1)*(JendC-Jbeg+1)
     ALLOCATE(Coefgs(neqnsgs,5),JCoefgs(5),Rhsgs(neqnsgs))
!

! poisson solver (for NSPCG use)
     neqns = (Iend-Ibeg+1)*(Jend-Jbeg+1)*(Kend-Kbeg+1)
     ALLOCATE(Coef(5*neqns,5*15),JCoef(5*15),Rhs(neqns))
!
     end subroutine allocate_variables
!
!--------------------------------------------------------------------------------------------------------------------------     
!
     subroutine generate_grid
!
!--------------------------------------------------------------------------------------------------------------------------
!
!    (5) generate_grid
!
!    This subroutine is used to generate grids
!
!    Called by: main
!
!    Last update: 20/12/2010, Gangfeng Ma
!--------------------------------------------------------------------------------------------------------------------------
!
     use global
     implicit none
     integer :: i,j,k
!
!    horizontal grid
!
     x(Ibeg) = npx*(Mloc-2*Nghost)*dx

     do i = Ibeg+1,Mloc1
     x(i) = x(i-1)+dx
     xc(i-1) = x(i-1)+0.5*dx
     enddo

     do i = Ibeg-1,Ibeg-Nghost,-1
     x(i) = x(i+1)-dx
     xc(i) = x(i+1)-0.5*dx
     enddo

     y(Jbeg) = npy*(Nloc-2*Nghost)*dy

     do j = Jbeg+1,Nloc1
     y(j) = y(j-1)+dy
     yc(j-1) = y(j-1)+0.5*dy
     enddo

     do j = Jbeg-1,Jbeg-Nghost,-1
     y(j) = y(j+1)-dy
     yc(j) = y(j+1)-0.5*dy
     enddo
!
!    vertical grid
!
     if (Ivgrd==1) then
     do k = 1,Kloc
     dsig(k) = 1.0/float(Kglob)
     enddo
     elseif(Ivgrd==2) then
     dsig(Kbeg) = (Grd_R-1.0)/(Grd_R**float(Kglob)-1.0)
     do k = Kbeg+1,Kend
     dsig(k) = dsig(k-1)*Grd_R
     enddo

     do k = 1,Nghost
     dsig(Kbeg-k) = dsig(Kbeg+k-1)
     enddo

     do k = 1,Nghost
     dsig(Kend+k) = dsig(Kend-k+1)
     enddo
     endif

     sig(Kbeg) = Zero
     do k = Kbeg+1,Kloc1
     sig(k) = sig(k-1)+dsig(k-1)
     sigc(k-1) = sig(k-1)+0.5*dsig(k-1)
     enddo
     do k = Kbeg-1,1,-1
     sig(k) = sig(k+1)-dsig(k)
     sigc(k) = sig(k+1)-0.5*dsig(k)
     enddo
!
     end subroutine generate_grid
!
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!  
     subroutine read_bathymetry
!
!------------------------------------------------------------------------------------------
!
!    (6) read_bathymetry
! 
!    This subroutine is used to read bathymetry 
!                                               
!    Called by: main  
!  
!    Gangfeng Ma, 21/12/2010   
!
!    Needs revision to make input depth grid filename arbitrary
!   
!-------------------------------------------------------------------------------------------
!
     use global
     implicit none
     integer :: i,j,m,n,iter,iglob,jglob
     integer :: Maskp(Mglob+1,Nglob+1)
     real(SP), dimension(Mglob+1,Nglob+1) :: HG
! added by Cheng for fluid slide
     real(SP), dimension(Mglob,Nglob) :: Dgstmp
     Hc = Zero ! added by Cheng for initialization
!
!    read bathymetry at grid points
!
     if(trim(adjustl(DEPTH_TYPE))=='CELL_GRID') then
       if(ANA_BATHY) then
         do j = 1,Nglob+1
         do i = 1,Mglob+1
           HG(i,j) = 30.0
         enddo
         enddo
       else
         open(5,file='depth.txt',status='old')
!
!    Need to make file name arbitrary  (kirby, 6/27/16)
!
         do j = 1,Nglob+1
           read(5,*) (HG(i,j),i=1,Mglob+1)
         enddo
       endif

! find permanent dry points (how would this be possible?)
       Maskp = 1
       do j = 1,Nglob+1
       do i = 1,Mglob+1
         if(HG(i,j)<-1000.0) Maskp(i,j) = 0
       enddo
       enddo
 
       ! interpolate depth into cell center
       do j = 1,Nglob
       do i = 1,Mglob
         HCG(i,j) = (HG(i,j)*Maskp(i,j)+HG(i+1,j)*Maskp(i+1,j)+  &
             HG(i,j+1)*Maskp(i,j+1)+HG(i+1,j+1)*Maskp(i+1,j+1))/  &
             (Maskp(i,j)+Maskp(i+1,j)+Maskp(i,j+1)+Maskp(i+1,j+1)+1.e-16)
       enddo
       enddo

       do j = Jbeg,Jend
       do i = Ibeg,Iend
         iglob = npx*(Mloc-2*Nghost)+i-Nghost
         jglob = npy*(Nloc-2*Nghost)+j-Nghost
         Hc(i,j) = HCG(iglob,jglob)
       enddo
       enddo

     elseif(trim(adjustl(DEPTH_TYPE))=='CELL_CENTER') then
       ! read bathymetry at cell center
       if(ANA_BATHY) then
         do j = 1,Nglob
         do i = 1,Mglob
           HCG(i,j) = 0.5
         enddo
         enddo
       else ! not analytical bathymetry, read from depth file
         open(5,file='depth.txt',status='old')
         do j = 1,Nglob
           read(5,*) (HCG(i,j),i=1,Mglob)
         enddo
       endif
	   
         do j = Jbeg,Jend
         do i = Ibeg,Iend
           iglob = npx*(Mloc-2*Nghost)+i-Nghost
           jglob = npy*(Nloc-2*Nghost)+j-Nghost
           Hc(i,j) = HCG(iglob,jglob)
         enddo
         enddo
     endif

     ! collect data into ghost cells 
     call phi_2D_coll(Hc)

     ! save the initial water depth
     Hc0 = Hc


     ! added by Cheng for fluid slide and deformable slide

     Dgs = zero
     Hgs  = zero
     open(5001,file=TRIM(Slide_File),status='old')
     do j = 1,Nglob
       read(5001,*) (Dgstmp(i,j),i=1,Mglob) !slide shape
     enddo
     close(5001)

     do i = Ibeg,Iend
     do j = Jbeg,Jend
       iglob = npx*(Mloc-2*Nghost)+i-Nghost
       jglob = npy*(Nloc-2*Nghost)+j-Nghost

       Dgs(i,j) = Dgstmp(iglob,jglob) !slide shape
       if(Dgs(i,j)-SLIDE_MINTHICK<=1.e-8) then
         Dgs(i,j) = SLIDE_MINTHICK
	   endif
       Hgs(i,j) = Hc0(i,j)-Dgs(i,j) !chi
     enddo
     enddo
	 
     ! wetting-drying mask for fluid slide
     Maskgs = 1
     do j = 1,Nloc
     do i = 1,Mloc
       if(Dgs(i,j)-SLIDE_MINTHICK<=1.e-8) then
         Maskgs(i,j) = 0
       endif
     enddo
     enddo

! interpolate into grid center
     do j=Jbeg,Jend
     do i=Ibeg,Iend
	   if(Dgs(i,j)-SLIDE_MINTHICK<=1.e-8) then
         Hc(i,j) = Hc0(i,j) 
	   else
	     Hc(i,j) = Hc0(i,j)-Dgs(i,j) 
	   endif
     enddo
     enddo
	 
! second derivatives of depth at cell center
     do j = 1,Nloc
     do i = 2,Mloc-1
       Delx2H0(i,j)=((Hc0(i+1,j)-Hc0(i,j))/dx-(Hc0(i,j)-Hc0(i-1,j))/dx)/dx;
     enddo
     enddo
     do j = 1,Nloc
       Delx2H0(1,j) = Delx2H0(2,j)
       Delx2H0(Mloc,j) = Delx2H0(Mloc-1,j)
     enddo
	 
     do j = 2,Nloc-1
     do i = 1,Mloc
        Dely2H0(i,j) = ((Hc0(i,j+1)-Hc0(i,j))/dy-(Hc0(i,j)-Hc0(i,j-1))/dy)/dy;
     enddo
     enddo
     do i = 1,Mloc
       Dely2H0(i,1) = Dely2H0(i,2)
       Dely2H0(i,Nloc) = Dely2H0(i,Nloc-1)
     enddo

     ! collect data into ghost cells
     call phi_int_exch(Maskgs)
	 call phi_2D_exch(Delx2H0)
	 call phi_2D_exch(Dely2H0)
     call phi_2D_coll(Hc)
     call phi_2D_coll(Dgs) 
     call phi_2D_coll(Hgs)  
	 
     ! reconstruct base depth below the fluid slide at x-y faces
     do j = 1,Nloc
     do i = 2,Mloc
       H0fx(i,j) = 0.5*(Hc0(i-1,j)+Hc0(i,j))
     enddo
     H0fx(1,j) = Hc0(1,j)
     H0fx(Mloc1,j) = Hc0(Mloc,j)
     enddo

     do i = 1,Mloc
     do j = 2,Nloc
       H0fy(i,j) = 0.5*(Hc0(i,j-1)+Hc0(i,j))
     enddo
     H0fy(i,1) = Hc0(i,1)
     H0fy(i,Nloc1) = Hc0(i,Nloc)
     enddo

     ! derivatives of base depth below the fluid slide at cell center
     do j = 1,Nloc
     do i = 1,Mloc
       DelxH0(i,j) = (H0fx(i+1,j)-H0fx(i,j))/dx
       DelyH0(i,j) = (H0fy(i,j+1)-H0fy(i,j))/dy
     enddo
     enddo
	call phi_2D_exch(DelxH0)
	call phi_2D_exch(DelyH0)


     ! find pernament dry cells
     Mask_Struct = 1
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Hc(i,j)<-1000.0) then
         Mask_Struct(i,j) = 0
       endif
     enddo
     enddo

     ! reconstruct depth at x-y faces
     do j = 1,Nloc
     do i = 2,Mloc
       Hfx(i,j) = 0.5*(Hc(i-1,j)+Hc(i,j))
     enddo
     Hfx(1,j) = Hc(1,j)
     Hfx(Mloc1,j) = Hc(Mloc,j)
     enddo
	 !hfx0 is flux of bed, which is different from hfx at wet-dry front (Cheng)
	 Hfx0 = Hfx	

     do i = 1,Mloc
     do j = 2,Nloc
       Hfy(i,j) = 0.5*(Hc(i,j-1)+Hc(i,j))
     enddo
     Hfy(i,1) = Hc(i,1)
     Hfy(i,Nloc1) = Hc(i,Nloc)
     enddo
     Hfy0 = Hfy

     ! derivatives of water depth at cell center
     do j = 1,Nloc
     do i = 1,Mloc
       DelxH(i,j) = (Hfx(i+1,j)-Hfx(i,j))/dx
       DelyH(i,j) = (Hfy(i,j+1)-Hfy(i,j))/dy
     enddo
     enddo

     end subroutine read_bathymetry


     subroutine read_bathymetry_comprehensive
!------------------------------------------------------ 
!    This subroutine is used to read bathymetry                                                
!    Called by  
!       main    
!    update: 21/12/2010, Gangfeng Ma    
!    update: 04/13/2012, Fengyan Shi
!    fyshi make a standard slide application    
!-----------------------------------------------------
     use global
     implicit none
     integer :: i,j,m,n,iter,iglob,jglob,Kslide
     integer :: Maskp(Mglob+1,Nglob+1)
     real(SP), dimension(Mglob+1,Nglob+1) :: HG
!     real(SP), dimension(Mglob,Nglob) :: HCG


     ! read bathymetry at grid points
     if(trim(adjustl(DEPTH_TYPE))=='CELL_GRID') then
       if(ANA_BATHY) then
         do j = 1,Nglob+1
         do i = 1,Mglob+1
           HG(i,j) = 0.35
         enddo
         enddo
       else
         open(5,file=TRIM(Depth_File),status='old')
         do j = 1,Nglob+1
           read(5,*) (HG(i,j),i=1,Mglob+1)
         enddo
       endif

       ! find pernament dry points
       Maskp = 1
       do j = 1,Nglob+1
       do i = 1,Mglob+1
         if(HG(i,j)<-1000.0) Maskp(i,j) = 0
       enddo
       enddo
 
       ! interpolate depth into cell center
       do j = 1,Nglob
       do i = 1,Mglob
         HCG(i,j) = (HG(i,j)*Maskp(i,j)+HG(i+1,j)*Maskp(i+1,j)+  &
             HG(i,j+1)*Maskp(i,j+1)+HG(i+1,j+1)*Maskp(i+1,j+1))/  &
             (Maskp(i,j)+Maskp(i+1,j)+Maskp(i,j+1)+Maskp(i+1,j+1)+1.e-16)
       enddo
       enddo

       do j = Jbeg,Jend
       do i = Ibeg,Iend
         iglob = npx*(Mloc-2*Nghost)+i-Nghost
         jglob = npy*(Nloc-2*Nghost)+j-Nghost
         Hc(i,j) = HCG(iglob,jglob)
       enddo
       enddo

     elseif(trim(adjustl(DEPTH_TYPE))=='CELL_CENTER') then
       ! read bathymetry at cell center

       if(ANA_BATHY) then

 
       else ! not analytical bathymetry, read from depth file

         open(5,file=TRIM(Depth_File),status='old')
         do j = 1,Nglob
           read(5,*) (HCG(i,j),i=1,Mglob)
         enddo


         do j = Jbeg,Jend
         do i = Ibeg,Iend
           iglob = npx*(Mloc-2*Nghost)+i-Nghost
           jglob = npy*(Nloc-2*Nghost)+j-Nghost
           Hc(i,j) = HCG(iglob,jglob)
         enddo
         enddo


     endif  ! end analytical 

     endif  ! end cell center

     ! find pernament dry cells
     Mask_Struct = 1
     do j = Jbeg,Jend
     do i = Ibeg,Iend
       if(Hc(i,j)<-1000.0) then
         Mask_Struct(i,j) = 0
       endif
     enddo
     enddo

     ! collect data into ghost cells
     call phi_2D_coll(Hc)

     ! reconstruct depth at x-y faces

! gfma fixed the bug 12/20/2011
    do j = 1,Nloc
    do i = 2,Mloc
      Hfx(i,j) = 0.5*(Hc(i-1,j)+Hc(i,j))
    enddo
    Hfx(1,j) = Hc(1,j)
    Hfx(Mloc1,j) = Hc(Mloc,j)
    enddo

    do i = 1,Mloc
    do j = 2,Nloc
      Hfy(i,j) = 0.5*(Hc(i,j-1)+Hc(i,j))
    enddo
    Hfy(i,1) = Hc(i,1)
    Hfy(i,Nloc1) = Hc(i,Nloc)
    enddo


     ! derivatives of water depth at cell center
     do j = 1,Nloc
     do i = 1,Mloc
       DelxH(i,j) = (Hfx(i+1,j)-Hfx(i,j))/dx
       DelyH(i,j) = (Hfy(i,j+1)-Hfy(i,j))/dy
     enddo
     enddo


     end subroutine read_bathymetry_comprehensive
	 
	 
	 
!
!------------------------------------------------------------------------------------------
!
     subroutine initial
!
!--------------------------------------------------------------------------------  
!
!    (7) initial
!
!    This subroutine is used to initialize model run 
!
!    Called by: main 
!
!    Last update: 21/12/2010, Gangfeng Ma  
!
!--------------------------------------------------------------------------------
!
     use global
     implicit none
     integer  :: i,j,k,n,m,nmax,iglob,jglob
     real(SP) :: xsol(80),zsol(80),zmax,xmax,xterp,zterp,tmp,zc(Kloc), &
                 utmp1,wtmp1,utmp2,wtmp2,xk(321,16),zk(321,16),  &
                 uk(321,16),wk(321,16)
     real(SP) :: Ufric,Zlev1,Zlev,mud_dens,Zslide(Mloc,Nloc),conc_slide, &
                 alpha0,L0,T,bl,hs0,ls0,ls1,ls2,lsx,lsx1,hslide,  &                   
                 eslide,xt,yt,zt,kb,kw,Slope,Xslide,SlideX1,SlideX2
     real(SP), dimension(Mglob,Nglob) :: EtaG
     real(SP), dimension(Mglob,Nglob,Kglob) :: UG,VG,WG,SaliG

     ! simulation time
     TIME = Zero
     RUN_STEP = 0
     dt = dt_ini     
     Screen_Count = Zero
     Plot_Count = Zero
     Plot_Count_Stat = Zero


     Icount = 0
     
     ! working arrays
     D = Zero
     U = Zero
     V = Zero
     W = Zero
     P = Zero
     Omega = Zero
     DU = Zero
     DV = Zero
     DW = Zero
     D0 = Zero
     Eta0 = Zero
     DU0 = Zero
     DV0 = Zero
     DW0 = Zero
     Uf = Zero
     Vf = Zero
     Wf = Zero
     Rho = Rho0
     
     ! source terms
     SourceC = Zero
     SourceX = Zero
     SourceY = Zero

     ! fluxes
     DxL = Zero
     DxR = Zero
     DyL = Zero
     DyR = Zero
     UxL = Zero
     UxR = Zero
     UyL = Zero
     UyR = Zero
     UzL = Zero
     UzR = Zero
     VxL = Zero
     VxR = Zero
     VyL = Zero
     VyR = Zero
     VzL = Zero
     VzR = Zero
     WxL = Zero
     WxR = Zero
     WyL = Zero
     WyR = Zero
     WzL = Zero
     WzR = Zero
     DUxL = Zero
     DUxR = Zero
     DUyL = Zero
     DUyR = Zero
     DVxL = Zero
     DVxR = Zero
     DVyL = Zero
     DVyR = Zero
     DWxL = Zero
     DWxR = Zero
     DWyL = Zero
     DWyR = Zero
     OzL = Zero
     OzR = Zero
     SxL = Zero
     SxR = Zero
     SxS = Zero
     SyL = Zero
     SyR = Zero
     SyS = Zero
     Ex = Zero
	 ExL = Zero
	 ExR = Zero
     Ey = Zero
     EyL = Zero
	 EyR = Zero
     Fx = Zero
	 FxL = Zero
	 FxR = Zero
     Fy = Zero
     FyL = Zero
	 FyR = Zero
     Fz = Zero
     Gx = Zero
	 GxL = Zero
	 GxR = Zero
     Gy = Zero
	 GyL = Zero
	 GyR = Zero
     Gz = Zero
     Hx = Zero
	 HxL = Zero
	 HxR = Zero
     Hy = Zero
	 HyL = Zero
	 HyR = Zero
     Hz = Zero
     EtaxL = Zero
     EtaxR = Zero
     EtayL = Zero
     EtayR = Zero
     DelxEta = Zero
     DelyEta = Zero
     DeltH = Zero
     DeltHo = Zero
     Delt2H = Zero
     DelxD = Zero
     DelyD = Zero
     DelxU = Zero
     DelyU = Zero
     DelzU = Zero
     DelxV = Zero
     DelyV = Zero
     DelzV = Zero
     DelxW = Zero
     DelyW = Zero
     DelzW = Zero
     DelxDU = Zero
     DelyDU = Zero
     DelxDV = Zero
     DelyDV = Zero
     DelxDW = Zero
     DelyDW = Zero
     DelzO = Zero
     Sponge = One
     Cmu = Visc
     CmuHt = Zero
     CmuVt = Zero
     CmuR  = Zero
     Richf = Zero

     Diffxx = Zero
     Diffxy = Zero
     Diffxz = Zero
     Diffyx = Zero
     Diffyy = Zero
     Diffyz = Zero
     Diffzx = Zero
     Diffzy = Zero
     Diffzz = Zero

     Uin_X0 = Zero
     Vin_X0 = Zero
     Win_X0 = Zero
     Ein_X0 = Zero
     Din_X0 = Zero
     Uin_Xn = Zero
	 Uin_Xni = Zero
     Uin_Xni0 = Zero
     Vin_Xn = Zero
     Win_Xn = Zero
     Ein_Xn = Zero
     Din_Xn = Zero   

     HeightMax = Zero ! added by cheng for recording Hmax
	 
     Setup = Zero
     WaveHeight = Zero
     Umean = Zero
     Vmean = Zero
     Num_Zero_Up = 0
     Emax = -1000.
     Emin = 1000.

     WdU = Zero
     WdV = Zero

     Lag_Umean = Zero
     Lag_Vmean = Zero
     Lag_Wmean = Zero
     Euler_Umean = Zero
     Euler_Vmean = Zero
     Euler_Wmean = Zero


     ExtForceX = Zero
     ExtForceY = Zero

     ! baroclinic terms
     DRhoX = Zero
     DRhoY = Zero







! add by Cheng for fluid slide

! add by Cheng for deformable slide
     do j = 1,Nloc
     do i = 1,Mloc
       if(Maskgs(i,j)==1) then
         Ugs(i,j) = SLIDE_INIU
		 Vgs(i,j) = SLIDE_INIV
		 Wgs(i,j) = SLIDE_INIW
	   else
	     Ugs(i,j) = Zero
		 Vgs(i,j) = Zero
		 Wgs(i,j) = Zero
       endif
       DUgs(i,j) = Ugs(i,j)*Dgs(i,j)
       DVgs(i,j) = Vgs(i,j)*Dgs(i,j)
	   DWgs(i,j) = Wgs(i,j)*Dgs(i,j)
     enddo
     enddo
     Ugs0=Ugs;Vgs0=Vgs;Wgs0=Wgs;DUgs0=DUgs;DVgs0=DVgs;DWgs0=DWgs;Wgs00=Wgs0
	 Wtgs=Zero;Wbgs=Zero;Qtgs=Zero;Qbgs=Zero
     DelxUgs=Zero;DelxVgs=Zero;DelxWgs=Zero
	 DelxDUgs=Zero;DelxDVgs=Zero;DelxDWgs=Zero
	 DelxDgs=Zero;DelxHgs=Zero
	 DelyUgs=Zero;DelyVgs=Zero;DelyWgs=Zero
	 DelyDUgs=Zero;DelyDVgs=Zero;DelyDWgs=Zero
	 DelyDgs=Zero;DelyHgs=Zero
     Kx=Zero;SrcmgsX=Zero;SrcpgsX=Zero;SrctgsX=Zero
     Ky=Zero;SrcmgsY=Zero;SrcpgsY=Zero;SrctgsY=Zero;SrctgsZ=Zero
	 UgsP=Zero;VgsP=Zero;WgsP=Zero;DwDt=Zero
     Cxx=Zero;Cxy=Zero;Cxz=Zero;Cyx=Zero;Cyy=Zero;Cyz=Zero;Czx=Zero;Czy=Zero;Czz=Zero;Czz0=Zero
	 Taxx=Zero;Tayx=Zero;Taxy=Zero;Tayy=Zero;Taxz=Zero;Tayz=Zero;
	 Tbxx=Zero;Tbxy=Zero;Tbxz=Zero;Tbyx=Zero;Tbyy=Zero;Tbyz=Zero;Tbzx=Zero;Tbzy=Zero;Tbzz=Zero
	 Tezz=Zero;Pss=Zero;Qsgs=Zero
     DgsxL=Zero;DgsxR=Zero;HSgsxL=Zero;HSgsxR=Zero
	 UgsxL=Zero;UgsxR=Zero;VgsxL=Zero;VgsxR=Zero;WgsxL=Zero;WgsxR=Zero
	 DUgsxL=Zero;DUgsxR=Zero;DVgsxL=Zero;DVgsxR=Zero;DWgsxL=Zero;DWgsxR=Zero
	 Cxgsx=Zero;SgsxL=Zero;SgsxR=Zero;Egsx=Zero;Fgsx=Zero;Ggsx=Zero;Hgsx=Zero
	 EgsxL=Zero;EgsxR=Zero;FgsxL=Zero;FgsxR=Zero
	 GgsxL=Zero;GgsxR=Zero;HgsxL=Zero;HgsxR=Zero
	 SrcmgsxL=Zero;SrcmgsxR=Zero;SrcpgsxL=Zero;SrcpgsxR=Zero
     DgsyL=Zero;DgsyR=Zero;HSgsyL=Zero;HSgsyR=Zero
	 UgsyL=Zero;UgsyR=Zero;VgsyL=Zero;VgsyR=Zero;WgsyL=Zero;WgsyR=Zero
	 DUgsyL=Zero;DUgsyR=Zero;DVgsyL=Zero;DVgsyR=Zero;DWgsyL=Zero;DWgsyR=Zero
	 Cygsy=Zero;SgsyL=Zero;SgsyR=Zero;Egsy=Zero;Fgsy=Zero;Ggsy=Zero;Hgsy=Zero
	 EgsyL=Zero;EgsyR=Zero;FgsyL=Zero;FgsyR=Zero
	 GgsyL=Zero;GgsyR=Zero;HgsyL=Zero;HgsyR=Zero
	 SrcmgsyL=Zero;SrcmgsyR=Zero;SrcpgsyL=Zero;SrcpgsyR=Zero
	 QbgsC=Zero;WgsC=Zero
	 if(trim(RHEO_OPT)=='VISCOUS') then
	    Rhogs = SLIDE_DENSITY
	 elseif(trim(RHEO_OPT)=='GRANULAR') then
	    Rhogs=Zero
		PhiInt=Zero;PhiBed=Zero
	    do j = Jbeg,Jend
        do i = Ibeg,Iend
!		  if(Mask(i,j)==1) then
	        Rhogs(i,j) = GRAIN_DENSITY*SLIDE_CONC + Rho0*(1.0-SLIDE_CONC)
!		  else
!		    Rhogs(i,j) = GRAIN_DENSITY*SLIDE_CONC
!		  endif
		  if(Mask(i,j)==1) then
	        PhiInt(i,j) = PhiInt_F
			PhiBed(i,j) = PhiBed_F
		  else
	        PhiInt(i,j) = PhiInt_A
			PhiBed(i,j) = PhiBed_A
		  endif
		enddo
		enddo
		call phi_2D_coll(Rhogs)
		call phi_2D_coll(PhiInt)
		call phi_2D_coll(PhiBed)
	 endif
	 Rhogs0=Rhogs
	 call flux_coeff_gs



     Tke = Zero
     Eps = Zero
     DTke = Zero
     DEps = Zero
     DTke0 = Zero
     DEps0 = Zero

     !added by m.derakhti            
     Wsx = Zero
     Wsy = Zero
	 

     ! wave breaking mask
     Brks = 0

     ! pressure boundary
     Bc_Prs = Zero

     ! initial surface elevation (user-specified)
     Eta = Zero


     if(IVturb==20.or.IVturb==30) then
       call RandomU
     endif

     if(INITIAL_EUVW) then

        ! initial condition for eta
        open(21,file='eta0.txt',status='old')
        do j = 1,Nglob
          read(21,*) (EtaG(i,j),i=1,Mglob)
        enddo
      
        ! initial condition for U, V, W
        open(22,file='uvw0.txt',status='old')
        do k = 1,Kglob
        do j = 1,Nglob
          read(22,*) (UG(i,j,k),i=1,Mglob)
        enddo
        enddo

        do k = 1,Kglob
        do j = 1,Nglob
          read(22,*) (VG(i,j,k),i=1,Mglob)
        enddo
        enddo

        do k = 1,Kglob
        do j = 1,Nglob
          read(22,*) (WG(i,j,k),i=1,Mglob)
        enddo
        enddo

       do j = Jbeg,Jend
       do i = Ibeg,Iend
         iglob = npx*(Mloc-2*Nghost)+i-Nghost
         jglob = npy*(Nloc-2*Nghost)+j-Nghost
         Eta(i,j) = EtaG(iglob,jglob)
         do k = Kbeg,Kend
           U(i,j,k) = UG(iglob,jglob,k-Nghost)
           V(i,j,k) = VG(iglob,jglob,k-Nghost)
           W(i,j,k) = WG(iglob,jglob,k-Nghost)
         enddo
       enddo
       enddo

!       ! solitary wave from Tanaka solution
!       open(21,file='soliton.dat')
!       do n = 1,80
!         read(21,*) i,xsol(n),zsol(n),tmp,tmp
!       enddo
!       close(21)
!
!       ! find the peak location                                                   
!       zmax = -1.0e+10
!       do n = 1,80
!         if(zsol(n)>zmax) then
!           zmax = zsol(n)
!           xmax = xsol(n)
!         endif
!       enddo
!
!       ! move the peak to x = 3.0m                                                
!       do n = 1,80
!         xsol(n) = xsol(n)+8.0-xmax
!       enddo
!
!       ! interpolate into computational grid
!       do j = Jbeg,Jend
!       do i = Ibeg,Iend
!         if(xc(i)>xsol(80).and.xc(i)<xsol(1)) cycle
!         do n = 2,80
!           if(xc(i)>=xsol(n-1).and.xc(i)<xsol(n)) then
!             xterp = (xc(i)-xsol(n-1))/(xsol(n)-xsol(n-1))
!             Eta(i,j) = (1.0-xterp)*zsol(n-1)+xterp*zsol(n)
!           endif
!         enddo
!       enddo
!       enddo
!
!       open(22,file='plotuv.dat')
!       do n = 1,16
!       do m = 1,321
!         read(22,*) k,xk(m,n),zk(m,n),tmp,uk(m,n),wk(m,n)
!         xk(m,n) = xk(m,n)+8.0
!       enddo
!       enddo
!       close(22)
!       
!       do j = Jbeg,Jend
!       do k = Kbeg,Kend
!       do i = Ibeg,Iend
!         if(xc(i)>xsol(80).and.xc(i)<xsol(1)) cycle
!         zc(k) = (1.0+Eta(i,j))*sigc(k)-1.0
!
!         do n = 2,16
!         do m = 2,321
!           if(xc(i)>=xk(m-1,n).and.xc(i)<xk(m,n).and.  &
!                    zc(k)>=zk(m,n-1).and.zc(k)<zk(m,n)) then
!             xterp = (xc(i)-xk(m-1,n))/(xk(m,n)-xk(m-1,n))
!             zterp = (zc(k)-zk(m,n-1))/(zk(m,n)-zk(m,n-1)) 
!             utmp1 = (1.0-xterp)*uk(m-1,n-1)+xterp*uk(m,n-1)
!             wtmp1 = (1.0-xterp)*wk(m-1,n-1)+xterp*wk(m,n-1)
!             utmp2 = (1.0-xterp)*uk(m-1,n)+xterp*uk(m,n)
!             wtmp2 = (1.0-xterp)*wk(m-1,n)+xterp*wk(m,n)
!             
!             U(i,j,k) = (1.0-zterp)*utmp1+zterp*utmp2
!             W(i,j,k) = (1.0-zterp)*wtmp1+zterp*wtmp2
!           endif
!         enddo
!         enddo
!       enddo
!       enddo
!       enddo
! 100   continue
     endif    

     ! wetting-drying mask
     ! Mask: 1 - wet; 0 - dry
     ! Mask_Struct: 0 - permanent dry point
     ! Mask9: mask for itself and 8 elements around
     Mask = 1
     do j = 1,Nloc
     do i = 1,Mloc
       if((Eta(i,j)+Hc(i,j))<=MinDep) then
         Mask(i,j) = 0
         Eta(i,j) = MinDep-Hc(i,j)
       else
         Mask(i,j) = 1
       endif
     enddo
     enddo
     Mask = Mask*Mask_Struct
	 
     ! collect data into ghost cells (move from before above loop to after it by Cheng)
     call phi_2D_coll(Eta)
     Eta0 = Eta

     ! collect mask into ghost cells
     call phi_int_exch(Mask)

     do j = Jbeg,Jend
     do i = Ibeg,Iend
      Mask9(i,j) = Mask(i,j)*Mask(i-1,j)*Mask(i+1,j)  &
                *Mask(i+1,j+1)*Mask(i,j+1)*Mask(i-1,j+1) &
                *Mask(i+1,j-1)*Mask(i,j-1)*Mask(i-1,j-1)
     enddo
     enddo

     ! total water depth and flux
     D = max(Hc+Eta, MinDep)

     call vel_bc
     call phi_3D_exch(U)
     call phi_3D_exch(V)
     call phi_3D_exch(W)

     do k = 1,Kloc
     do j = 1,Nloc
     do i = 1,Mloc
       DU(i,j,k) = D(i,j)*U(i,j,k)*Mask(i,j)
       DV(i,j,k) = D(i,j)*V(i,j,k)*Mask(i,j)
       DW(i,j,k) = D(i,j)*W(i,j,k)*Mask(i,j)
     enddo
     enddo
     enddo

     if(VISCOUS_FLOW) then
       ! initial seeding values for turbulence
!       Tke_min = 0.5*(1.4e-3)**2
!       Eps_min = 0.09*Tke_min**2/(0.1*Visc)
!       Tke_min = 1.e-12
!       Eps_min = 0.09*Tke_min**2/(1.e-4*Visc)
       Tke_min = 1.e-9
       Eps_min = 1.e-9
       if (RNG) then
          Cmut_min = 8.5e-2*Tke_min**2/Eps_min
       else
          Cmut_min = 9.0e-2*Tke_min**2/Eps_min
       endif
       do k = 1,Kloc
       do j = 1,Nloc
       do i = 1,Mloc
         Tke(i,j,k) = Tke_min
         Eps(i,j,k) = Eps_min
         CmuHt(i,j,k) = Cmut_min
         CmuVt(i,j,k) = Cmut_min
         DTke(i,j,k) = D(i,j)*Tke(i,j,k)*Mask(i,j)
         DEps(i,j,k) = D(i,j)*Eps(i,j,k)*Mask(i,j)
       enddo
       enddo
       enddo
     endif



    ! added by Cheng for hot start
     if(HOTSTART)then
       call hot_start
     endif

!
!    SSP Runge-Kutta method parameters
!
     if(TIME_ORDER(1:3)=='THI') then
       It_Order = 3
       ALPHA(1) = 0.0
       ALPHA(2) = 3.0/4.0
       ALPHA(3) = 1.0/3.0
       BETA(1) = 1.0
       BETA(2) = 1.0/4.0
       BETA(3) = 2.0/3.0
     elseif(TIME_ORDER(1:3)=='SEC') then
       It_Order = 2
       ALPHA(1) = 0.0
       ALPHA(2) = 1.0/2.0
       BETA(1) = 1.0
       BETA(2) = 1.0/2.0
     else
       It_Order = 1
       ALPHA(1) = 0.0
       BETA(1) = 1.0
     endif
!
!    sponge layer
!
     if(SPONGE_ON) then
       call calculate_sponge
     endif
!    
     end subroutine initial
!
!------------------------------------------------------------------------------------------------------
!
