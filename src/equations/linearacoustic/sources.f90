!=================================================================================================================================
! Copyright (c) 2010-2024  Prof. Claus-Dieter Munz
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
!
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
#include "flexi.h"
#include "eos.h"
!==================================================================================================================================
!> This module contains routines necessary for the Equation of State. For the linear scalar advection diffusion equation,
!> no EOS is needed. This module provides some dummy routines due to compatibility issues with the Navier-Stokes equation system.
!==================================================================================================================================
MODULE MOD_AcSources
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

INTEGER,PARAMETER      :: ACMASKSHAPE_RAMP        = 1
INTEGER,PARAMETER      :: ACMASKSHAPE_RAMPYPOS    = 11
INTEGER,PARAMETER      :: ACMASKSHAPE_RAMPYNEG    = 12
INTEGER,PARAMETER      :: ACMASKSHAPE_CYLINDRICAL = 2

INTEGER,PARAMETER      :: ACMASKBASEFLOW_CONSTANT  = 1

INTEGER,PARAMETER      :: PRM_ACSRC_PRESSURETIMEDERIV          = 0
INTEGER,PARAMETER      :: PRM_ACSRC_LAMBVEC                    = 1
INTEGER,PARAMETER      :: PRM_ACSRC_PRESSUREMATDERIV           = 2

INTERFACE InitAcSources
  MODULE PROCEDURE InitAcSources
END INTERFACE

INTERFACE PrepAcSources
  MODULE PROCEDURE PrepAcSources
END INTERFACE

INTERFACE ApplyAcSources
  MODULE PROCEDURE ApplyAcSources
END INTERFACE

INTERFACE FinalizeAcSources
  MODULE PROCEDURE FinalizeAcSources
END INTERFACE


PUBLIC::DefineParametersAcSources, InitAcSources, PrepAcSources, ApplyAcSources, FinalizeAcSources,DefineParametersAcMask
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersAcSources()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
USE MOD_ReadInTools ,ONLY: addStrListEntry
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("AcousticSources")
CALL prms%CreateIntFromStringOption('AcSourceVariable', "Type of acoustic source term. pressuretimederiv (0), lambvector (1),"// &
                                                        " pressurematderiv (2)")
CALL addStrListEntry('AcSourceVariable','pressuretimederiv',      PRM_ACSRC_PRESSURETIMEDERIV)
CALL addStrListEntry('AcSourceVariable','lambvector',             PRM_ACSRC_LAMBVEC)
CALL addStrListEntry('AcSourceVariable','pressurematderiv',       PRM_ACSRC_PRESSUREMATDERIV)
CALL prms%CreateStringOption( 'AcSourceProjectName',   "FLEXI project name indicating the source data base to be used.")
CALL prms%CreateStringOption( 'AcSourcePath',          "Source path indicating the source data base to be used.", "./")
CALL prms%CreateStringOption( 'LambVecAvgFile',        "File containing the mean lamb vector field required to calculate it's fluctuation values.")
CALL prms%CreateLogicalOption( 'AcSourceMask',         "Use mask to only use specific parts of the source data.")
CALL prms%CreateLogicalOption( 'AcSourceTemporalRamp', "Use temporal ramping function to smoothly ramp in source data over time.")
CALL prms%CreateRealOption(    'AcSourceMaskRadius',   "Radius of source mask. Everything outside is leveled to zero.")
CALL prms%CreateRealOption(    'AcSourceMaskWidth',    "Width of ramping for source mask.")
CALL prms%CreateRealOption(    't0p5',                 "Parameter for temporal ramping: Time at which 50% of the source amplitude is reached.")
CALL prms%CreateRealOption(    't1',                   "Parameter for temporal ramping: Time at which 99.9% of the source amplitude is reached.")
CALL prms%CreateRealOption(    'TEndAc',               "End time of the acoustic simulation data (mandatory).")
CALL prms%CreateLogicalOption( 'AcSourceDebugOut',     "Output additional fielddata for debug purpouse.")
END SUBROUTINE DefineParametersAcSources


!==================================================================================================================================
!> Initialize the base flow and its gradients for the linearized Euler equations
!==================================================================================================================================
SUBROUTINE InitAcSources()
! MODULES
USE MOD_Globals
USE MOD_IO_HDF5,           ONLY: File_ID,nDims,HSize
USE MOD_HDF5_input,        ONLY: OpenDataFile,CloseDataFile,ReadArray
USE MOD_HDF5_Input,        ONLY: ISVALIDHDF5FILE,ReadAttribute,GetDataSize
USE MOD_ReadInTools,       ONLY: GETINTFROMSTR,GETREALARRAY,GETSTR,GETINT,GETLOGICAL,GETREAL
USE MOD_Preproc
USE MOD_Globals
USE MOD_Mesh_Vars,         ONLY:nGlobalElems
USE MOD_Equation_Vars,     ONLY:Kappa,StrVarNamesBase
USE MOD_Restart_Vars,      ONLY:RestartTime
USE MOD_AcSources_Vars
USE MOD_Equation_Vars,     ONLY: readAcSourcesFromFile
USE MOD_Mesh_Vars,         ONLY:nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=4)               :: txtBuf
CHARACTER(LEN=255)             :: buffer,FileType,listName,LambVecAvgFile
INTEGER                        :: iFile,iEndFile,iStartFile,stat,ioUnit,iVar,iVar2
INTEGER                        :: nSrcFilesTot
INTEGER                        :: nElems_HDF5,N_HDF5 
CHARACTER(LEN=255),ALLOCATABLE :: VarNamesSrc_HDF5(:),VarNamesFirst(:),SrcFileNamesTmp(:)
REAL                           :: dt,dt_min,dt_max
REAL,ALLOCATABLE               :: Time_acsrcTmp(:)
!==================================================================================================================================
SWRITE(UNIT_StdOut,'(A)')' INIT ACOUSTIC SOURCES...'
AcSourceVariable=GETINTFROMSTR('AcSourceVariable') ! pressuretimederiv (0), lambvector (1), pressurematderiv (2)
SELECT CASE(AcSourceVariable)
CASE(PRM_ACSRC_PRESSURETIMEDERIV)
  nVarSrc=2
  nVarSrcRead=1
  ALLOCATE(VarNamesSrc(nVarSrcRead),mapSrcToVar(nVarSrc),inputMapRead(nVarSrcRead))
  VarNamesSrc=VarNamesSrcPtime
  mapSrcToVar(1)=1
  mapSrcToVar(2)=5
  inputMapRead(1)=2
CASE(PRM_ACSRC_LAMBVEC)
  nVarSrc=PP_dim
  nVarSrcRead=PP_dim
  ALLOCATE(VarNamesSrc(nVarSrcRead),mapSrcToVar(nVarSrc),inputMapRead(nVarSrcRead))
  VarNamesSrc=VarNamesSrcLamb
  mapSrcToVar(1) = 2
  mapSrcToVar(2) = 3
  inputMapRead(1)=1
  inputMapRead(2)=2
#if PP_dim==3
  mapSrcToVar(3) = 4
  inputMapRead(3)=3
#endif
  ! We need the fluctuation of the lamb vector, therefore read in a timeavg file to subtract
  LambVecAvgFile=GETSTR('LambVecAvgFile')
  CALL ReadAcSourceMeanfile(LambVecAvgFile)
CASE(PRM_ACSRC_PRESSUREMATDERIV)
!  nVarSrc=2
  nVarSrc=1
  nVarSrcRead=1
  ALLOCATE(VarNamesSrc(nVarSrcRead),mapSrcToVar(nVarSrc),inputMapRead(nVarSrcRead))
  VarNamesSrc=VarNamesSrcPMat
!  mapSrcToVar(1) = 1
!  mapSrcToVar(2) = 5
!  inputMapRead(1)=2
  mapSrcToVar(1) = 5
  inputMapRead(1)=1
END SELECT

!TODO init equation which calles this init is called before initTimeDisc which reads TEnd
! prevent reading twice, is init timedisc earlyer callabel?
! Read the end time TEnd from ini file
!TEndAc     = GETREAL('TEnd')
TEndAc     = GETREAL('TEndAc')

! source masking
CALL DefineParametersAcMask()
CALL InitAcMask()
! doMaskSource=GETLOGICAL('AcSourceMask','F')
! IF(doMaskSource) THEN
!  rMask    =GETREAL('AcSourceMaskRadius')
!  widthMask=GETREAL('AcSourceMaskWidth')
! END IF

! temporal source ramping
doTemporalRamp=GETLOGICAL('AcSourceTemporalRamp','F')
IF(doTemporalRamp) THEN
  t0p5=GETREAL('t0p5')
  t1  =GETREAL('t1')
END IF

IF(readAcSourcesFromFile) THEN
  NTime=4 ! size of stencil for interpolation
  
  AcSourceProjectName=GETSTR('AcSourceProjectName')
  AcSourcePath=GETSTR('AcSourcePath','./')
  
  IF(MPIROOT) THEN
    ! get list of all source files
    ioUnit=GETFREEUNIT()
    WRITE(txtBuf,'(I4)')myRank
    listName=TRIM('filelist'//TRIM(ADJUSTL(TRIM(txtBuf)))//'.txt')
    buffer="ls " // TRIM(AcSourcePath) // "/" // TRIM(AcSourceProjectName) // "_AcSrc_*.h5" // " > "//listName
    call system(buffer)
    call system(buffer)
    call system(buffer) !do three times, gives more stability on HAWK (in theory only once needed, likely a lustre issue)
    open(ioUnit,FILE=TRIM(listName),action="read")
    iFile = 0
    DO
     READ(ioUnit,FMT='(a)',iostat=stat)
     IF (stat/=0) EXIT
     iFile = iFile+1
    END DO
    nSrcFilesTot = iFile
    SWRITE(UNIT_StdOut,'(a,I0)') "  Total number of source files found: " , nSrcFilesTot
  END IF
#if USE_MPI
  CALL MPI_BCAST(nSrcFilesTot,1,MPI_INTEGER,0,MPI_COMM_FLEXI,iError)
#endif

  ALLOCATE(SrcFileNamesTmp(nSrcFilesTot))
  IF(MPIROOT) THEN
    REWIND(ioUnit)
    DO iFile = 1,nSrcFilesTot
     READ(ioUnit,'(a)') SrcFileNamesTmp(iFile)
    END DO
    CLOSE(ioUnit)
    OPEN(unit=ioUnit, iostat=stat, file=TRIM(listName), status='old')
    IF (stat == 0) CLOSE(ioUnit, status='delete') 
  END IF
#if USE_MPI
  CALL MPI_BCAST(SrcFileNamesTmp,LEN(SrcFileNamesTmp)*nSrcFilesTot,MPI_CHARACTER,0,MPI_COMM_FLEXI,iError)
#endif
   
  
  ! get time span
  ALLOCATE(Time_acsrcTmp(nSrcFilesTot))
  dt_min=HUGE(1.)
  dt_max=TINY(1.)
  iStartFile=1
  iEndFile=nSrcFilesTot
  DO iFile = 1,nSrcFilesTot
    CALL OpenDataFile(SrcFileNamesTmp(iFile),create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
    ! read time
    CALL ReadAttribute(File_ID,'Time',1,RealScalar=Time_acsrcTmp(iFile))
  
    IF(Time_acsrcTmp(iFile).LT.RestartTime-2.*dt_max) THEN
      iStartFile=iFile
      CALL CloseDataFile()
      CYCLE
    END IF
  
    ! check  monotonicity
    IF(iFile.GT.1) THEN
      dt=Time_acsrcTmp(iFile)-Time_acsrcTmp(iFile-1)
      IF(dt.LE.0) THEN
        CALL CollectiveStop(__STAMP__,&
          "Non-monotonic time series in source database!")
      END IF
      dt_min=MIN(dt_min,dt)
      dt_max=MAX(dt_max,dt)
    END IF

    CALL CloseDataFile()

    ! Stop if time is sufficient
    IF(Time_acsrcTmp(iFile).GT.TEndAc+2.*dt_max) THEN
      iEndFile=iFile
      EXIT
    END IF
  END DO

  nSrcFiles=iEndFile-iStartFile+1
  ALLOCATE(Time_acsrc(nSrcFiles),SrcFileNames(nSrcFiles))
  Time_acsrc(1:nSrcFiles)=Time_acsrcTmp(iStartFile:iEndFile)
  SrcFileNames=SrcFileNamesTmp(iStartFile:iEndFile)
  DEALLOCATE(Time_acsrcTmp,SrcFileNamesTmp)

  SWRITE(UNIT_StdOut,'(A,I4,A,F9.4,A,F9.4,A)')'  Using source database from ', nSrcFiles,' files ranging from t= '&
                                               ,Time_acsrc(1),' to t= ',Time_acsrc(nSrcFiles),'.'
  SWRITE(UNIT_StdOut,'(A,ES18.9)')'  Min. sampling step: ', dt_min
  SWRITE(UNIT_StdOut,'(A,ES18.9)')'  Max. sampling step: ', dt_max
  SWRITE(UNIT_StdOut,'(A)')'  Performing sanity checks on source database...'

  ! perform checks on source database and create variable mappings
  DO iFile = 1,nSrcFiles
    IF(.NOT.ISVALIDHDF5FILE(SrcFileNames(iFile))) CYCLE

    CALL OpenDataFile(SrcFileNames(iFile),create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
    CALL ReadAttribute(File_ID,'File_Type',1,StrScalar=FileType)
    IF(TRIM(FileType).NE.TRIM('AcSrc')) THEN
      CALL CollectiveStop(__STAMP__,&
        'Provided files are not of file type "AcSrc"!')
    END IF

    CALL GetDataSize(File_ID,'DG_Solution',nDims,HSize)
    ! check size and N
    nElems_HDF5=INT(HSize(5),4)
    N_HDF5=INT(HSize(2),4)-1
    IF((N_HDF5.NE.PP_N).OR.(nElems_HDF5.NE.nGlobalElems)) THEN
        CALL CollectiveStop(__STAMP__,&
          "Source database does not have same polynomial degree or mesh file!")
    END IF

    ! read variable names
    nVarSrc_HDF5=INT(HSize(1),4)
    ALLOCATE(VarNamesSrc_HDF5(nVarSrc_HDF5))
    CALL ReadAttribute(File_ID,'VarNames',nVarSrc_HDF5,StrArray=VarNamesSrc_HDF5)
    IF(iFile.EQ.1) THEN ! use first file to generate variable name mappings and check if required source variables are available
      ! generate mapping
      ALLOCATE(varMapRead(nVarSrcRead))
      varMapRead=-999
      DO iVar=1,nVarSrcRead
         DO iVar2=1,nVarSrc_HDF5
           IF(TRIM(VarNamesSrc_HDF5(iVar2)).EQ.TRIM(VarNamesSrc(iVar))) THEN
             varMapRead(iVar)=iVar2
           END IF
         END DO
      END DO
      IF(ANY(varMapRead.LE.0)) THEN
        CALL CollectiveStop(__STAMP__,&
          "Source database does not contain necessary source variables!")
      END IF
      ALLOCATE(VarNamesFirst(nVarSrc_HDF5))
      VarNamesFirst=VarNamesSrc_HDF5
      DEALLOCATE(VarNamesSrc_HDF5)
    ELSE                ! make sure mappings remain valid throughout the file list
      CALL GetDataSize(File_ID,'DG_Solution',nDims,HSize)
      nVarSrc_HDF5=INT(HSize(1),4)
      CALL ReadAttribute(File_ID,'VarNames',nVarSrc_HDF5,StrArray=VarNamesSrc_HDF5)
      IF(ANY(VarNamesSrc_HDF5.NE.VarNamesFirst)) THEN
        CALL CollectiveStop(__STAMP__,&
          "Variable name mismatch in source database!")
      END IF
      DEALLOCATE(VarNamesSrc_HDF5)
    END IF

    CALL CloseDataFile()

    IF(MOD(iFile,100).EQ.0) THEN
      SWRITE(UNIT_StdOut,'(A,I4,A,I4,A)')'  ',iFile,' of ',nSrcFiles,' checked...'
    END IF
  END DO


  ! check if time series covers Tstart-Tend
  IF(Time_acsrc(1).GT.RestartTime) THEN
    SWRITE(UNIT_StdOut,'(A,F9.4,A)')'  Warning: Source database begins ',Time_acsrc(1)-RestartTime,' after t=tRestart!'
  END IF
  IF(Time_acsrc(nSrcFiles).LT.TEndAc) THEN
    SWRITE(UNIT_StdOut,'(A,F9.4,A)')'  Warning: Source database ends ',TEndAc-Time_acsrc(nSrcFiles),' before t=tEnd!'
  END IF

  SWRITE(UNIT_StdOut,'(A)')'  done!'

  ! Set up first acoustic sources
  CALL ReadInitAcSources(RestartTime)

END IF


SWRITE(UNIT_StdOut,'(A)')' INIT ACOUSTIC SOURCES DONE!'

END SUBROUTINE InitAcSources



!==================================================================================================================================
!> Prepares time interpolation polynomial of source data for initialization 
!=================================================================================================================================
SUBROUTINE ReadInitAcSources(t)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Basis,          ONLY: BarycentricWeights
uSE MOD_Mesh_Vars,      ONLY: nElems
USE MOD_AcSources_Vars, ONLY: nVarSrc,NTime,Time_acsrc,acSrcBuffer,SrcFileNames,iBuffer,nextGlob,nextBuffer,Stencils,iGlob
USE MOD_AcSources_Vars, ONLY: wBaryAcSrc

IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: t  !< current time
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: iC
INTEGER                  :: iRead(4)
INTEGER                  :: iFile
!==================================================================================================================================
! get 4 source samples closest to t
iC=MINLOC(ABS(Time_acsrc(:)-t),DIM=1)
IF(Time_acsrc(iC).GT.t) THEN
  iRead(:)=(/-2,-1,0,1/)
ELSE
  iRead(:)=(/-1,0,1,2/)
END IF
iRead=iRead+iC
IF(MINVAL(iRead).LT.1) THEN ! crop to available range
  iRead=iRead-MINVAL(iRead)+1
END IF

ALLOCATE(acSrcBuffer(nVarSrc,0:PP_N,0:PP_N,0:PP_NZ,nElems,NTime))

! Initialize time index mapping
nextGlob=MAXVAL(iRead)+1
nextBuffer=1
Stencils(1,:)=(/1,2,3,4/)
Stencils(2,:)=(/2,3,4,1/)
Stencils(3,:)=(/3,4,1,2/)
Stencils(4,:)=(/4,1,2,3/)
ALLOCATE(iBuffer(NTime),iGlob(NTime))
iBuffer(:)=Stencils(1,:)
iGlob(:)=Stencils(1,:)

! Read the source data
DO iFile=1,NTime
  CALL ReadAcSourceFile(SrcFileNames(iRead(iFile)),Time_acsrc(iGlob(iFile)),acSrcBuffer(:,:,:,:,:,iFile))
END DO

! prepare Lagrange interpolation polynomials for initial stencil
ALLOCATE(wBaryAcSrc(NTime))
CALL BarycentricWeights(NTime-1,Time_acsrc(iGlob),wBaryAcSrc)
END SUBROUTINE ReadInitAcSources



!==================================================================================================================================
!> Prepares time interpolation polynomial of source data and reads in source data if necessary
!=================================================================================================================================
SUBROUTINE PrepAcSources(t,dt)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Basis,         ONLY: BarycentricWeights
USE MOD_Mesh_Vars,     ONLY: nElems
USE MOD_AcSources_Vars,ONLY: nextGlob,nextBuffer,iBuffer,iGlob,NTime,Stencils,SrcFileNames,nSrcFiles
USE MOD_AcSources_Vars,ONLY: Time_acsrc,acSrcBuffer,wBaryAcSrc,TendAc
!USE MOD_Timedisc_Vars, ONLY: TEnd
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: t  !< current time
REAL,INTENT(IN)  :: dt !< time step
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: iElem,i,j,k,l
REAL                     :: tCheck
!==================================================================================================================================
tCheck=MIN(t+dt,TendAc)
IF(tCheck.GT.Time_acsrc(iGlob(3))) THEN ! check if new source data must be read
  IF(nextGlob.LE.nSrcFiles) THEN
    ! read in new source file
    CALL ReadAcSourceFile(SrcFileNames(nextGlob),Time_acsrc(nextGlob),acSrcBuffer(:,:,:,:,:,nextBuffer))
    ! update indices
    nextGlob=nextGlob+1                       ! entry of next source file in global list
    nextBuffer=nextBuffer+1                   ! Buffer entry to read into in acSrcBuffer
    IF(nextBuffer.EQ.NTime+1) nextBuffer=1
    iGlob=iGlob+1                             ! Mapping from stencil to current global time index
    iBuffer=Stencils(nextBuffer,:)            ! Mapping from stencil to current buffer
  ELSE
    CALL abort(__STAMP__,&
        'No more acoustic source data left!')
  END IF
  ! prepare Lagrange interpolation polynomials for current stencil
  CALL BarycentricWeights(NTime-1,Time_acsrc(iGlob),wBaryAcSrc)
END IF
END SUBROUTINE PrepAcSources


!==================================================================================================================================
!> Read in a single source file and possibly calculate dependent source variables
!=================================================================================================================================
SUBROUTINE ReadAcSourceFile(FileName,t,acSrcBuffer_out)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,      ONLY: nElems,offsetElem
USE MOD_Mesh_Vars,      ONLY: Elem_xGP
USE MOD_AcSources_Vars, ONLY: nVarSrc,nVarSrc_HDF5,varMapRead,inputMapRead,AcSourceVariable
USE MOD_AcSources_Vars, ONLY: acSrcMean
USE MOD_AcSources_Vars, ONLY: doMaskSource,widthMask,rMask
USE MOD_AcSources_Vars, ONLY: doTemporalRamp,t0p5,t1
USE MOD_BaseFlow_Vars,  ONLY: UBase
USE MOD_IO_HDF5,        ONLY: File_ID,nDims,HSize
USE MOD_HDF5_input,     ONLY: OpenDataFile,CloseDataFile,ReadArray
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: FileName                                               !< filename to load
REAL,INTENT(IN)             :: t                                                      !< buffer time    
REAL,INTENT(OUT)            :: acSrcBuffer_out(nVarSrc,0:PP_N,0:PP_N,0:PP_NZ,nElems)  !< output buffer (one sample)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: tmpBuffer(nVarSrc_HDF5,0:PP_N,0:PP_N,0:PP_NZ,nElems)
INTEGER          :: iElem,i,j,k
REAL             :: Time1,Time2
!==================================================================================================================================
IF(MPIROOT)THEN
  WRITE(UNIT_stdOut,'(a)',ADVANCE='NO')' READ ACOUSTIC SOURCE TERMS FROM HDF5 FILE...'
  GETTIME(Time1)
END IF
CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL ReadArray('DG_Solution',5,(/nVarSrc_HDF5,PP_N+1,PP_N+1,PP_NZ+1,nElems/),OffsetElem,5,RealArray=tmpBuffer)
acSrcBuffer_out(inputMapRead(:),:,:,:,:)=tmpBuffer(varMapRead(:),:,:,:,:)
CALL CloseDataFile()

! Calculate dependent source data
SELECT CASE(AcSourceVariable)
CASE(PRM_ACSRC_PRESSURETIMEDERIV,PRM_ACSRC_PRESSUREMATDERIV)
  acSrcBuffer_out=-acSrcBuffer_out
!  acSrcBuffer_out(1,:,:,:,:)=acSrcBuffer_out(2,:,:,:,:)/UBase(SOSP,:,:,:,:)/UBase(SOSP,:,:,:,:) ! rho'=p'/c² (isentropic case)
CASE(PRM_ACSRC_LAMBVEC)
! subtract mean field
  acSrcBuffer_out(:,:,:,:,:)=acSrcBuffer_out(:,:,:,:,:)-acSrcMean(:,:,:,:,:)  ! S = -L'= -(omega x u - mean(omega x u))
END SELECT

! Spatial mask function ("Gaussian shape filter" Ewert & Schroeder 2003)
IF(doMaskSource)THEN
  DO iElem=1,nElems
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      acSrcBuffer_out(:,i,j,k,iElem)=acSrcBuffer_out(:,i,j,k,iElem)*0.5*(1.-TANH(3.8002*(NORM2(Elem_xGP(1:2,i,j,k,iElem))-rMask)/widthMask))
      !TODO add more generall sourc Masks
      acSrcBuffer_out(:,i,j,k,iElem)=acSrcBuffer_out(:,i,j,k,iElem)*0.5*(1.+TANH(3.8002*(Elem_xGP(1,i,j,k,iElem)-0.5)/0.5))
      IF(Elem_xGP(1,i,j,k,iElem) .LT. 0.2) acSrcBuffer_out(:,i,j,k,iElem) = 0
    END DO; END DO; END DO ! i,j,k
  END DO ! iElem
END IF ! doMaskSource
! Temporal ramping
IF(doTemporalRamp)THEN
  DO iElem=1,nElems
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      acSrcBuffer_out(:,i,j,k,iElem)=acSrcBuffer_out(:,i,j,k,iElem)*0.5*(1.+TANH(3.8002*(t-t0p5)/(2.*(t1-t0p5)))) 
    END DO; END DO; END DO ! i,j,k
  END DO ! iElem
END IF ! doTemporalRamp
IF(MPIRoot)THEN
  GETTIME(Time2)
  WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')'DONE  [',Time2-Time1,'s]'
END IF
END SUBROUTINE ReadAcSourceFile



!==================================================================================================================================
!> Read in a mean source file, needed for lamb vector fluctuations
!=================================================================================================================================
SUBROUTINE ReadAcSourceMeanfile(FileName)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,      ONLY: nElems,offsetElem,nGlobalElems
USE MOD_AcSources_Vars, ONLY: nVarSrcRead,nVarSrc,VarNamesSrc
USE MOD_AcSources_Vars, ONLY: acSrcMean 
USE MOD_IO_HDF5,        ONLY: File_ID,nDims,HSize
USE MOD_HDF5_input,        ONLY: OpenDataFile,CloseDataFile,ReadArray,GetDataProps,ReadAttribute
USE MOD_Interpolation_Vars,ONLY: NodeType
USE MOD_Interpolation,     ONLY: GetVandermonde
USE MOD_ChangeBasis,       ONLY: ChangeBasis3D
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*) :: FileName
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE               :: tmpBuffer(:,:,:,:,:)
INTEGER                        :: N_HDF5,nVar_HDF5,nElems_HDF5,nVarMean_HDF5,iVar,iVar2
INTEGER,ALLOCATABLE            :: readMap(:)
CHARACTER(LEN=255),ALLOCATABLE :: VarNamesMean_HDF5(:)
CHARACTER(LEN=255)             :: NodeType_HDF5
REAL,ALLOCATABLE               :: UTmp(:,:,:,:,:),Vdm_NHDF5_N(:,:)
INTEGER                        :: iElem
INTEGER                        :: N_HDF5Z
!==================================================================================================================================
IF(MPIRoot)THEN
  WRITE(UNIT_StdOut,'(A,A,A)',ADVANCE='NO')'  Read mean file ',TRIM(FileName),'...'
END IF
IF(TRIM(FileName).NE.'none') THEN
  CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  CALL GetDataProps(nVar_HDF5,N_HDF5,nElems_HDF5,NodeType_HDF5,'Mean')
  
  IF((nElems_HDF5.NE.nGlobalElems).OR.(nVar_HDF5.LT.PP_nVar))THEN
    CALL abort(__STAMP__,&
               'Source file does not match solution. Elements,nVar',nElems_HDF5,REAL(nVar_HDF5))
  ENDIF

  ! read variable names
  !nVarMean_HDF5=INT(HSize(1),4)
  ALLOCATE(VarNamesMean_HDF5(nVar_HDF5))
  CALL ReadAttribute(File_ID,'VarNames_Mean',nVar_HDF5,StrArray=VarNamesMean_HDF5)
  ! generate variable mapping
  ALLOCATE(readMap(nVarSrcRead))
  readMap=-999
  DO iVar=1,nVarSrcRead
     DO iVar2=1,nVar_HDF5
       IF(TRIM(VarNamesMean_HDF5(iVar2)).EQ.TRIM(VarNamesSrc(iVar))) THEN
         readMap(iVar)=iVar2
       END IF
     END DO
  END DO
  IF(ANY(readMap.LE.0)) THEN
    CALL CollectiveStop(__STAMP__,&
      "Source database does not contain necessary source variables!")
  END IF
  ALLOCATE(tmpBuffer(nVar_HDF5,0:PP_N,0:PP_N,0:PP_NZ,nElems))
  ! Read in state
  IF((N_HDF5.EQ.PP_N).AND.(TRIM(NodeType_HDF5).EQ.TRIM(NodeType)))THEN
    ! No interpolation needed, read solution directly from file
    CALL ReadArray('Mean',5,(/nVar_HDF5,PP_N+1,PP_N+1,PP_NZ+1,nElems/),OffsetElem,5,RealArray=tmpBuffer)
  ELSE
    ! We need to interpolate the solution to the new computational grid
    SWRITE(UNIT_stdOut,*)'Interpolating source from file with N_HDF5=',N_HDF5,' to N=',PP_N
#if PP_dim==3
    N_HDF5Z=N_HDF5
#else
    N_HDF5Z=0
#endif
    ALLOCATE(Vdm_NHDF5_N(0:N_HDF5,0:PP_N))
    ALLOCATE(UTmp(nVar_HDF5,0:N_HDF5,0:N_HDF5,0:N_HDF5Z,nElems))
    CALL GetVandermonde(N_HDF5,NodeType_HDF5,PP_N,NodeType,Vdm_NHDF5_N,modal=.TRUE.)
    CALL ReadArray('Mean',5,(/nVar_HDF5,N_HDF5+1,N_HDF5+1,N_HDF5Z+1,nElems/),OffsetElem,5,RealArray=UTmp)
    DO iElem=1,nElems
      CALL ChangeBasis3D(nVar_HDF5,N_HDF5,PP_N,Vdm_NHDF5_N,UTmp(:,:,:,:,iElem),tmpBuffer(:,:,:,:,iElem))
    END DO
    DEALLOCATE(UTmp,Vdm_NHDF5_N)
  END IF
  CALL CloseDataFile()

  ALLOCATE(acSrcMean(nVarSrc,0:PP_N,0:PP_N,0:PP_NZ,nElems))

  acSrcMean(:,:,:,:,:)=tmpBuffer(readMap(:),:,:,:,:)
  DEALLOCATE(readMap,VarNamesMean_HDF5,tmpBuffer)
ELSE
  SWRITE(UNIT_StdOut,'(A)',ADVANCE='YES')' no baseflow specified, set acSrcMean to zero!'
  ALLOCATE(acSrcMean(nVarSrc,0:PP_N,0:PP_N,0:PP_NZ,nElems))
  acSrcMean(:,:,:,:,:)=0.
END IF

IF(MPIRoot)THEN
  WRITE(UNIT_StdOut,'(A)',ADVANCE='YES')' done!'
END IF
END SUBROUTINE ReadAcSourceMeanfile



!==================================================================================================================================
!> Apply acoustic source data to time operator
!==================================================================================================================================
SUBROUTINE ApplyAcSources(source,t)
! MODULES
USE MOD_Mesh_Vars,     ONLY:nElems
USE MOD_PreProc
USE MOD_AcSources_Vars,ONLY:NTime,nVarSrc,AcSourceVariable,iGlob,wBaryAcSrc,Time_acsrc,acSrcBuffer,mapSrcToVar,iBuffer
USE MOD_Basis,         ONLY:LagrangeInterpolationPolys
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)  :: source(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems)  !< source array
REAL,INTENT(IN)     :: t                                             !< current time
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: iElem,i,j,k,iVar
REAL                     :: L(NTime)
!==================================================================================================================================
CALL LagrangeInterpolationPolys(t,NTime-1,Time_acsrc(iGlob),wBaryAcSrc,L)

DO iElem=1,nElems
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    DO iVar=1,nVarSrc
      source(mapSrcToVar(iVar),i,j,k,iElem)=source(mapSrcToVar(iVar),i,j,k,iElem) &
                                             + SUM(acSrcBuffer(iVar,i,j,k,iElem,iBuffer(:))*L(:))
    END DO
  END DO; END DO; END DO ! i,j,k
END DO ! iElem
END SUBROUTINE ApplyAcSources


!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersAcMask()
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("AcMask")
CALL prms%CreateLogicalOption('AcMaskLayer',    "Turn on to use AcMask regions for reducing reflections at boundaries.",'.FALSE.')
CALL prms%CreateRealOption(   'damping',        "Damping factor of AcMask. U_t=U_t-damping*(U-U_base) in fully damped "//&
                                                "regions.", '1.')
CALL prms%CreateIntFromStringOption( 'AcMaskShape',    "Set shape of AcMask: (1) ramp : cartesian / vector-aligned, (2) "//&
                                                       " cylindrical", multiple=.TRUE.)
CALL addStrListEntry('AcMaskShape','ramp',       AcMaskSHAPE_RAMP)
CALL addStrListEntry('AcMaskShape','rampypos',   ACMASKSHAPE_RAMPYPOS)
CALL addStrListEntry('AcMaskShape','rampyneg',   ACMASKSHAPE_RAMPYNEG)
CALL addStrListEntry('AcMaskShape','cylindrical',AcMaskSHAPE_CYLINDRICAL)
CALL prms%CreateRealOption(   'AcMaskDistance', "Length of AcMask ramp. The AcMask will have maximum strength at the end "//&
                                                "of the ramp and after that point.", multiple=.TRUE.)
CALL prms%CreateRealArrayOption('xStart',       "Coordinates of start position of AcMask ramp (AcMaskShape=ramp) "//&
                                                "or center (AcMaskShape=cylindrical).", multiple=.TRUE.)
CALL prms%CreateRealArrayOption('AcMaskDir',    "Direction vector of the AcMask ramp (AcMaskShape=ramp)", multiple=.TRUE.)
CALL prms%CreateRealOption(   'AcMaskRadius',   "Radius of the AcMask zone (AcMaskShape=cylindrical)", multiple=.TRUE.)
#if (PP_dim==3)
CALL prms%CreateRealArrayOption('AcMaskAxis',   "Axis vector of cylindrical AcMask (AcMaskShape=cylindrical)", multiple=.TRUE.)
#endif
CALL prms%CreateLogicalOption('AcMaskViz',      "Turn on to write a visualization file of AcMask region and strength.",'.FALSE.')
CALL prms%CreateIntFromStringOption( 'AcMaskBaseFlow', "Type of baseflow to be used for AcMask. (1) constant: fixed state,"//&
                                                "(2) exactfunction: exact function, (3) file: read baseflow file, (4) pruett: "//&
                                                "temporally varying, solution adaptive Pruett baseflow",'1')
CALL addStrListEntry('AcMaskBaseFlow','constant',     AcMaskBASEFLOW_CONSTANT)
CALL prms%CreateIntOption(    'AcMaskRefState', "Index of refstate in ini-file (AcMaskBaseFlow=constant)")
CALL prms%CreateIntOption(    'AcMaskExactFunc',"Index of exactfunction (AcMaskBaseFlow=exactfunction)")
CALL prms%CreateStringOption( 'AcMaskBaseFlowFile',"FLEXI solution (e.g. TimeAvg) file from which baseflow is read.")
CALL prms%CreateRealOption(   'tempFilterWidth',"Temporal filter width used to advance Pruett baseflow in time.)")

CALL prms%CreateLogicalOption(  'SkipElemByPoly'   , "Turn on in order to specify the region to skipp elements", '.FALSE.')
CALL prms%CreateIntOption(      'nSkipZones' , "Defines the number of exisiting skipping zones")
CALL prms%CreateIntOption(      'nVertices'     , "Define number of vertices per Skipping Zone defining the Polygon", MULTIPLE=.TRUE.)
CALL prms%CreateRealArrayOption('SkipZoneVertex' , "Focing Vertex that defines polygon", MULTIPLE=.TRUE.)
END SUBROUTINE DefineParametersAcMask


!==================================================================================================================================
!> \brief Initialize AcMask region (get parameters, allocate arrays).
!>
!> Important parameters are:
!>  - damping: Masking strength
!>  - MaskShape: Either a ramp (1) or a cylinder (2).
!==================================================================================================================================
SUBROUTINE InitAcMask
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ReadInTools
USE MOD_AcSources_Vars
!USE MOD_Exactfunc,    ONLY:ExactFunc
USE MOD_Equation_Vars,ONLY:RefStateCons
USE MOD_Mesh_Vars,    ONLY:Elem_xGP,nElems
USE MOD_Output_Vars,  ONLY:ProjectName
USE MOD_PruettDamping,ONLY:InitPruettDamping
USE MOD_Restart_Vars, ONLY:DoRestart,RestartTime,RestartFile
USE MOD_Equation_Vars,ONLY:IniExactFunc
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem,i,j,k,iZone
INTEGER             :: AcMaskExactFunc,AcMaskRefState,AcBaseFlowType
CHARACTER(LEN=255)  :: BaseFlowFile
LOGICAL             :: validBaseFlowFile,elemInside
!==================================================================================================================================
doMaskSource=GETLOGICAL('AcSourceMask','F')
IF(.NOT.doMaskSource) RETURN
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT ACMASK...'

AcMaskViz=GETLOGICAL('AcMaskViz','.FALSE.')

CALL CalcAcMaskRamp()

! Readin of Baseflow parameters
AcBaseFlowType = GETINTFROMSTR('AcMaskBaseFlow')
SELECT CASE(AcBaseflowType)
CASE(AcMaskBASEFLOW_CONSTANT) ! constant baseflow from refstate
  AcMaskRefState  = GETINT('AcMaskRefState')
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,&
    "Undefined AcMaskBaseFlow!")
END SELECT

! Preparation of the baseflow on each Gauss Point
SWRITE(UNIT_StdOut,'(A)') '  Initialize AcMask Base Flow...'
ALLOCATE(AcMaskBaseFlow(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))
SELECT CASE(AcBaseflowType)
CASE(AcMaskBASEFLOW_CONSTANT) ! constant baseflow from refstate
  DO iElem=1,nElems
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      AcMAskBaseFlow(:,i,j,k,iElem)=RefStateCons(:,AcMaskRefState)
    END DO; END DO; END DO
  END DO
END SELECT


!Select Elements which should be skipped completely (defined by polygone)
SkipElemByPoly = GETLOGICAL('SkipElemByPoly','.FALSE.')
IF(SkipElemByPoly) THEN
  nSkipZones = GETINT('nSkipZones')
  ALLOCATE(nVertices(nSkipZones))
  DO iZone=1,nSkipZones
    nVertices(iZone) = GETINT('nVertices')
  END DO
  ALLOCATE(SkipZoneVertex(nSkipZones,MAXVAL(nVertices(:)),2))
  DO iZone=1,nSkipZones
    DO j=1,nVertices(iZone)
      SkipZoneVertex(iZone,j,:) = GETREALARRAY('SkipZoneVertex',2)
    END DO
  END DO

  nSkippedElems=0.
  DO iElem=1,nElems
    DO iZone=1,nSkipZones
      InnerElem: DO k=0,PP_NZ
        DO j=0,PP_N; DO i=0,PP_N
          CALL PointInPoly(Elem_xGP(1,i,j,k,iElem),Elem_xGP(2,i,j,k,iElem),SkipZoneVertex(iZone,1:nVertices(iZone),1), &
                           SkipZoneVertex(iZone,1:nVertices(iZone),2),nVertices(iZone),elemInside)
          IF(elemInside) THEN
            nSkippedElems=nSkippedElems+1
            EXIT InnerElem
          END IF
        END DO; END DO
      END DO InnerElem
    END DO
  END DO
  
  ALLOCATE(skippedElems(nSkippedElems))
  nSkippedElems=0
  DO iElem=1,nElems
    DO iZone=1,nSkipZones
      InnerElem2: DO k=0,PP_NZ
        DO j=0,PP_N; DO i=0,PP_N
          CALL PointInPoly(Elem_xGP(1,i,j,k,iElem),Elem_xGP(2,i,j,k,iElem),SkipZoneVertex(iZone,1:nVertices(iZone),1), &
                           SkipZoneVertex(iZone,1:nVertices(iZone),2),nVertices(iZone),elemInside)
          IF(elemInside) THEN
            nSkippedElems=nSkippedElems+1
            skippedElems(nSkippedElems)=iElem
            EXIT InnerElem2
          END IF
        END DO; END DO
      END DO InnerElem2
    END DO
  END DO
END IF

SWRITE(UNIT_StdOut,'(A)')' INIT ACMASK DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitAcMask



!==================================================================================================================================
!> \brief Compute AcMask shape and strength at each solution point. Visualize AcMask.
!>
!> First, depending on the shape (linear or cylindrical), the strength  of the shape without the damping factor (x_star) is
!> calculated on the solution points.
!> From this, a mapping is built which contains only the elements with x_star > 0 somewhere, which is used to later apply
!> the AcMask only to regions where it is needed.
!> In this AcMask region, the final strength of the AcMask is built by limiting x_star to [0,1] and mutiply it by the damping.
!> If set in the parameter file, a visualization of the AcMask strength is written as .vtu files.
!> At the end, the AcMask is pre-multiplied by the Jacobian since we need to do this anyway when the AcMask is applied.
!==================================================================================================================================
SUBROUTINE CalcAcMaskRamp()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ReadInTools
USE MOD_AcSources_Vars
USE MOD_Output_Vars       ,ONLY:ProjectName
USE MOD_Mesh_Vars         ,ONLY:Elem_xGP
USE MOD_Interpolation_Vars,ONLY:NodeTypeCL,NodeType
USE MOD_Interpolation     ,ONLY:GetVandermonde
USE MOD_Output_Vars       ,ONLY:NVisu,Vdm_GaussN_NVisu
USE MOD_ChangeBasisByDim  ,ONLY:ChangeBasisVolume
USE MOD_Mesh_Vars         ,ONLY:sJ,nElems
USE MOD_VTK               ,ONLY:WriteDataToVTK
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                                 :: applyAcMask(nElems)
INTEGER                                 :: iElem,iAcMaskElem,i,j,k,iRamp
CHARACTER(LEN=255)                      :: FileString,VarNameAcMask(1)
REAL,DIMENSION(  0:PP_N,0:PP_N,0:PP_NZ) :: sigma, x_star
REAL                                    :: r_vec(PP_dim)
REAL,ALLOCATABLE,TARGET                 :: AcDummy(:,:,:,:)
REAL,ALLOCATABLE,TARGET                 :: AcMaskMat_NVisu(:,:,:,:,:)
REAL,ALLOCATABLE,TARGET                 :: Coords_NVisu(:,:,:,:,:)
REAL,POINTER                            :: AcMaskMat_NVisu_p(:,:,:,:,:)
REAL,POINTER                            :: Coords_NVisu_p(:,:,:,:,:)
INTEGER                                 :: nAcMaskRamps
INTEGER,ALLOCATABLE                     :: AcMaskShape(:)
REAL,ALLOCATABLE                        :: xStart(:,:)                 ! Starting Point for AcMask Ramp
REAL,ALLOCATABLE                        :: SpVec(:,:)                  ! Vector defining the ramp direction
REAL,ALLOCATABLE                        :: SpDistance(:)               ! Distance of the AcMask layer
REAL,ALLOCATABLE                        :: SpRadius(:)                 ! Radius of the cylindrical (3D) / radial (2D) AcMask layer
#if(PP_dim==3)
REAL,ALLOCATABLE                        :: SpAxis(:,:)                 ! Axis of the cylindrical AcMask layer (only 3D)
#endif
!==================================================================================================================================
SWRITE(UNIT_StdOut,'(A)') '  Initialize AcMask Ramping Function...'

! Precalculation of the AcMask strength on the whole domain to determine actual AcMask region

nAcMaskRamps= CountOption('AcMaskShape')
ALLOCATE(AcMaskShape(nAcMaskRamps))
ALLOCATE(SpDistance(nAcMaskRamps))
ALLOCATE(xStart(3,nAcMaskRamps))
ALLOCATE(SpVec(3,nAcMaskRamps))
ALLOCATE(SpRadius(nAcMaskRamps))
#if(PP_dim==3)
ALLOCATE(SpAxis(3,nAcMaskRamps))
#endif

DO iRamp=1,nAcMaskRamps
  ! readin geometrical parameters of the AcMask ramp
  AcMaskShape(iRamp)=GETINTFROMSTR('AcMaskShape')
  ! Readin of the AcMask Ramp thickness
  SpDistance(iRamp) = GETREAL('AcMaskDistance')
  ! start AcMask Ramp at xStart
  xStart(:,iRamp)= GETREALARRAY('xStart',3,'(/0.,0.,0./)')
#if PP_dim==2
    IF(xStart(3,iRamp).NE.0) THEN
      CALL CollectiveStop(__STAMP__,'You are computing in 2D! Please set xStart(3) = 0!')
    END IF
#endif
  ! Readin of geometrical parameters for different AcMask shapes
  SELECT CASE(AcMaskShape(iRamp))
  CASE(AcMaskShape_RAMP,ACMASKSHAPE_RAMPYPOS,ACMASKSHAPE_RAMPYNEG) ! ramp aligned with a vector
    SpVec(:,iRamp)= GETREALARRAY('AcMaskDir',3,'(/1.,0.,0./)')
#if PP_dim==2
    IF(SpVec(3,iRamp).NE.0) THEN
      CALL CollectiveStop(__STAMP__,'You are computing in 2D! Please set SpVec(3) = 0!')
    END IF
#endif
SpVec(:,iRamp)=SpVec(:,iRamp)/SQRT(DOT_PRODUCT(SpVec(:,iRamp),SpVec(:,iRamp))) ! Normalize SpVec
  CASE(AcMaskShape_CYLINDRICAL) ! circular AcMask
    SpRadius(iRamp)=GETREAL('AcMaskRadius')
!    SpRadius(iRamp)=0
#if PP_dim==3
    SpAxis(:,iRamp)=GETREALARRAY('AcMaskAxis',3,'(/0.,0.,1./)')
#endif
  END SELECT
END DO!iRamp

applyAcMask=.FALSE.
DO iElem=1,nElems
  x_star=0.
  DO iRamp=1,nAcMaskRamps
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    SELECT CASE(AcMaskShape(iRamp))
      CASE(AcMaskShape_RAMP) ! ramp aligned with a vector
      x_star(i,j,k) =       SUM((Elem_xGP(1:PP_dim,i,j,k,iElem)-xStart(1:PP_dim,iRamp))*SpVec(1:PP_dim,iRamp))/SpDistance(iRamp)
      CASE(ACMASKSHAPE_RAMPYPOS) ! ramp aligned with a vector
      IF(Elem_xGP(2,i,j,k,iElem).GE.0) THEN
        x_star(i,j,k) =       SUM((Elem_xGP(1:PP_dim,i,j,k,iElem)-xStart(1:PP_dim,iRamp))*SpVec(1:PP_dim,iRamp))/SpDistance(iRamp)
      ELSE
        x_star(i,j,k) = 0.
      END IF
      CASE(ACMASKSHAPE_RAMPYNEG) ! ramp aligned with a vector
      IF(Elem_xGP(2,i,j,k,iElem).LT.0) THEN
        x_star(i,j,k) =       SUM((Elem_xGP(1:PP_dim,i,j,k,iElem)-xStart(1:PP_dim,iRamp))*SpVec(1:PP_dim,iRamp))/SpDistance(iRamp)
      ELSE
        x_star(i,j,k) = 0.
      END IF
      CASE(AcMaskShape_CYLINDRICAL) ! cylindrical AcMask
      r_vec(:) = Elem_xGP(:,i,j,k,iElem)-xStart(1:PP_dim,iRamp)
#if(PP_dim==3)
      r_vec = r_vec  -SUM((Elem_xGP(:,i,j,k,iElem)-xStart(:,iRamp))*SpAxis(:,iRamp))*SpAxis(:,iRamp)
#endif
      x_star(i,j,k) = (SQRT(SUM(r_vec*r_vec))-SpRadius(iRamp))/SpDistance(iRamp)
    END SELECT
  END DO; END DO; END DO
  IF(ANY(x_star.GT.0.)) THEN
    applyAcMask(iElem)=.TRUE.
     CYCLE
  END IF
  END DO !iRamp
END DO !iElem=1,nElems

! Get AcMask count and build AcMask mappings
nAcMaskElems=COUNT(applyAcMask)
ALLOCATE(AcMaskMat(0:PP_N,0:PP_N,0:PP_NZ,nAcMaskElems))
ALLOCATE(AcMaskMap(nAcMaskElems))
iAcMaskElem=0
DO iElem=1,nElems
  IF(applyAcMask(iElem))THEN
    iAcMaskElem = iAcMaskElem+1
    AcMaskMap(iAcMaskElem) = iElem
  END IF
END DO

! Calculate the final AcMask strength in the AcMask region
AcMaskMat=0.
DO iAcMaskElem=1,nAcMaskElems
  iElem=AcMaskMap(iAcMaskElem)
  sigma=1.
  DO iRamp=1,nAcMaskRamps
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    SELECT CASE(AcMaskShape(iRamp))
      CASE(AcMaskShape_RAMP) ! ramp aligned with a vector
      x_star(i,j,k) =       SUM((Elem_xGP(1:PP_dim,i,j,k,iElem)-xStart(1:PP_dim,iRamp))*SpVec(1:PP_dim,iRamp))/SpDistance(iRamp)
      CASE(ACMASKSHAPE_RAMPYPOS) ! ramp aligned with a vector
      IF(Elem_xGP(2,i,j,k,iElem).GE.0) THEN
        x_star(i,j,k) =       SUM((Elem_xGP(1:PP_dim,i,j,k,iElem)-xStart(1:PP_dim,iRamp))*SpVec(1:PP_dim,iRamp))/SpDistance(iRamp)
      ELSE
        x_star(i,j,k) = 0.
      END IF
      CASE(ACMASKSHAPE_RAMPYNEG) ! ramp aligned with a vector
      IF(Elem_xGP(2,i,j,k,iElem).LT.0) THEN
        x_star(i,j,k) =       SUM((Elem_xGP(1:PP_dim,i,j,k,iElem)-xStart(1:PP_dim,iRamp))*SpVec(1:PP_dim,iRamp))/SpDistance(iRamp)
      ELSE
        x_star(i,j,k) = 0.
      END IF
      CASE(AcMaskShape_CYLINDRICAL) ! cylindrical AcMask
      r_vec(:) = Elem_xGP(:,i,j,k,iElem)-xStart(1:PP_dim,iRamp)
#if(PP_dim==3)
      r_vec = r_vec  -SUM((Elem_xGP(:,i,j,k,iElem)-xStart(:,iRamp))*SpAxis(:,iRamp))*SpAxis(:,iRamp)
#endif
      x_star(i,j,k) = (SQRT(SUM(r_vec*r_vec))-SpRadius(iRamp))/SpDistance(iRamp)
    END SELECT
  END DO; END DO; END DO
  ! Limit to [0,1]
  x_star = MAX(0.,x_star)
  x_star = MIN(1.,x_star)
  sigma  = MAX(0.,sigma-6.*x_star**5. + 15.*x_star**4. - 10.*x_star**3.)
  END DO !iRamp
  ! Apply damping factor
  AcMaskMat(:,:,:,iAcMaskElem) = sigma(:,:,:)
END DO !iAcMaskElem=1,nAcMaskElems

!Skipp complete elements from beeing an acoustic source element,edfined by an polygon

DEALLOCATE(AcMaskShape)
DEALLOCATE(SpDistance)
DEALLOCATE(xStart)
DEALLOCATE(SpVec)
DEALLOCATE(SpRadius)

! Visualize the AcMask Ramp - until now only 3D visualization!
IF(AcMaskViz) THEN
  FileString=TRIM(INTSTAMP(TRIM(ProjectName),myRank))//'_AcMaskRamp.vtu'
  ALLOCATE(Coords_NVisu(1:3, 0:NVisu,0:NVisu,0:ZDIM(NVisu),nElems))
  ALLOCATE(AcMaskMat_NVisu(1,0:NVisu,0:NVisu,0:ZDIM(NVisu),nElems))
  ALLOCATE(AcDummy(1,0:PP_N,0:PP_N,0:PP_NZ))
  ! Create coordinates of visualization points
  DO iElem=1,nElems
    CALL ChangeBasisVolume(3,PP_N,NVisu,Vdm_GaussN_NVisu,Elem_xGP(1:3,:,:,:,iElem),Coords_NVisu(1:3,:,:,:,iElem))
  END DO
  ! Interpolate solution onto visu grid
  AcMaskMat_NVisu=1.
  DO iAcMaskElem=1,nAcMaskElems
    iElem=AcMaskMap(iAcMaskElem)
    AcDummy(1,:,:,:)=AcMaskMat(:,:,:,iAcMaskElem)
    CALL ChangeBasisVolume(1,PP_N,NVisu,Vdm_GaussN_NVisu,AcDummy(1:1,:,:,:),AcMaskMat_NVisu(1:1,:,:,:,iElem))
  END DO !AcMaskElem=1,nAcMaskElems
  VarNameAcMask(1)='dAcMask'
  Coords_NVisu_p => Coords_NVisu
  AcMaskMat_NVisu_p => AcMaskMat_NVisu
  CALL WriteDataToVTK(1,NVisu,nElems,VarNameAcMask,Coords_NVisu_p,AcMaskMat_NVisu_p,TRIM(FileString),dim=PP_dim)
  DEALLOCATE(Coords_NVisu)
  DEALLOCATE(AcMaskMat_NVisu)
  DEALLOCATE(AcDummy)
END IF !AcMaskViz

! ! Finally add the contribution of the Jacobian to AcMaskMat (JU_src = (U-UBase)*SpMat)
! DO iAcMaskElem=1,nAcMaskElems
!   iElem=AcMaskMap(iAcMaskElem)
!   DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
!     AcMaskMat(i,j,k,iAcMaskElem) = AcMaskMat(i,j,k,iAcMaskElem)/sJ(i,j,k,iElem,0)
!   END DO; END DO; END DO
! END DO

END SUBROUTINE CalcAcMaskRamp


!----------------------------------------------------------------------------------------------------------------------------------!
! Check if a point is inside a polygon(PolyX,PolyY)
!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE PointInPoly(PointX,PointY,PolyX,PolyY,PolyN,Inside) 
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
! insert modules here
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN)         :: PointX,PointY,PolyX(PolyN),PolyY(PolyN)
INTEGER,INTENT(IN)      :: PolyN
LOGICAL,INTENT(INOUT)   :: Inside
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: i,j,InOrOut
REAL                    :: xj,yj,xi,yi
LOGICAL                 :: ix , iy , jx , jy , EOR
!===================================================================================================================================

! EXCLUSIVE OR STATEMENT FUNCTION.
EOR(ix,iy) = (ix .OR. iy) .AND. .NOT.(ix .AND. iy)

Inside = .FALSE.
InOrOut = -1

DO i=1,PolyN
  xi = PolyX(i) - PointX
  yi = PolyY(i) - PointY

  ! CHECK WHETHER THE POINT IN QUESTION IS AT THIS VERTEX.
  IF ( xi.EQ.0.0 .AND. yi.EQ.0.0 ) THEN
     InOrOut = 0 
     RETURN
  ENDIF

  ! j IS NEXT VERTEX NUMBER OF POLYGON.
  j = 1 + MOD(i,PolyN) !Check if last point in polygon
  xj = PolyX(j) - PointX
  yj = PolyY(j) - PointY

  ! IS THIS LINE OF 0 LENGTH ?
  IF ( xi.EQ.xj .AND. yi.EQ.yj ) CYCLE
  ix = (xi.GE.0.0)
  iy = (yi.GE.0.0)
  jx = (xj.GE.0.0)
  jy = (yj.GE.0.0)

  ! CHECK WHETHER (PointX,PointY) IS ON VERTICAL SIDE OF POLYGON.
  IF ( xi.EQ.0.0 .AND. xj.EQ.0.0 .AND. EOR(iy,jy) ) THEN
    InOrOut = 0
    RETURN
  ENDIF
  ! CHECK WHETHER (PointX,PointY) IS ON HORIZONTAL SIDE OF POLYGON.
  IF ( yi.EQ.0.0 .AND. yj.EQ.0.0 .AND. EOR(ix,jx) ) THEN
    InOrOut = 0
    RETURN
  ENDIF

  ! CHECK WHETHER BOTH ENDS OF THIS SIDE ARE COMPLETELY 1) TO RIGHT
  ! OF, 2) TO LEFT OF, OR 3) BELOW (PointX,PointY).
  IF ( .NOT.((iy .OR. jy) .AND. EOR(ix,jx)) ) CYCLE

  ! DOES THIS SIDE OBVIOUSLY CROSS LINE RISING VERTICALLY FROM (PointX,PointY)
  IF ( .NOT.(iy .AND. jy .AND. EOR(ix,jx)) ) THEN
    IF ( (yi*xj-xi*yj)/(xj-xi).LT.0.0 ) THEN
      CYCLE
    ELSEIF ( (yi*xj-xi*yj)/(xj-xi).EQ.0.0 ) THEN
      InOrOut = 0
      RETURN
    ELSE
      InOrOut = -InOrOut
    ENDIF
  ELSE
    InOrOut = -InOrOut
  ENDIF
END DO
IF (InOrOut .GE. 0) Inside = .TRUE.

END SUBROUTINE PointInPoly


!==================================================================================================================================
!> Finalizes the acoustic source module
!==================================================================================================================================
SUBROUTINE FinalizeAcSources()
! MODULES
USE MOD_AcSources_Vars
IMPLICIT NONE
!==================================================================================================================================
SDEALLOCATE(Time_acsrc)
SDEALLOCATE(SrcFileNames)
SDEALLOCATE(VarNamesSrc)
SDEALLOCATE(varMapRead)
SDEALLOCATE(inputMapRead)
SDEALLOCATE(mapSrcToVar)
SDEALLOCATE(iBuffer)
SDEALLOCATE(iGlob)
SDEALLOCATE(acSrcBuffer)
SDEALLOCATE(acSrcMean)
SDEALLOCATE(wBaryAcSrc)
END SUBROUTINE FinalizeAcSources

END MODULE MOD_AcSources
