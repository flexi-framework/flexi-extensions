#include "flexi.h"

MODULE MOD_MLMC_Input
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE ReadSums
  MODULE PROCEDURE ReadSums
END INTERFACE

INTERFACE ReadSumsBF
  MODULE PROCEDURE ReadSumsBF
END INTERFACE

INTERFACE GetParams
  MODULE PROCEDURE GetParams
END INTERFACE

INTERFACE ReadStateFile
  MODULE PROCEDURE ReadStateFile
END INTERFACE

INTERFACE ReadStateFileBF
  MODULE PROCEDURE ReadStateFileBF
END INTERFACE

PUBLIC::ReadSums
PUBLIC::ReadSumsBF
PUBLIC::GetParams
PUBLIC::ReadStateFile
PUBLIC::ReadStateFileBF
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Retrieves relevant header and dateset parameters from Flexi files and stores them in a type
!===================================================================================================================================
SUBROUTINE GetParams(filename)
USE MOD_Globals
USE MOD_MLMC_Vars
USE MOD_IO_HDF5
USE MOD_HDF5_Input,    ONLY: GetDataSize
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: filename !< input filename
!===================================================================================================================================
CALL OpenDataFile(TRIM(filename),create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL GetDataSize(File_ID,'DG_Solution',nDims,HSize)
CHECKSAFEINT(MAXVAL(HSize),4)
CHECKSAFEINT(MINVAL(HSize),4)
nVal=INT(HSize)
DEALLOCATE(HSize)
CALL CloseDataFile()
END SUBROUTINE

!===================================================================================================================================
!> Open a state file, read the old state and store the information later needed to write a new state.
!===================================================================================================================================
SUBROUTINE ReadStateFile(StateFile,DataSetName,offset)
! MODULES                                                                                                                          !
USE MOD_HDF5_Input,    ONLY: OpenDataFile,CloseDataFile,ReadArray,ReadAttribute
USE MOD_IO_HDF5,       ONLY: File_ID
USE MOD_Swapmesh_Vars, ONLY: nVar_State,NState,nElemsOld,Time_State,UOld
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN) :: StateFile !< State file to be read
CHARACTER(LEN=255),INTENT(IN) :: DataSetName
INTEGER,INTENT(IN)            :: offset
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL OpenDataFile(StateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL ReadArray(TRIM(DataSetName),6,&
               (/nVar_State,NState+1,NState+1,NState+1,nElemsOld,1/),offset,6,RealArray=UOld)
CALL ReadAttribute(File_ID,'Time',1,RealScalar=Time_State)
CALL CloseDataFile()
END SUBROUTINE ReadStateFile

!===================================================================================================================================
!> Open a state file, read the old state and store the information later needed to write a new state.
!===================================================================================================================================
SUBROUTINE ReadStateFileBF(StateFile,DataSetName,offset)
! MODULES                                                                                                                          !
USE MOD_HDF5_Input,    ONLY: OpenDataFile,CloseDataFile,ReadArray
USE MOD_IO_HDF5,       ONLY: File_ID
USE MOD_MLMC_Vars,     ONLY: BodyForces
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN) :: StateFile !< State file to be read
CHARACTER(LEN=255),INTENT(IN) :: DataSetName
INTEGER,INTENT(IN)            :: offset
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL OpenDataFile(StateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL ReadArray(TRIM(DataSetName),2,(/9,1/),offset,1,RealArray=BodyForces)
CALL CloseDataFile()
END SUBROUTINE ReadStateFileBF


SUBROUTINE ReadSums(FileName)
!===================================================================================================================================
! Reads in Solution from given HDF5 State File
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_MLMC_Vars
USE MOD_IO_HDF5
USE MOD_HDF5_Input
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN) :: FileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL OpenDataFile(TRIM(FileName),create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL ReadArray('UFineSum',    5,nVal,0,5,RealArray=UFineSum    )
CALL ReadArray('UCoarseSum',  5,nVal,0,5,RealArray=UCoarseSum  )
CALL ReadArray('UFineSqSum',  5,nVal,0,5,RealArray=UFineSqSum  )
CALL ReadArray('UCoarseSqSum',5,nVal,0,5,RealArray=UCoarseSqSum)
CALL ReadArray('DUSqSum',     5,nVal,0,5,RealArray=DUSqSum     )
CALL ReadAttribute(File_ID,'nSamples',1,IntScalar=nSamples_Sums)
CALL ReadAttribute(File_ID,'Time',1,RealScalar=Time_Sums)
CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar=MeshFile_Sums)
CALL CloseDataFile()
END SUBROUTINE ReadSums

SUBROUTINE ReadSumsBF(FileName)
!===================================================================================================================================
! Reads in Solution from given HDF5 State File
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_MLMC_Vars
USE MOD_IO_HDF5
USE MOD_HDF5_Input
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN) :: FileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL OpenDataFile(TRIM(FileName),create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL ReadArray('BodyForcesFineSum',    1,nValBodyForces,0,1,RealArray=BodyForcesFineSum    )
CALL ReadArray('BodyForcesCoarseSum',  1,nValBodyForces,0,1,RealArray=BodyForcesCoarseSum  )
CALL ReadArray('BodyForcesFineSqSum',  1,nValBodyForces,0,1,RealArray=BodyForcesFineSqSum  )
CALL ReadArray('BodyForcesCoarseSqSum',1,nValBodyForces,0,1,RealArray=BodyForcesCoarseSqSum)
CALL ReadArray('DBodyForcesSqSum',     1,nValBodyForces,0,1,RealArray=DBodyForcesSqSum     )
CALL ReadAttribute(File_ID,'nSamples',1,IntScalar=nSamples_Sums)
CALL CloseDataFile()
END SUBROUTINE ReadSumsBF


END MODULE MOD_MLMC_Input
