#include "flexi.h"

MODULE MOD_EstimateSigma_ReadState
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
INTERFACE ReadOldSums
  MODULE PROCEDURE ReadOldSums
END INTERFACE

PUBLIC::ReadOldSums
!===================================================================================================================================

CONTAINS


SUBROUTINE ReadOldSums()
!===================================================================================================================================
! Reads in Solution from given HDF5 State File
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_EstimateSigma_Vars
USE MOD_Swapmesh_Vars
USE MOD_IO_HDF5
USE MOD_HDF5_Input

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE     :: tmp(:,:,:,:,:)
!===================================================================================================================================
ALLOCATE(tmp(5*nVarTotal,0:NNew,0:NNew,0:NNew,nElemsNew))

CALL OpenDataFile(FileNameSums,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL ReadArray('DG_Solution',5,(/5*nVarTotal, NNew+1, NNew+1, NNew+1, nElemsNew/),0,5,RealArray=tmp )
CALL CloseDataFile()

UFineSum     = tmp(1            :  nVarTotal,:,:,:,:)
UCoarseSum   = tmp(  nVarTotal+1:2*nVarTotal,:,:,:,:)
UFineSqSum   = tmp(2*nVarTotal+1:3*nVarTotal,:,:,:,:)
UCoarseSqSum = tmp(3*nVarTotal+1:4*nVarTotal,:,:,:,:)
DUSqSum      = tmp(4*nVarTotal+1:           ,:,:,:,:)
SDEALLOCATE(tmp)
END SUBROUTINE ReadOldSums 


END MODULE MOD_EstimateSigma_ReadState
