# =========================================================================
# EstimateSigma
# =========================================================================
SET(TOOL_NAME "mlmc_combinelevels_rp")
SET(TOOL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/posti/mlmc/recordpoints/combinelevels/")

MESSAGE(STATUS "Posti: building ${TOOL_NAME} tool.")
INCLUDE_DIRECTORIES( ${TOOL_DIR}) 

FILE(GLOB ${TOOL_NAME}F90 ${TOOL_DIR}/combinelevels_rp.f90
                          ${TOOL_DIR}/combinelevels_rp_input.f90
                          ${TOOL_DIR}/combinelevels_rp_vars.f90
                          ${TOOL_DIR}/combinelevels_rp_output.f90
                          #${TOOL_DIR}/../estimatesigma/estimatesigma_RP_utput.f90
                          ${TOOL_DIR}/../../../recordpoints/visualize/rpdata/*.f90
                          ${TOOL_DIR}/../../../recordpoints/visualize/parameters.f90
                          ${TOOL_DIR}/../../../recordpoints/visualize/rpset/*.f90
                          ${TOOL_DIR}/../../../recordpoints/visualize/output/*.f90
                          ${TOOL_DIR}/../../../recordpoints/visualize/fft/spec_vars.f90
                          ${TOOL_DIR}/../../../recordpoints/visualize/navierstokes/*.f90
                          ${TOOL_DIR}/../../../recordpoints/prepare/spline.f90)

                          

ADD_LIBRARY(${TOOL_NAME}libF90 OBJECT ${${TOOL_NAME}F90})
SET_TARGET_PROPERTIES(${TOOL_NAME}libF90 PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS}) 
ADD_LIBRARY(${TOOL_NAME}lib STATIC $<TARGET_OBJECTS:${TOOL_NAME}libF90>)
ADD_DEPENDENCIES(${TOOL_NAME}libF90 libflexishared)


# Standalone Posti executable
add_exec(              ${TOOL_NAME} ${TOOL_DIR}/combinelevels_rp.f90)
TARGET_LINK_LIBRARIES( ${TOOL_NAME} ${TOOL_NAME}lib libflexishared ${LAPACK_LIBRARIES} ${CGNS_LIBRARY})

ADD_DEPENDENCIES(      ${TOOL_NAME} ${TOOL_NAME}lib libflexishared)
SET_TARGET_PROPERTIES( ${TOOL_NAME} PROPERTIES LINKER_LANGUAGE Fortran)


