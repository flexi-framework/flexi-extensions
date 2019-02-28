# =========================================================================
# EstimateSigma
# =========================================================================
SET(TOOL_NAME "mlmc_estimatesigma_rp")
SET(TOOL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/posti/mlmc/recordpoints/estimatesigma/")

MESSAGE(STATUS "Posti: building ${TOOL_NAME} tool.")
INCLUDE_DIRECTORIES( ${TOOL_DIR}) 

FILE(GLOB ${TOOL_NAME}F90 ${TOOL_DIR}/estimatesigma_rp.f90
                          ${TOOL_DIR}/estimatesigma_rp_vars.f90
                          ${TOOL_DIR}/estimatesigma_rp_output.f90
                          ${TOOL_DIR}/estimatesigma_rp_input.f90
                          ${TOOL_DIR}/../../../recordpoints/visualize/rpdata/*.f90
                          ${TOOL_DIR}/../../../recordpoints/visualize/rpset/*.f90
                          ${TOOL_DIR}/../../../recordpoints/visualize/parameters.f90
                          ${TOOL_DIR}/../../../recordpoints/visualize/fft/*.f90
                          ${TOOL_DIR}/../../../recordpoints/visualize/output/*.f90
                          ${TOOL_DIR}/../../../recordpoints/visualize/interpolation/*.f90
                          ${TOOL_DIR}/../../../recordpoints/visualize/navierstokes/*.f90
                          ${TOOL_DIR}/../../../common/fftw/fftw.f90
                          ${TOOL_DIR}/../../../recordpoints/prepare/spline.f90)

ADD_LIBRARY(${TOOL_NAME}libF90 OBJECT ${${TOOL_NAME}F90})
SET_TARGET_PROPERTIES(${TOOL_NAME}libF90 PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS}) 
ADD_LIBRARY(${TOOL_NAME}lib STATIC $<TARGET_OBJECTS:${TOOL_NAME}libF90>)
ADD_DEPENDENCIES(${TOOL_NAME}libF90 flexilibF90 HDF5 fftw)
TARGET_INCLUDE_DIRECTORIES(${TOOL_NAME}libF90 PRIVATE ${FFTW_INCLUDE_DIR})


# Standalone Posti executable
add_exec(              ${TOOL_NAME} ${TOOL_DIR}/estimatesigma_rp.f90)
TARGET_LINK_LIBRARIES( ${TOOL_NAME} ${TOOL_NAME}lib flexilib userblocklib_dummy ${linkedlibs} ${LAPACK_LIBRARIES} ${CGNS_LIBRARY} ${FFTW_LIBRARIES})

ADD_DEPENDENCIES(      ${TOOL_NAME} ${TOOL_NAME}lib flexilib userblocklib_dummy ${linkedlibs} HDF5 fftw)
SET_TARGET_PROPERTIES( ${TOOL_NAME} PROPERTIES LINKER_LANGUAGE Fortran)
