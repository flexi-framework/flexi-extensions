# =========================================================================
# Nisp Recordpoints
# =========================================================================
SET(TOOL_NAME "posti_nisp_rp")
SET(TOOL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/posti/nisp/recordpoints/")

MESSAGE(STATUS "Posti: building ${TOOL_NAME} tool.")
INCLUDE_DIRECTORIES( ${TOOL_DIR}) 

FILE(GLOB ${TOOL_NAME}F90 ${TOOL_DIR}/posti_nisp_rp.f90
                          ${TOOL_DIR}/nisp_rp.f90
                          ${TOOL_DIR}/nisp_rp_vars.f90
                          ${TOOL_DIR}/nisp_rp_output.f90
                          ${TOOL_DIR}/../../common/vtkOutput/vtkStructuredOutput.f90
                          ${TOOL_DIR}/../../recordpoints/visualize/visualizerecordpoints.f90
                          ${TOOL_DIR}/../../recordpoints/visualize/rpdata/*.f90
                          ${TOOL_DIR}/../../recordpoints/visualize/rpset/*.f90
                          ${TOOL_DIR}/../../recordpoints/visualize/parameters.f90
                          ${TOOL_DIR}/../../recordpoints/visualize/fft/*.f90
                          ${TOOL_DIR}/../../recordpoints/visualize/output/*.f90
                          ${TOOL_DIR}/../../recordpoints/visualize/interpolation/*.f90
                          ${TOOL_DIR}/../../recordpoints/visualize/navierstokes/*.f90
                          ${TOOL_DIR}/../../common/fftw/fftw.f90
                          ${TOOL_DIR}/../../recordpoints/prepare/spline.f90)

ADD_LIBRARY(${TOOL_NAME}libF90 OBJECT ${${TOOL_NAME}F90})
SET_TARGET_PROPERTIES(${TOOL_NAME}libF90 PROPERTIES COMPILE_FLAGS ${FLEXI_COMPILE_FLAGS}) 
ADD_LIBRARY(${TOOL_NAME}lib STATIC $<TARGET_OBJECTS:${TOOL_NAME}libF90>)
ADD_DEPENDENCIES(${TOOL_NAME}libF90 libflexishared)
TARGET_INCLUDE_DIRECTORIES(${TOOL_NAME}libF90 PRIVATE ${FFTW_INCLUDE_DIR})

#check whether these two lines are correct (why link externals to lib AND exec??)
TARGET_LINK_LIBRARIES(${TOOL_NAME}lib ${FFTW_LIBRARIES})

# Standalone Posti executable
add_exec(              ${TOOL_NAME} ${TOOL_DIR}/nisp_rp.f90)
TARGET_LINK_LIBRARIES( ${TOOL_NAME} ${TOOL_NAME}lib libflexishared)

ADD_DEPENDENCIES(      ${TOOL_NAME} ${TOOL_NAME}lib libflexishared)
#SET_TARGET_PROPERTIES( ${TOOL_NAME} PROPERTIES LINKER_LANGUAGE Fortran)
