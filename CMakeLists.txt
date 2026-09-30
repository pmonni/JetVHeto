cmake_minimum_required(VERSION 4.0)
set(LIB_MAJOR_VERSION "3")
set(LIB_MINOR_VERSION "0")
set(LIB_PATCH_VERSION "0")
set(LIB_VERSION_STRING "${LIB_MAJOR_VERSION}.${LIB_MINOR_VERSION}.${LIB_PATCH_VERSION}")
if (CMAKE_VERSION VERSION_LESS 4.0)
  project(RadISH)
  set(CMAKE_MACOSX_RPATH 1)
  set(PROJECT_VERSION_MAJOR "${LIB_MAJOR_VERSION}")
  set(PROJECT_VERSION_MINOR "${LIB_MINOR_VERSION}")
  set(PROJECT_VERSION_PATCH"${LIB_PATCH_VERSION}")
  set(PROJECT_VERSION "${LIB_VERSION_STRING}")
else()
  cmake_policy(SET CMP0048 NEW)
  set(CMAKE_MACOSX_RPATH 1)
  PROJECT(RadISH VERSION "${LIB_VERSION_STRING}")
endif()

enable_language(Fortran)

# CMake modules
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

# Add our local modules to the module path
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# A touch of colour
include(${CMAKE_MODULE_PATH}/SetColours.cmake)

if (CMAKE_VERSION VERSION_LESS 3.0.0)
  message ("${Yellow}It appears you are using a version of cmake older than 3.0.0. Please upgrade it.${ColourReset}")
endif()

####################
# Compiler options #
####################

# Error message if Fortran90 is not supported
if(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
  message(FATAL_ERROR "Fortran compiler does not support F90")
endif(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

# Fortran flags
include(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake)

# C++ flags
set(CMAKE_C_FLAGS "-O2 -Wall")
set(CMAKE_CXX_FLAGS "-O2 -Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -Wall -ggdb")

set(SRC ${CMAKE_SOURCE_DIR}/src)
set(EXTERNAL_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/external)
set(LIB ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(SHARE ${CMAKE_CURRENT_BINARY_DIR}/share)
set(BIN ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(PYTHON ${CMAKE_CURRENT_BINARY_DIR}/python)
set(TMP ${CMAKE_CURRENT_BINARY_DIR}/tmp)
set(INCLUDE ${CMAKE_CURRENT_BINARY_DIR}/include)
set(TMD_PTJ_DATA_DIR ${CMAKE_SOURCE_DIR}/data/TMDs_ptj)
set(JETVHETO_INPUTS_DIR ${CMAKE_SOURCE_DIR}/inputs)
set(JETVHETO_SCRIPTS_DIR ${CMAKE_SOURCE_DIR}/scripts)

######################
# External libraries #
######################

set(RadISH_dependencies chaplin)

# CHAPLIN
message ("${BoldBlue}Chaplin will be downloaded and built before building RadISH.${ColourReset}")
add_subdirectory (chaplin)
if (chaplin_FOUND)
  # if CHAPLIN was found we can point location
  get_filename_component (chaplin_in} PATH)
  set (make_command ${make_commandinlin_DIR})
elseif (TARGET chaplin)
  # assume we'll want to install the dependence to force
  # building it before RadISH and installation
  set (make_command ${make_commandinix_dir})
  set (depends ${depends} chaplin)
else()
  message ("${BoldRed}Chaplin is required to install RadISH, but it was not found and there is no rule to install it.${ColourReset}")
  set (fail True)
endif ()
set(CHAPLIN_PATH ${EXTERNAL_BUILD_DIR})
set(CHAPLIN_LD_FLAGS "-L${CHAPLIN_PATH}/lib -lchaplin")
set(CHAPLINLIBDIR ${CHAPLIN_PATH}/lib)

# LHAPDF
find_program(LHAPDF_CONFIG lhapdf-config)
if (LHAPDF_CONFIG)
  execute_process(
    COMMAND "${LHAPDF_CONFIG}" --incdir
    OUTPUT_VARIABLE LHAPDF_INCLUDEDIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  set(LHAPDF_FFLAGS "-I${LHAPDF_INCLUDEDIR}" CACHE STRING INTERNAL)
  execute_process(
    COMMAND "${LHAPDF_CONFIG}" --libs
    OUTPUT_VARIABLE LHAPDF_FLAGS
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  set(LHAPDF_LD_FLAGS ${LHAPDF_FLAGS} CACHE STRING INTERNAL)
  execute_process(
    COMMAND "${LHAPDF_CONFIG}" --libdir
    OUTPUT_VARIABLE LHAPDFLIBDIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
else ()
  message ("${BoldBlue}LHAPDF not found, will be downloaded and built before building RadISH${ColourReset}")
  set(RadISH_dependencies ${RadISH_dependencies} lhapdf)
  add_subdirectory (lhapdf)
  if (lhapdf_FOUND)
    # if LHAPDF was found we can point location
    get_filename_component (lhapdf_in} PATH)
    set (make_command ${make_commandinlin_DIR})
  elseif (TARGET lhapdf)
    # assume we'll want to install the dependence to force
    # building it before RadISH and installation
    set (make_command ${make_commandinix_dir})
    set (depends ${depends} lhapdf)
  else()
    message ("${BoldRed}LHAPDF is required to install RadISH, but it was not found and there is no rule to install it.${ColourReset}")
  set (fail True)
  endif ()
  set (LHAPDF_PATH ${EXTERNAL_BUILD_DIR})
  set (LHAPDF_INCLUDEDIR ${LHAPDF_PATH}/include)
  set (LHAPDF_FFLAGS "-I${LHAPDF_INCLUDEDIR}")
  set (LHAPDFLIBDIR  ${LHAPDF_PATH}/lib)
  set (LHAPDF_LD_FLAGS "-L${LHAPDF_PATH}/lib -lLHAPDF")
endif(LHAPDF_CONFIG)

# HOPPET
find_program(HOPPET_CONFIG hoppet-config)
if (HOPPET_CONFIG)
  execute_process(
    COMMAND "${HOPPET_CONFIG}" --ldflags
    OUTPUT_VARIABLE HOPPET_LD_FLAGS
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  set(HOPPET_LD_FLAGS ${HOPPET_LD_FLAGS} CACHE STRING INTERNAL)
  execute_process(
    COMMAND "${HOPPET_CONFIG}" --fflags
    OUTPUT_VARIABLE HOPPET_FFLAGS
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  set(HOPPET_FFLAGS ${HOPPET_FFLAGS} CACHE STRING INTERNAL)
  include_directories(SYSTEM ${HOPPET_INCLUDES})
else ()
  message ("${BoldBlue}HOPPET not found, will be downloaded and built before building RadISH${ColourReset}")
  set(RadISH_dependencies ${RadISH_dependencies} hoppet)
  add_subdirectory (hoppet)
  if (hoppet_FOUND)
    # if LHAPDF was found we can point location
    get_filename_component (hoppet_in} PATH)
    set (make_command ${make_commandinlin_DIR})
  elseif (TARGET hoppet)
    # assume we'll want to install the dependence to force
    # building it before RadISH and installation
    set (make_command ${make_commandinix_dir})
    set (depends ${depends} hoppet)
  else()
    message ("${BoldRed}HOPPET is required to install RadISH, but it was not found and there is no rule to install it.${ColourReset}")
  set (fail True)
  endif ()
  set (HOPPET_PATH ${EXTERNAL_BUILD_DIR})
  set (HOPPET_INCLUDEDIR ${HOPPET_PATH}/include/hoppet)
  set (HOPPET_FFLAGS "-I${HOPPET_INCLUDEDIR}")
  set (HOPPET_LD_FLAGS "-L${HOPPET_PATH}/lib -lhoppet_v1")
endif(HOPPET_CONFIG)

# Fortran Flags for compilation
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${HOPPET_FFLAGS}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${LHAPDF_FFLAGS}")

# src files
file(GLOB_RECURSE sources_resummation src/*.f90 src/*.f)
# The main program is compiled only for the executable target below.
list(REMOVE_ITEM sources_resummation "${SRC}/jetvheto.f90")

# Have the .mod files placed in the include folder
set(CMAKE_Fortran_MODULE_DIRECTORY ${INCLUDE})

# Compile the shared sources once for the JetVHeto executable.
add_library(objlib OBJECT ${sources_resummation})
target_include_directories(objlib PRIVATE
    ${INCLUDE}
    ${hoppet_INCLUDE}
    ${TMD_PTJ_DATA_DIR}
)
add_dependencies(objlib ${RadISH_dependencies})

add_executable(jetvheto ${SRC}/jetvheto.f90 $<TARGET_OBJECTS:objlib>)
target_include_directories(jetvheto PRIVATE ${INCLUDE} ${hoppet_INCLUDE})
target_link_libraries(jetvheto ${CHAPLIN_LD_FLAGS} ${LHAPDF_LD_FLAGS} ${HOPPET_LD_FLAGS})
add_dependencies(jetvheto ${RadISH_dependencies})
set_target_properties(jetvheto
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)
add_custom_command(TARGET jetvheto POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${TMD_PTJ_DATA_DIR} ${CMAKE_BINARY_DIR}/data/TMDs_ptj
    COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${JETVHETO_INPUTS_DIR} ${CMAKE_BINARY_DIR}/inputs
    COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${JETVHETO_SCRIPTS_DIR} ${CMAKE_BINARY_DIR}/scripts
    COMMENT "Copying N3LL data, input cards, and figure scripts"
)

################
# Install step #
################

install (TARGETS jetvheto DESTINATION bin)
install (DIRECTORY ${TMD_PTJ_DATA_DIR} DESTINATION data)
