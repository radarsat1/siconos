#
# Some convenience macros
#

# Basic list manipulation
MACRO(CAR var)
  SET(${var} ${ARGV1})
ENDMACRO(CAR)

MACRO(CDR var junk)
  SET(${var} ${ARGN})
ENDMACRO(CDR)

# LIST(APPEND ...) is not correct on <COMPILER>_FLAGS 
MACRO(APPEND_FLAGS)
  CAR(_V ${ARGV})
  CDR(_F ${ARGV})
  SET(${_V} "${${_V}} ${_F}")
ENDMACRO(APPEND_FLAGS)

# The use of ADD_DEFINITION results in a warning with Fortran compiler
MACRO(APPEND_C_FLAGS)
  APPEND_FLAGS(CMAKE_C_FLAGS ${ARGV})
ENDMACRO(APPEND_C_FLAGS)

MACRO(APPEND_CXX_FLAGS)
  APPEND_FLAGS(CMAKE_CXX_FLAGS ${ARGV})
ENDMACRO(APPEND_CXX_FLAGS)

MACRO(APPEND_Fortran_FLAGS)
  APPEND_FLAGS(CMAKE_Fortran_FLAGS ${ARGV})
ENDMACRO(APPEND_Fortran_FLAGS)

# Do them once and remember the values for other projects (-> tests)
MACRO(REMEMBER_INCLUDE_DIRECTORIES _DIRS)
  FOREACH(_D ${_DIRS})
    IF(NOT ${PROJECT_NAME}_REMEMBER_INC_${_D})
      SET(${PROJECT_NAME}_REMEMBER_INC_${_D} TRUE)
      LIST(APPEND ${PROJECT_NAME}_INCLUDE_DIRECTORIES ${_D})
      INCLUDE_DIRECTORIES(${_DIRS})
    ENDIF(NOT ${PROJECT_NAME}_REMEMBER_INC_${_D})
  ENDFOREACH(_D ${_DIRS})
ENDMACRO(REMEMBER_INCLUDE_DIRECTORIES _DIRS)

MACRO(REMEMBER_LINK_DIRECTORIES _DIRS)
  FOREACH(_D ${_DIRS})
    IF(NOT ${PROJECT_NAME}_REMEMBER_LINK_${_D})
      SET(${PROJECT_NAME}_REMEMBER_LINK_${_D} TRUE)
      LIST(APPEND ${PROJECT_NAME}_LINK_DIRECTORIES ${_D})
      LINK_DIRECTORIES(${_D})
    ENDIF(NOT ${PROJECT_NAME}_REMEMBER_LINK_${_D})
  ENDFOREACH(_D ${_DIRS})
ENDMACRO(REMEMBER_LINK_DIRECTORIES _DIRS)

MACRO(REMEMBER_LINK_LIBRARIES _LIBS)
  FOREACH(_LIB ${_LIBS})
    IF(NOT ${PROJECT_NAME}_REMEMBER_LINK_LIBRARIES_${_LIB})
      SET(${PROJECT_NAME}_REMEMBER_LINK_LIBRARIES_${_LIB} TRUE)
      LIST(APPEND ${PROJECT_NAME}_LINK_LIBRARIES ${_LIB})
    ENDIF(NOT ${PROJECT_NAME}_REMEMBER_LINK_LIBRARIES_${_LIB})
  ENDFOREACH(_LIB ${_LIBS})
ENDMACRO(REMEMBER_LINK_LIBRARIES _LIBS)

# link/include is done in the BLAS, LAPACK macros, but not in others.
MACRO(COMPILE_WITH)
  CAR(_NAME ${ARGV})
  CDR(_REST ${ARGV})
  CAR(_REQ ${_REST})
  IF(_REST)
    CDR(_RREST ${_REST})
    CAR(_COMP ${_RREST})
  ENDIF(_REST)
  STRING(TOUPPER ${_NAME} _UNAME)
  LIST(APPEND _NAMES ${_NAME})
  LIST(APPEND _NAMES ${_UNAME})
  SET(_FOUND)
  IF(_REQ STREQUAL STANDARD)
  ELSE(_REQ STREQUAL STANDARD)
    IF(_REQ STREQUAL REQUIRED)
      FIND_PACKAGE(${_NAME} REQUIRED)
    ELSE(_REQ STREQUAL REQUIRED)
      FIND_PACKAGE(${ARGV})
    ENDIF(_REQ STREQUAL REQUIRED)
  ENDIF(_REQ STREQUAL STANDARD)
  FOREACH(_N ${_NAMES})
    IF(${_N}_FOUND)
      SET(_FOUND TRUE)
      IF(${_N}_INCLUDE_DIRS)
		REMEMBER_INCLUDE_DIRECTORIES("${${_N}_INCLUDE_DIRS}")
      ENDIF(${_N}_INCLUDE_DIRS)
      IF(${_N}_INCLUDE_DIR)
        REMEMBER_INCLUDE_DIRECTORIES("${${_N}_INCLUDE_DIR}")
      ENDIF(${_N}_INCLUDE_DIR)
      IF(${_N}_INCLUDE_PATH)
        REMEMBER_INCLUDE_DIRECTORIES("${${_N}_INCLUDE_PATH}")
      ENDIF(${_N}_INCLUDE_PATH)
      IF(${_N}_LIBRARY_DIRS)
        REMEMBER_LINK_DIRECTORIES("${${_N}_LIBRARY_DIRS}")
      ENDIF(${_N}_LIBRARY_DIRS)
      IF(${_N}_LIBRARIES)
        REMEMBER_LINK_LIBRARIES("${${_N}_LIBRARIES}")
      ENDIF(${_N}_LIBRARIES)
      IF(_COMP STREQUAL COMPLETE)
        IF(COMPLETE_${_N}_LIBRARIES)	
          REMEMBER_LINK_LIBRARIES("${COMPLETE_${_N}_LIBRARIES}")
        ENDIF(COMPLETE_${_N}_LIBRARIES)
      ENDIF(_COMP STREQUAL COMPLETE)
      IF(${_N}_DEFINITIONS) # not Fortran is supposed
        APPEND_C_FLAGS(${${_N}_DEFINITIONS})
        APPEND_CXX_FLAGS(${${_N}_DEFINITIONS})
      ENDIF(${_N}_DEFINITIONS)
    ENDIF(${_N}_FOUND)
  ENDFOREACH(_N ${_NAME} ${_UNAME})

  IF(_REQ STREQUAL REQUIRED)
    IF(_FOUND)
    ELSE(_FOUND)
      MESSAGE(FATAL_ERROR "${_NAME} NOT FOUND")
    ENDIF(_FOUND)
  ENDIF(_REQ STREQUAL REQUIRED)

  # update NumericsConfig.h/KernelConfig.h
  IF(NOT CONFIG_H_${_NAME}_CONFIGURED)
    SET(CONFIG_H_${_NAME}_CONFIGURED 1 CACHE BOOL 
      "${PROJECT_SHORT_NAME}Config.h generation for package ${_NAME}")
    CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/config.h.cmake 
      ${CMAKE_BINARY_DIR}/${PROJECT_SHORT_NAME}Config.h)
  ENDIF(NOT CONFIG_H_${_NAME}_CONFIGURED)
  SET(_N)
  SET(_NAME) 
  SET(_UNAME)
  SET(_NAMES)
ENDMACRO(COMPILE_WITH)

# Tests
MACRO(BEGIN_TEST _D)
  SET(_CURRENT_TEST_DIRECTORY ${_D})
  FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${_D})
  
  # find and copy data files : *.mat, *.dat and *.xml, and etc.
  FILE(GLOB_RECURSE _DATA_FILES 
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/${_D}
    *.mat 
    *.dat
    *.xml
    *.DAT
    *.INI)
    
  FOREACH(_F ${_DATA_FILES})
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/${_D}/${_F} ${CMAKE_CURRENT_BINARY_DIR}/${_D}/${_F} COPYONLY)
  ENDFOREACH(_F ${_DATA_FILES})
  
  # configure test CMakeLists.txt (needed for a chdir before running test)
  CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/cmake/CMakeListsForTests.cmake 
    ${CMAKE_CURRENT_BINARY_DIR}/${_CURRENT_TEST_DIRECTORY}/CMakeLists.txt @ONLY)

  SET(_EXE_LIST_${_CURRENT_TEST_DIRECTORY})
ENDMACRO(BEGIN_TEST _D)

# Declaration of a siconos test
MACRO(NEW_TEST)
  CAR(_EXE ${ARGV})
  CDR(_SOURCES ${ARGV})
  LIST(APPEND _EXE_LIST_${_CURRENT_TEST_DIRECTORY} ${_EXE})
  SET(${_EXE}_FSOURCES)
  FOREACH(_F ${_SOURCES})
    LIST(APPEND ${_EXE}_FSOURCES ${CMAKE_CURRENT_SOURCE_DIR}/${_CURRENT_TEST_DIRECTORY}/${_F})
  ENDFOREACH(_F ${_SOURCES})
 
  IF(TEST_MAIN)
    LIST(APPEND ${_EXE}_FSOURCES ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_MAIN})
  ENDIF(TEST_MAIN)
  
  # pb env in ctest, see http://www.vtk.org/Bug/view.php?id=6391#bugnotes
  CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/cmake/ldwrap.c.in 
    ${CMAKE_CURRENT_BINARY_DIR}/${_CURRENT_TEST_DIRECTORY}/${_EXE}.ldwrap.c)
  
ENDMACRO(NEW_TEST)


MACRO(NEW_FC_TEST)
  SET(TEST_SOLVER ${ARGV0})
  SET(TEST_DATA ${ARGV1})
  
  SET(TEST_TOLERANCE ${ARGV2})
  IF(NOT DEFINED TEST_TOLERANCE)
    SET(TEST_TOLERANCE 0)
  ENDIF(NOT DEFINED TEST_TOLERANCE)
  
  SET(TEST_MAXITER ${ARGV3})
  IF(NOT DEFINED TEST_MAXITER)
    SET(TEST_MAXITER 0)
  ENDIF(NOT DEFINED TEST_MAXITER)
  
  SET(TEST_INTERNAL_SOLVER ${ARGV4})
  IF(NOT DEFINED TEST_INTERNAL_SOLVER)
    SET(TEST_INTERNAL_SOLVER 0)
  ENDIF(NOT DEFINED TEST_INTERNAL_SOLVER)
  
  SET(TEST_INTERNAL_SOLVER_TOLERANCE ${ARGV5})
  IF(NOT DEFINED TEST_INTERNAL_SOLVER_TOLERANCE)
    SET(TEST_INTERNAL_SOLVER_TOLERANCE 0)
  ENDIF(NOT DEFINED TEST_INTERNAL_SOLVER_TOLERANCE)
  
  SET(TEST_INTERNAL_SOLVER_MAXITER ${ARGV6})
  IF(NOT DEFINED TEST_INTERNAL_SOLVER_MAXITER)
    SET(TEST_INTERNAL_SOLVER_MAXITER 0)
  ENDIF(NOT DEFINED TEST_INTERNAL_SOLVER_MAXITER)
  
  STRING(REGEX REPLACE SICONOS_FRICTION_ "" TEST_SOLVER_NAME ${TEST_SOLVER})
  STRING(REGEX REPLACE SICONOS_FRICTION_ "" TEST_INTERNAL_SOLVER_NAME ${TEST_INTERNAL_SOLVER})
  STRING(REGEX REPLACE "0" "" TEST_INTERNAL_SOLVER_NAME ${TEST_INTERNAL_SOLVER_NAME})
  STRING(REGEX REPLACE "\\.dat" "" TEST_DATA_NAME ${TEST_DATA})

  SET(TEST_NAME "test-${TEST_SOLVER_NAME}${TEST_INTERNAL_SOLVER_NAME}-${TEST_DATA_NAME}")

  CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/cmake/fctest.c.in 
    ${CMAKE_CURRENT_BINARY_DIR}/${_CURRENT_TEST_DIRECTORY}/${TEST_NAME}.c)
  
  LIST(APPEND _EXE_LIST_${_CURRENT_TEST_DIRECTORY} ${TEST_NAME})
  LIST(APPEND ${TEST_NAME}_FSOURCES ${CMAKE_CURRENT_BINARY_DIR}/${_CURRENT_TEST_DIRECTORY}/${TEST_NAME}.c)

ENDMACRO(NEW_FC_TEST)

MACRO(NEW_LCP_TEST)
  SET(TEST_SOLVER ${ARGV0})
  SET(TEST_DATA ${ARGV1})
  
  SET(TEST_SBM ${ARGV2})
  SET(TEST_SBM_C "SBM")
  IF(NOT DEFINED TEST_SBM)
    SET(TEST_SBM 0)
    SET(TEST_SBM_C "")
  ENDIF(NOT DEFINED TEST_SBM)
    
  STRING(REGEX REPLACE SICONOS_ "" TEST_SOLVER_NAME ${TEST_SOLVER})
  STRING(REGEX REPLACE "\\.dat" "" TEST_DATA_NAME ${TEST_DATA})

  SET(TEST_NAME "test-${TEST_SOLVER_NAME}-${TEST_DATA_NAME}${TEST_SBM_C}")

  CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/cmake/lcptest.c.in 
    ${CMAKE_CURRENT_BINARY_DIR}/${_CURRENT_TEST_DIRECTORY}/${TEST_NAME}.c)
  
  LIST(APPEND _EXE_LIST_${_CURRENT_TEST_DIRECTORY} ${TEST_NAME})
  LIST(APPEND ${TEST_NAME}_FSOURCES ${CMAKE_CURRENT_BINARY_DIR}/${_CURRENT_TEST_DIRECTORY}/${TEST_NAME}.c)

ENDMACRO(NEW_LCP_TEST)

MACRO(NEW_LS_TEST)
  SET(TEST_SOLVER ${ARGV0})
  SET(TEST_DATA ${ARGV1})
  
  SET(TEST_SBM ${ARGV2})
  SET(TEST_SBM_C "SBM")
  IF(NOT DEFINED TEST_SBM)
    SET(TEST_SBM 0)
    SET(TEST_SBM_C "")
  ENDIF(NOT DEFINED TEST_SBM)
    
  STRING(REGEX REPLACE SICONOS_ "" TEST_SOLVER_NAME ${TEST_SOLVER})
  STRING(REGEX REPLACE "\\.dat" "" TEST_DATA_NAME ${TEST_DATA})

  SET(TEST_NAME "test-${TEST_SOLVER_NAME}-${TEST_DATA_NAME}${TEST_SBM_C}")

  CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/cmake/lstest.c.in 
    ${CMAKE_CURRENT_BINARY_DIR}/${_CURRENT_TEST_DIRECTORY}/${TEST_NAME}.c)
  
  LIST(APPEND _EXE_LIST_${_CURRENT_TEST_DIRECTORY} ${TEST_NAME})
  LIST(APPEND ${TEST_NAME}_FSOURCES ${CMAKE_CURRENT_BINARY_DIR}/${_CURRENT_TEST_DIRECTORY}/${TEST_NAME}.c)

ENDMACRO(NEW_LS_TEST)




MACRO(NEW_GMP_TEST)
  SET(TEST_SOLVER ${ARGV0})
  SET(TEST_DATA ${ARGV1})
  SET(TEST_GMP_REDUCED 1)
  
  SET(TEST_TOLERANCE ${ARGV2})
  IF(NOT DEFINED TEST_TOLERANCE)
    SET(TEST_TOLERANCE 0)
  ENDIF(NOT DEFINED TEST_TOLERANCE)
  
  SET(TEST_MAXITER ${ARGV3})
  IF(NOT DEFINED TEST_MAXITER)
    SET(TEST_MAXITER 0)
  ENDIF(NOT DEFINED TEST_MAXITER)
  
  SET(TEST_INTERNAL_SOLVER ${ARGV4})
  IF(NOT DEFINED TEST_INTERNAL_SOLVER)
    SET(TEST_INTERNAL_SOLVER 0)
  ENDIF(NOT DEFINED TEST_INTERNAL_SOLVER)
  
  SET(TEST_INTERNAL_SOLVER_TOLERANCE ${ARGV5})
  IF(NOT DEFINED TEST_INTERNAL_SOLVER_TOLERANCE)
    SET(TEST_INTERNAL_SOLVER_TOLERANCE 0)
  ENDIF(NOT DEFINED TEST_INTERNAL_SOLVER_TOLERANCE)
  
  SET(TEST_INTERNAL_SOLVER_MAXITER ${ARGV6})
  IF(NOT DEFINED TEST_INTERNAL_SOLVER_MAXITER)
    SET(TEST_INTERNAL_SOLVER_MAXITER 0)
  ENDIF(NOT DEFINED TEST_INTERNAL_SOLVER_MAXITER)

  SET(TEST_GMP_REDUCED ${ARGV7})
  IF(NOT DEFINED TEST_GMP_REDUCED)
    SET(TEST_GMP_REDUCED 0)
  ENDIF(NOT DEFINED TEST_GMP_REDUCED)
  
  STRING(REGEX REPLACE SICONOS_FRICTION_ "" TEST_SOLVER_NAME "REDUCED"${TEST_GMP_REDUCED}_ ${TEST_SOLVER})
  STRING(REGEX REPLACE SICONOS_FRICTION_ "" TEST_INTERNAL_SOLVER_NAME ${TEST_INTERNAL_SOLVER})
  STRING(REGEX REPLACE "0" "" TEST_INTERNAL_SOLVER_NAME ${TEST_INTERNAL_SOLVER_NAME})
  STRING(REGEX REPLACE "\\.dat" "" TEST_DATA_NAME ${TEST_DATA})

  SET(TEST_NAME "test-GMP-${TEST_SOLVER_NAME}${TEST_INTERNAL_SOLVER_NAME}-${TEST_DATA_NAME}")

  CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/cmake/gmptest.c.in 
    ${CMAKE_CURRENT_BINARY_DIR}/${_CURRENT_TEST_DIRECTORY}/${TEST_NAME}.c)
  
  LIST(APPEND _EXE_LIST_${_CURRENT_TEST_DIRECTORY} ${TEST_NAME})
  LIST(APPEND ${TEST_NAME}_FSOURCES ${CMAKE_CURRENT_BINARY_DIR}/${_CURRENT_TEST_DIRECTORY}/${TEST_NAME}.c)

ENDMACRO(NEW_GMP_TEST)

# add subdirs (i.e. CMakeLists.txt generated for tests) to the build
MACRO(END_TEST)
  ADD_SUBDIRECTORY(${CMAKE_CURRENT_BINARY_DIR}/${_CURRENT_TEST_DIRECTORY} ${CMAKE_CURRENT_BINARY_DIR}/${_CURRENT_TEST_DIRECTORY})
ENDMACRO(END_TEST)

# to prevent the reference of inside sources directories 
MACRO(CHECK_INSTALL_INCLUDE_DIRECTORIES)
  SET(CHECKED_${PROJECT_NAME}_INCLUDE_DIRECTORIES)
  FOREACH(_D ${${PROJECT_NAME}_INCLUDE_DIRECTORIES})
    IF(_D MATCHES "${CMAKE_SOURCE_DIR}")
    ELSE(_D MATCHES "${CMAKE_SOURCE_DIR}")
      LIST(APPEND CHECKED_${PROJECT_NAME}_INCLUDE_DIRECTORIES ${_D})
    ENDIF(_D MATCHES "${CMAKE_SOURCE_DIR}")
  ENDFOREACH(_D ${${PROJECT_NAME}_INCLUDE_DIRECTORIES})
ENDMACRO(CHECK_INSTALL_INCLUDE_DIRECTORIES)

# debug
MACRO(PRINT_VAR V)
  MESSAGE(STATUS "${V} = ${${V}}")
ENDMACRO(PRINT_VAR V)

# copy directory
MACRO(COPY_DIR SRC DST)
  FILE(GLOB_RECURSE FILES ${SRC} * *.*)
  FOREACH(_F ${FILES})
    MESSAGE(STATUS "COPY_DIR: COPYING ${_F} to ${DST}")
    FILE(RELATIVE_PATH _BF ${SRC} ${_F})
    CONFIGURE_FILE(${_F} ${DST}/${_BF} COPYONLY)
  ENDFOREACH(_F ${FILES})
ENDMACRO(COPY_DIR SRC DST)

# copy pattern files
MACRO(COPY_PATTERN SRC DST PATTERN)
  FILE(GLOB_RECURSE FILES ${SRC} ${PATTERN})
  FOREACH(_F ${FILES})
    MESSAGE(STATUS "COPY_PATTERN: COPYING ${_F} to ${DST}")
    FILE(RELATIVE_PATH _BF ${SRC} ${_F})
    CONFIGURE_FILE(${_F} ${DST}/${_BF} COPYONLY)
  ENDFOREACH(_F ${FILES})
ENDMACRO(COPY_PATTERN SRC DST PATTERN)

# check avaibility of python modules
# adapted from http://www.cmake.org/pipermail/cmake/2011-January/041666.html
function(find_python_module module)
	string(TOUPPER ${module} module_upper)
	if(NOT PY_${module_upper})
		if(ARGC GREATER 1 AND ARGV1 STREQUAL "REQUIRED")
			set(${module}_FIND_REQUIRED TRUE)
		endif()
		# A module's location is usually a directory, but for binary modules
		# it's a .so file.
		execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c" 
			"import re, ${module}; print re.compile('/__init__.py.*').sub('',${module}.__file__)"
			RESULT_VARIABLE _${module}_status 
			OUTPUT_VARIABLE _${module}_location
			ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
		if(NOT _${module}_status)
			set(PY_${module_upper} ${_${module}_location} CACHE STRING 
				"Location of Python module ${module}")
		endif(NOT _${module}_status)
	endif(NOT PY_${module_upper})
	find_package_handle_standard_args(PY_${module} DEFAULT_MSG PY_${module_upper})
endfunction(find_python_module)

