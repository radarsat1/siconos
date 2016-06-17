# ================================================================
# All the default values for siconos cmake parameters
#
# Usage:
# cmake path-to-sources
#  --> to keep default value
# 
# cmake path-to-sources -DWITH_PYTHON_WRAPPER=ON
#  --> to enable (ON), or disable (OFF) the concerned option.
#
# For details about all these options check siconos install guide.
# ================================================================

set(CMAKE_INSTALL_PREFIX "/home/sinclairs/.local" CACHE PATH "prefix" FORCE)

# --------- User-defined options ---------
# Use cmake -DOPTION_NAME=some-value ... to modify default value.
set(WITH_DOCUMENTATION OFF CACHE BOOL "Build Documentation. Default = OFF" FORCE)
set(WITH_PYTHON_WRAPPER ON CACHE BOOL "Build python bindings using swig. Default = ON" FORCE)
set(WITH_DOXYGEN_WARNINGS OFF CACHE BOOL "Explore doxygen warnings. Default = OFF" FORCE)
set(WITH_DOXY2SWIG OFF CACHE BOOL "Build swig docstrings from doxygen xml output. Default = OFF." FORCE)
set(WITH_SYSTEM_INFO OFF CACHE BOOL "Verbose mode to get some system/arch details. Default = OFF." FORCE)
set(WITH_TESTING OFF CACHE BOOL "Enable 'make test' target" FORCE)
set(WITH_GIT ON CACHE BOOL "Consider sources are under GIT" FORCE)
set(WITH_SERIALIZATION OFF CACHE BOOL "Compilation of serialization functions. Default = OFF" FORCE)
set(WITH_GENERATION OFF CACHE BOOL "Generation of serialization functions with gccxml. Default = OFF" FORCE)
set(WITH_CXX ON CACHE BOOL "Enable CXX compiler for numerics. Default = ON" FORCE)
set(WITH_UNSTABLE OFF CACHE BOOL "Enable this to include all 'unstable' sources. Default=OFF" FORCE)
set(BUILD_SHARED_LIBS ON CACHE BOOL "Building of shared libraries. Default = ON" FORCE)
set(DEV_MODE OFF CACHE BOOL "Compilation flags setup for developpers. Default = OFF" FORCE)
set(WITH_BULLET ON CACHE BOOL "compilation with Bullet Bindings. Default = OFF" FORCE)
set(WITH_OCC OFF CACHE BOOL "compilation with OpenCascade Bindings. Default = OFF" FORCE)
set(WITH_MUMPS OFF CACHE BOOL "Compilation with MUMPS solver. Default = OFF" FORCE)
set(WITH_FCLIB OFF CACHE BOOL "link with fclib when this mode is enable. Default = OFF" FORCE)
set(WITH_FREECAD OFF CACHE BOOL "Use FreeCAD. Default = OFF" FORCE)
set(WITH_MECHANISMS OFF CACHE BOOL "Generation of bindings for Saladyn Mechanisms toolbox. Default = OFF" FORCE)
set(WITH_XML ON CACHE BOOL "Enable xml files i/o. Default = ON" FORCE)
set(WITH_DOCKER OFF CACHE BOOL "Build inside a docker container. Default = OFF" FORCE)

#####################################
# Use ccache if available
find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)

# Set python install mode:
# - user --> behave as 'python setup.py install --user'
# - standard --> install in python site-package (ie behave as python setup.py install)
# - prefix --> install in python CMAKE_INSTALL_PREFIX (ie behave as python setup.py install --prefix=CMAKE_INSTALL_PREFIX)
IF(UNIX)
  # on unix, there is no reason to use the standard option. By default, CMAKE_INSTALL_PREFIX is set to /usr/local and therefore,
  # the python packages should be installed in /usr/local/...
  set(siconos_python_install "prefix" CACHE STRING "Install mode for siconos python package")
ELSE(UNIX)
  set(siconos_python_install "standard" CACHE STRING "Install mode for siconos python package")
ENDIF(UNIX)


# List of components to build and installed
# List of siconos component to be installed
# complete list = externals numerics kernel control mechanics io
set(COMPONENTS externals numerics kernel control mechanics io CACHE INTERNAL "List of siconos components to build and install")

# set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
# set(CMAKE_C_FLAGS_DEBUG "-g -O0")
set(CMAKE_BUILD_TYPE "RELEASE")
set(CMAKE_CXX_FLAGS "-Werror=delete-non-virtual-dtor -Werror=unreachable-code -Wall -Wuninitialized -Wextra -Wno-unused-parameter -Werror=implicit-function-declaration -Werror=logical-not-parentheses -Werror=sizeof-array-argument -Werror=array-bounds -Werror=type-limits -Werror=return-type -Wodr -Werror=overloaded-virtual -Wformat=2 -Werror=format-security -Werror=non-virtual-dtor")
