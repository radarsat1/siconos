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
option(WITH_DOCUMENTATION "Build Documentation. Default = OFF" OFF)
option(WITH_PYTHON_WRAPPER "Build python bindings using swig. Default = ON" ON)
option(WITH_DOXYGEN_WARNINGS "Explore doxygen warnings. Default = OFF" OFF)
option(WITH_DOXY2SWIG "Build swig docstrings from doxygen xml output. Default = OFF." OFF)
option(WITH_SYSTEM_INFO "Verbose mode to get some system/arch details. Default = OFF." OFF)
option(WITH_TESTING "Enable 'make test' target" ON)
option(WITH_GIT "Consider sources are under GIT" ON)
option(WITH_SERIALIZATION "Compilation of serialization functions. Default = OFF" OFF)
option(WITH_GENERATION "Generation of serialization functions with gccxml. Default = OFF" OFF)
option(WITH_CXX "Enable CXX compiler for numerics. Default = ON" ON)
option(WITH_UNSTABLE "Enable this to include all 'unstable' sources. Default=OFF" OFF)
option(BUILD_SHARED_LIBS "Building of shared libraries. Default = ON" ON)
option(DEV_MODE "Compilation flags setup for developpers. Default = OFF" OFF)
option(WITH_BULLET "compilation with Bullet Bindings. Default = OFF" ON)
option(WITH_OCC "compilation with OpenCascade Bindings. Default = OFF" OFF)
option(WITH_MUMPS "Compilation with MUMPS solver. Default = OFF" OFF)
option(WITH_FCLIB "link with fclib when this mode is enable. Default = OFF" OFF)
option(WITH_FREECAD "Use FreeCAD. Default = OFF" OFF)
option(WITH_MECHANISMS "Generation of bindings for Saladyn Mechanisms toolbox. Default = OFF" OFF)
option(WITH_XML "Enable xml files i/o. Default = ON" ON)
option(WITH_DOCKER "Build inside a docker container. Default = OFF" OFF)

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

set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_C_FLAGS_DEBUG "-g -O0")
set(CMAKE_BUILD_TYPE "DEBUG")
