########################################################
# cmake file for building Marlin ACTSTracking Package
# @author Karol Krizka
CMAKE_MINIMUM_REQUIRED(VERSION 3.15 FATAL_ERROR)
########################################################

# project name
PROJECT( ACTSTracking )

# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 1 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )

# load default settings from ILCSOFT_CMAKE_MODULES
#INCLUDE( ilcsoft_default_settings )


FIND_PACKAGE( Marlin 1.0 REQUIRED ) # minimum required Marlin version
#INCLUDE_DIRECTORIES( SYSTEM ${Marlin_INCLUDE_DIRS} )
#LINK_LIBRARIES( ${Marlin_LIBRARIES} )
#ADD_DEFINITIONS( ${Marlin_DEFINITIONS} )

FIND_PACKAGE( ROOT REQUIRED)
#INCLUDE_DIRECTORIES( SYSTEM ${ROOT_INCLUDE_DIRS} )
#LINK_LIBRARIES( ${ROOT_LIBRARIES} )
#ADD_DEFINITIONS( ${ROOT_DEFINITIONS} )

### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ENDIF()

### DATA ###########################################################
ADD_SUBDIRECTORY(example)

### LIBRARY #################################################################

# sources
ADD_LIBRARY(ACTSTracking SHARED
  ./src/ACTSDetectorLoaderProc.cxx
  )
TARGET_INCLUDE_DIRECTORIES(ACTSTracking PRIVATE ./ACTSTracking ${Marlin_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(ACTSTracking ${ROOT_LIBRARIES} ActsCore ActsPluginTGeo ActsPluginJson)
