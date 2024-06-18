cmake_minimum_required(VERSION 3.1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(BFPPlugin)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_FILES_DIRECTORY ${CMAKE_BINARY_DIR}/CMakeFiles)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Set a default build type if none was specified
set(default_build_type "Release")

find_package(PkgConfig)
# Find FEBio SDK or git repo automatically
if(WIN32)
	set(TEMP_PATHS ${CMAKE_SOURCE_DIR}/.. ${CMAKE_SOURCE_DIR}/../.. $ENV{HOMEPATH}/ $ENV{HOMEPATH}/source/repos $ENV{HOMEPATH}/*)
else()
    set(TEMP_PATHS ${CMAKE_SOURCE_DIR}/.. ${CMAKE_SOURCE_DIR}/../.. $ENV{HOME}/ $ENV{HOME}/*)
endif()

find_path(FEBio_SDK FECore/Archive.h
    PATHS ${TEMP_PATHS}
    PATH_SUFFIXES FEBio
    DOC "Path to the FEBio SDK, or git repo.")

if(NOT FEBio_SDK)
    if(WIN32)
        set(TEMP_PATHS $ENV{PROGRAMFILES}/* $ENV{HOMEPATH}/*)
    elseif(APPLE)
        set(TEMP_PATHS /Applications/* $ENV{HOME}/*)
    else()
        set(TEMP_PATHS $ENV{HOME}/*)
    endif()

    find_path(FEBio_SDK "include/FECore/Archive.h"
        PATHS ${TEMP_PATHS}
        PATH_SUFFIXES sdk
        DOC "Path to the FEBio SDK, or git repo.")
endif()

if(NOT FEBio_SDK)
    set(FEBio_SDK "" CACHE PATH "Path to the FEBio SDK, or git repo.")
    message(FATAL_ERROR "Unable to find path to FEBio SDK or git repo automatically. Please set FEBio_SDK to the path to your FEBio SDK or git repo.")
endif()

# Only update the include and lib directories if the FEBio_SDK path has been changed.
if(NOT OLD_SDK)
    set(NEWPATH TRUE)
else()
    string(COMPARE NOTEQUAL ${FEBio_SDK} ${OLD_SDK} NEWPATH)
endif()

if(NEWPATH)
    message("NEWPATH")
    # Is this the SDK?
    string(REGEX MATCH "sdk" IS_SDK ${FEBio_SDK})

    set(LIB_SUFFIXES "")
    if(IS_SDK)
        set(FEBio_INC "${FEBio_SDK}/include" CACHE PATH "Path to FEBio include directory." FORCE)

        if(WIN32)
            list(APPEND LIB_SUFFIXES "vs2017/Release" "vs2017/Debug")
        else()
            list(APPEND LIB_SUFFIXES "lib")
        endif()
    else()
        set(FEBio_INC ${FEBio_SDK}/../include CACHE PATH "Path to FEBio include directory." FORCE)

        if(WIN32)
            list(APPEND LIB_SUFFIXES "cmbuild/lib/Release" "cmbuild/lib/Debug" "cbuild/lib/Release" "cbuild/lib/Debug" "build/lib/Release" "build/lib/Debug")
        else()
            list(APPEND LIB_SUFFIXES "cbuild/lib" "cmbuild/lib" "build/lib" "cbuild/Release/lib" "cmbuild/Release/lib" "build/Release/lib" "cbuild/Debug/lib" "cmbuild/Debug/lib" "build/Debug/lib")
        endif()
    endif()

    mark_as_advanced(FEBio_INC)

    # Find lib path
    find_library(FECORE
        NAMES FECore fecore fecore_gcc64 fecore_lnx64
        PATHS ${FEBio_SDK}
        PATH_SUFFIXES ${LIB_SUFFIXES}
        DOC "FEBio library path")

    if(FECORE)
        get_filename_component(FECORE_TEMP ${FECORE} DIRECTORY)
        set(FEBio_LIB_DIR ${FECORE_TEMP} CACHE PATH "Path to the FEBio lib directory." FORCE)
        mark_as_advanced(FEBio_LIB_DIR)
        unset(FECORE_TEMP)
        unset(FECORE CACHE)
    else()
        set(FEBio_LIB_DIR CACHE PATH "Path to the FEBio lib directory." FORCE)
        message(SEND_ERROR "Unable to find FEBio Library path automatically. Set FEBio_LIB_DIR.")
        unset(FECORE CACHE)
    endif()
endif()

set(OLD_SDK ${FEBio_SDK} CACHE PATH "Don't edit. Old SDK path used to automatically make changes." FORCE)
mark_as_advanced(OLD_SDK)

#### Check if paths are valid and find libs ####
function(findLib libName libDir libOut)
    string(TOLOWER ${libName} lname)

    find_library(TEMP
        NAMES ${libName} ${lname} ${lname}_gcc64 ${lname}_lnx64
        PATHS ${${libDir}}
        NO_DEFAULT_PATH)

    if(TEMP)
        set(${libOut} ${TEMP} PARENT_SCOPE)
        unset(TEMP CACHE)
    else()
        if(WIN32)
            message(SEND_ERROR "Could not find ${libName}.lib. Check ${libName}.")
        elseif(APPLE)
            message(SEND_ERROR "Could not find ${libName}.so, ${libName}.a, or ${libName}.dylib Check ${libDir}.")
        else()
            message(SEND_ERROR "Could not find ${libName}.so, or ${libName}.a. Check ${libDir}")
        endif()
        unset(TEMP CACHE)
    endif()
endfunction()

if(FEBio_LIB_DIR)
    findLib(FECore FEBio_LIB_DIR FECORE)
    findLib(FEBioMech FEBio_LIB_DIR FEBIOMECH)
    findLib(FEBioMix  FEBio_LIB_DIR FEBIOMIX)
    findLib(NumCore FEBio_LIB_DIR NUMCORE)
    findLib(FEBioFluid  FEBio_LIB_DIR FEBIOFLUID)
endif()

find_package(precice)

##### Set appropriate defines and includes #####

if(WIN32)
    add_definitions(-DWIN32 /MP /openmp)
elseif(APPLE)
    add_definitions(-D__APPLE__)
    set(CMAKE_OSX_DEPLOYMENT_TARGET 10.13)
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -undefined dynamic_lookup")
else()
	add_definitions(-DLINUX)
    add_compile_options(-static-libstdc++ -static-libgcc -w -Wall)

    set(CMAKE_BUILD_RPATH_USE_LINK_PATH FALSE)
    set(CMAKE_BUILD_RPATH $ORIGIN/../lib/)
endif()

# Extra compiler flags for intel compiler
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
	set(CMAKE_SHARED_LINKER_FLAGS ${CMAKE_SHARED_LINKER_FLAGS} -static-intel)
endif()

include_directories(${FEBio_INC})

##### Add library #####

file(GLOB HDR_BFPPlugin "*.h")
file(GLOB SRC_BFPPlugin "*.cpp")

add_library(BFPPlugin SHARED ${HDR_BFPPlugin} ${SRC_BFPPlugin})
set_property(TARGET BFPPlugin PROPERTY AUTOGEN_BUILD_DIR ${CMAKE_BINARY_DIR}/CMakeFiles/AutoGen/BFPPlugin_autogen)

##### Link Libraries #####

target_link_libraries(BFPPlugin ${FECORE} ${FEBIOMIX} ${FEBIOFLUID}  ${FEBIOMECH} ${NUMCORE} ${FECORE} precice::precice ${FEBioTPM})
