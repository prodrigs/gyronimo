cmake_minimum_required(VERSION 3.16)  # CMake version check

# Create project gyronimo
project(GYRONIMO VERSION 0.1 
                 DESCRIPTION "gyromotion for the people, by the people - An object-oriented library for gyromotion applications in plasma physics."
				 HOMEPAGE_URL https://github.com/prodrigs/gyronimo
				 LANGUAGES CXX
)

include(GNUInstallDirs)
include(ExternalProject)

# build options
option(BUILD_VMEC_SUPPORT "Add support to VMEC" OFF)
option(BUILD_NETCDF_CXX4 "Build netcdf-cxx4 from repository" OFF)
option(BUILD_WITH_CLANG "Build with Clang++-11" OFF)
option(BUILD_WITH_DPCPP "Build with Intel's dpcpp" OFF)
option(BUILD_WITH_OPENMP "Build with OpenMP" OFF)


# detect platform
string(COMPARE EQUAL "Darwin"  ${CMAKE_SYSTEM_NAME} OS_X)
string(COMPARE EQUAL "Linux"   ${CMAKE_SYSTEM_NAME} LINUX)
string(COMPARE EQUAL "Windows" ${CMAKE_SYSTEM_NAME} WINDOWS)
if(LINUX)
    set(OS_STRING "linux")
elseif(WINDOWS)
    set(OS_STRING "windows")
elseif(OS_X)
    set(OS_STRING "macOS")
else()
    set(OS_STRING "Unknown")
endif()

# clone external submodules
if (BUILD_VMEC_SUPPORT)
	find_package(Git QUIET)
	if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
	# Update submodules as needed
		option(GIT_SUBMODULE "Check submodules during build" ON)
		if(GIT_SUBMODULE)
			message(STATUS "Submodule update")
			execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
							WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
							RESULT_VARIABLE GIT_SUBMOD_RESULT)
			if(NOT GIT_SUBMOD_RESULT EQUAL "0")
				message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
			endif()
		endif()
	endif()
	# verify if splinter & netcdf-cxx4 were correctly cloned
	if(NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/splinter/CMakeLists.txt" OR 
	NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/netcdf-cxx4/CMakeLists.txt" OR
	NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/SplineLib/CMakeLists.txt")
		message(FATAL_ERROR "The submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
	endif()
endif()

# we need gcc 10 with c++20
if(BUILD_WITH_CLANG)
	# doesn't compile yet (eg. issues with ranges)
        # check: https://wg21.link/p0896r4
	set(CMAKE_C_COMPILER "clang-11")        	                               
	set(CMAKE_CXX_COMPILER "clang++-11")
	#set(CMAKE_CXX_STANDARD 20)            # require c++20 standard
	#set(CMAKE_CXX_STANDARD_REQUIRED True)
	add_definitions(-std=c++20) 
elseif(BUILD_WITH_DPCPP)
	# doesn't compile yet (eg. issues with ranges)
        # check: https://wg21.link/p0896r4
	set(CMAKE_C_COMPILER "dpcpp")        	                               
	set(CMAKE_CXX_COMPILER "dpcpp")
	set(CMAKE_CXX_STANDARD 20)            # require c++20 standard
	set(CMAKE_CXX_STANDARD_REQUIRED True)
else()
	set(CMAKE_C_COMPILER "gcc-10")
	set(CMAKE_CXX_COMPILER "g++-10")
	set(CMAKE_CXX_STANDARD 20)            # require c++20 standard
	set(CMAKE_CXX_STANDARD_REQUIRED True)
endif()
add_definitions(-Wfatal-errors)

if(BUILD_WITH_OPENMP)
	find_package(OpenMP)
	if (OPENMP_FOUND)
		set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
		set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	endif()
endif()

# add netcdf-cxx4 library
if (BUILD_VMEC_SUPPORT)
	ExternalProject_Add(SPLINTER
	SOURCE_DIR ${PROJECT_SOURCE_DIR}/extern/splinter
	PREFIX ${CMAKE_CURRENT_BINARY_DIR}/splinter
	CMAKE_ARGS 
				#  -DCMAKE_CXX_COMPILER=${CMAKE_C_COMPILER}		# splinter doens't compile with g++10
				#  -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}   # TODO: find out why
				-DCMAKE_CXX_FLAGS=-Wno-deprecated-declarations
				-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
	)
	if (BUILD_NETCDF_CXX4)
		# add netcdf-cxx4 library as an external project
		ExternalProject_Add(NETCDF_CXX4
		SOURCE_DIR ${PROJECT_SOURCE_DIR}/extern/netcdf-cxx4
		PREFIX ${CMAKE_CURRENT_BINARY_DIR}/netcdf-cxx4
		CMAKE_ARGS 
					-DCMAKE_CXX_COMPILER=${CMAKE_C_COMPILER}		# splinter doens't compile with g++10
					-DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}   # TODO: find out why
					-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
		)
		if(${OS_STRING} STREQUAL "macOS") 
			set(GYRO_LIBS ${GYRO_LIBS} -L${CMAKE_CURRENT_BINARY_DIR}/netcdf-cxx4/lib -lnetcdf-cxx4)
		elseif(${OS_STRING} STREQUAL "linux")
			set(GYRO_LIBS ${CMAKE_CURRENT_BINARY_DIR}/netcdf-cxx4/lib -lnetcdf-cxx4)
		else()
			set(GYRO_LIBS ${CMAKE_CURRENT_BINARY_DIR}/netcdf-cxx4/lib/libnetcdf-cxx4.a)
		endif()
	else()
		#  - netcdf-cxx4: for supporting VMEC equilibria)
			include(FindPkgConfig)
			pkg_check_modules (NETCDF REQUIRED netcdf-cxx4)
			if(NOT NETCDF_FOUND)
				message(FATAL_ERROR "install netcdf-cxx4 libraries to build VMEC support.")
			endif()
			include_directories(${NETCDF_INCLUDE_DIRS})
			set(GYRO_LIBS ${APPS_LIBS} ${NETCDF_LIBRARIES})
			set(APPS_LIBS ${APPS_LIBS} ${NETCDF_LIBRARIES})
	endif()
endif()

file(GLOB gyronimo_sources ${PROJECT_SOURCE_DIR}/src/*.cc)
if(NOT BUILD_VMEC_SUPPORT)
	list(REMOVE_ITEM gyronimo_sources ${PROJECT_SOURCE_DIR}/src/parser_vmec.cc)
	list(REMOVE_ITEM gyronimo_sources ${PROJECT_SOURCE_DIR}/src/metric_vmec.cc)
endif()

# required libraries

#  - GSL: used by gyronimo
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
set(GYRO_LIBS ${GYRO_LIBS} ${GSL_LIBRARIES})

# - Boost: used by gyronimo and Apps
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost REQUIRED COMPONENTS headers)
#find_package(Boost COMPONENTS *boost libraries here*)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    set(GYRO_LIBS ${GYRO_LIBS} ${Boost_LIBRARIES})
    set(APPS_LIBS ${APPS_LIBS} ${Boost_LIBRARIES})
endif()

# build gyronimo library
# tbi: add_subdirectory("${PROJECT_SOURCE_DIR}/src")
add_library(gyronimo SHARED ${gyronimo_sources})
#file(GLOB gyronimo_headers ${PROJECT_SOURCE_DIR}/include/gyronimo/*/*.hh)
set_target_properties(gyronimo PROPERTIES 
	VERSION ${PROJECT_VERSION} 
	SOVERSION 1)
#	PUBLIC_HEADER ${PROJECT_SOURCE_DIR}/include)


if(BUILD_VMEC_SUPPORT)
	# add gyronimo external dependencies 
	# - splinter library
	add_dependencies(gyronimo SPLINTER)
	# - netcdf-cxx4 library for VMEC parsing
	if(BUILD_NETCDF_CXX4)
		add_dependencies(gyronimo NETCDF_CXX4)
	# else()
	# 	add_dependencies(gyronimo netcdf)
	endif()
	if(${OS_STRING} STREQUAL "macOS") 
		set(SPLINTER_LIB -L${CMAKE_CURRENT_BINARY_DIR}/splinter/lib -lsplinter-3-0)
	elseif(${OS_STRING} STREQUAL "linux")
		set(SPLINTER_LIB -L${CMAKE_CURRENT_BINARY_DIR}/splinter/lib -lsplinter-3-0)
	else()
		set(SPLINTER_LIB ${CMAKE_CURRENT_BINARY_DIR}/splinter/lib/libsplinter-static-3-0.a)
	endif()
    set(GYRO_LIBS ${GYRO_LIBS} ${SPLINTER_LIB})
	set(GYRO_INCLUDE ${GYRO_INCLUDE} ${CMAKE_CURRENT_BINARY_DIR}/splinter/include/SPLINTER)
endif()

configure_file(${PROJECT_SOURCE_DIR}/include/gyronimo/core/version.hh.in version.hh)
target_include_directories(gyronimo PUBLIC 
	                   ${CMAKE_CURRENT_BINARY_DIR}
	                   ${PROJECT_SOURCE_DIR}/include
					   ${GYRO_INCLUDE}
 )


target_link_libraries(gyronimo PUBLIC ${GYRO_LIBS}
						
)
list(APPEND INCLUDES "${PROJECT_BINARY_DIR}")
list(APPEND INCLUDES "${PROJECT_SOURCE_DIR}/include")
list(APPEND INCLUDES "${PROJECT_SOURCE_DIR}/apps/include")

# add the executables
file (GLOB targets ${PROJECT_SOURCE_DIR}/apps/*.cc)
if(NOT BUILD_VMEC_SUPPORT)
	list(REMOVE_ITEM targets ${PROJECT_SOURCE_DIR}/apps/vmecdump.cc)
endif()

foreach (target_file IN LISTS targets)
	get_filename_component (target ${target_file} NAME_WE)
	add_executable (${target}.x ${target_file})
	target_include_directories(${target}.x PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           ${INCLUDES}
                           )
	target_link_libraries (${target}.x PRIVATE gyronimo ${APPS_LIBS})
	install (TARGETS ${target}.x DESTINATION bin)
endforeach()

# build documentation
find_package(Doxygen)
option(BUILD_DOCUMENTATION "builds API documentation (requires Doxygen)" ${DOXYGEN_FOUND})
if(BUILD_DOCUMENTATION)
    if(NOT DOXYGEN_FOUND)
         message(FATAL_ERROR "Doxygen is needed to build the documentation.")
    endif()
    set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/doc/config.dox)
    set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/doc/config.dox)
    configure_file(${doxyfile_in} ${doxyfile} @ONLY)
    message("Doxygen build started.")
    add_custom_target(doc #ALL
                      COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile_in}
		      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
                      COMMENT "Generating API documentation with Doxygen"
                      VERBATIM)
endif()

# install gyronimo library
install (TARGETS gyronimo
                LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
		# PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install (DIRECTORY ${CMAKE_SOURCE_DIR}/include/ 
	 DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
         FILES_MATCHING PATTERN "*.hh")

# install documentation
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/gyronimo DESTINATION share/doc)

# pkg-config support
configure_file(${PROJECT_SOURCE_DIR}/src/gyronimo.pc.in gyronimo.pc @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/gyronimo.pc DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)
