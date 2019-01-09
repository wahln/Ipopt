cmake_minimum_required(VERSION 3.10)

project(IpOpt)


option(IPOPT_HAS_EIGEN			   	   	"Enable Eigen solvers"  ON)
option(IPOPT_USE_OPENBLAS			   	"Use OpenBLAS instead of reference BLAS" ON)
option(IPOPT_HAS_MUMPS					"Enable MUMPS solver"	OFF)
option(IPOPT_USE_CLAPACK				"Use CLAPACK" ON)

#Third Party
include(ExternalProject)
set(THIRDPARTY_DIR ${CMAKE_BINARY_DIR}/ThirdParty)
set(CMAKE_MSVCIDE_RUN_PATH  "${THIRDPARTY_DIR}/bin")

#-----------------------------------------------------------------------------
# OpenBLAS
#-----------------------------------------------------------------------------

set(IPOPT_DEPENDS )

# set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libs with openBLAS")



# Sanity checks
if(IPOPT_USE_OPENBLAS)
	if(DEFINED OpenBLAS_DIR AND NOT EXISTS ${OpenBLAS_DIR})
		message(FATAL_ERROR "OpenBLAS_DIR variable is defined but corresponds to non-existing directory")
	endif()
	
#	find_package(Git)
#	if(NOT GIT_FOUND)
#		message(FATAL_ERROR "OpenBLAS requires Git")
#	endif()

  if(NOT DEFINED OpenBLAS_DIR)

    ExternalProject_Add(OpenBLAS
		LIST_SEPARATOR ${sep}
      	URL https://github.com/xianyi/OpenBLAS/archive/v0.3.5.tar.gz
 	   	URL_MD5 579BDA57F68EA6E9074BF5780E8620BB
#		GIT_REPOSITORY https://github.com/wahln/OpenBLAS.git
#		GIT_TAG 3ff64abbc8faef52006f5f19d5615c5e677fffa2
#		GIT_SHALLOW	ON
		PREFIX ${THIRDPARTY_DIR}/
		CMAKE_CACHE_ARGS 
			-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
#			-DBUILD_SHARED_LIBS:BOOL=OFF
			-DCMAKE_RUNTIME_OUTPUT_DIRECTORY:PATH=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
	  )

    # set(Eigen_DIR ${ep_prefix})
    # mitkFunctionInstallExternalCMakeProject(${proj})
	list(APPEND IPOPT_DEPENDS OpenBLAS)
	
  else()

    # mitkMacroEmptyExternalProject(${proj} "${proj_DEPENDENCIES}")

  endif()
	ExternalProject_Get_Property(OpenBLAS INSTALL_DIR)
  message(STATUS "OpenBLAS Install dir: ${INSTALL_DIR}")
  
  	set(IPOPT_BLAS_LIBRARY ${THIRDPARTY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}openblas${CMAKE_STATIC_LIBRARY_SUFFIX} CACHE FILEPATH "Path to static blas library")
	mark_as_advanced(IPOPT_BLAS_LIBRARY)
endif()

if(IPOPT_USE_CLAPACK)
	if(DEFINED CLAPACK_DIR AND NOT EXISTS ${CLAPACK_DIR})
		message(FATAL_ERROR "CLAPACK_DIR variable is defined but corresponds to non-existing directory")
	endif()
	
	if(NOT DEFINED CLAPACK_DIR)
		ExternalProject_Add(CLAPACK
			LIST_SEPARATOR ${sep}
			URL http://icl.cs.utk.edu/lapack-for-windows/clapack/clapack-3.2.1-CMAKE.tgz
			URL_MD5 EE4947C1FE86D4114F2559003711027A
			PREFIX ${THIRDPARTY_DIR}/
			CMAKE_CACHE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
			)
	else()
	
	endif()
	
	ExternalProject_Get_Property(CLAPACK INSTALL_DIR)
	message(STATUS "CLAPACK Install dir: ${INSTALL_DIR}")
	
	set(IPOPT_LAPACK_LIBRARY ${THIRDPARTY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}lapack${CMAKE_STATIC_LIBRARY_SUFFIX} CACHE FILEPATH "Path to static lapack library")
	mark_as_advanced(IPOPT_LAPACK_LIBRARY)
	
	list(APPEND IPOPT_DEPENDS CLAPACK)
	
endif()
#-----------------------------------------------------------------------------
# Eigen
#-----------------------------------------------------------------------------
if(IPOPT_HAS_EIGEN)

  # Sanity checks
  if(DEFINED Eigen_DIR AND NOT EXISTS ${Eigen_DIR})
    message(FATAL_ERROR "Eigen_DIR variable is defined but corresponds to non-existing directory")
  endif()

  if(NOT DEFINED Eigen_DIR)

    ExternalProject_Add(Eigen
      LIST_SEPARATOR ${sep}
      URL http://bitbucket.org/eigen/eigen/get/3.3.7.tar.bz2
	  URL_MD5 05B1F7511C93980C385EBE11BD3C93FA
	  DOWNLOAD_NAME eigen3.tar.bz2
	  PREFIX ${THIRDPARTY_DIR}/
#	  SOURCE_DIR <PREFIX>
      CMAKE_ARGS
 #       ${ep_common_args}
        -DBUILD_TESTING:BOOL=ON
        -DEIGEN_BUILD_PKGCONFIG:BOOL=OFF
      CMAKE_CACHE_ARGS
	  	-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
		-DBLAS_openblas_LIBRARY:FILEPATH=${IPOPT_BLAS_LIBRARY}
 #       ${ep_common_cache_args}
#      CMAKE_CACHE_DEFAULT_ARGS
 #       ${ep_common_cache_default_args}
	  DEPENDS OpenBLAS
    )
	ExternalProject_Get_Property(Eigen SOURCE_DIR)
	# message("Source dir of Eigen: ${SOURCE_DIR}")
    set(Eigen_DIR ${SOURCE_DIR} CACHE PATH "Path to Eigen template library")

	# ExternalProject_Get_Property(Eigen INSTALL_DIR)
	# message(STATUS "Eigen Install dir: ${INSTALL_DIR}")
	

  else()
	ExternalProject_Add(Eigen
	  SOURCE_DIR ${Eigen_DIR}	
	  INSTALL_DIR ${THIRDPARTY_DIR}/
	  BINARY_DIR ${THIRDPARTY_DIR}/src/Eigen-build
      CMAKE_ARGS
 #       ${ep_common_args}
        -DBUILD_TESTING:BOOL=ON
        -DEIGEN_BUILD_PKGCONFIG:BOOL=OFF
      CMAKE_CACHE_ARGS
	  	-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
		-DBLAS_openblas_LIBRARY:FILEPATH=${IPOPT_BLAS_LIBRARY}
 #       ${ep_common_cache_args}
#      CMAKE_CACHE_DEFAULT_ARGS
 #       ${ep_common_cache_default_args}
	  DEPENDS OpenBLAS
    )

  endif()
  mark_as_advanced(Eigen_DIR)
  list (APPEND IPOPT_DEPENDS Eigen)
endif()


#-----------------------------------------------------------------------------
# OpenMP
#-----------------------------------------------------------------------------

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

#-----------------------------------------------------------------------------
# AMPL / ASL
#-----------------------------------------------------------------------------

# if(IPOPT_USE_ASL)



# ExternalProject_Add(ASL
      # LIST_SEPARATOR ${sep}
      # URL https://github.com/ampl/mp/archive/1.3.0.tar.gz
	  # URL_MD5 05B1F7511C93980C385EBE11BD3C93FA
	  # PREFIX ${THIRDPARTY_DIR}/

	  # DEPENDS OpenBLAS
    # )

# endif()


subdirs(Ipopt)

