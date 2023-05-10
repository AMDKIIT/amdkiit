#**************************************************************************
# INSTRUCTIONS
#**************************************************************************
# 1. Always use CMake Version 3.4 or higher
# 2. Check path where library(MPI,LAPACK,BLAS,LIBXC) is installed, give that path within " " 
# EXAMPLE: LIBexample is installed in /home/usrname/LIBexample
# then change, 
# set(LIBexample_ROOT "/home/usrname/LIBexample")
#**************************************************************************

cmake_minimum_required(VERSION 3.4)

project(amdkiit VERSION 1.0.0 LANGUAGES CXX C Fortran) 
message(STATUS "Project version ${CMAKE_PROJECT_VERSION}")

string(TIMESTAMP CURRENT_TIME "%d-%m-%Y %H:%M:%S")
message("Current time is ${CURRENT_TIME}")

configure_file(version.h.in ../version.f90)


#by default GNU compiler will be used, to use INTEL uncomment the following
add_compile_definitions(PROJECT_VERSION="${CMAKE_PROJECT_VERSION}")
set(CMAKE_Fortran_COMPILER_ID "GNU")

# Set "/path/to/directory" for libraries
set(MATH_ROOT "/opt/nnair.iitk/libs") # CHECK PATH
set(FFTW_ROOT "/home/nnair.iitk/libs") # CHECK PATH
set(LIBXC_ROOT "/home/nnair.iitk/libs/") # CHECK PATH

#**************************************************************************
# Intel Fortran
#**************************************************************************
if( CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" )
     message(STATUS "FORTRAN COMPILER: " ${CMAKE_Fortran_COMPILER_ID})

     #set "/path/to/directory" for  MKL and other libraries
     set(MATH_ROOT $ENV{MKLROOT})
     set(LIBXC_ROOT "/home/paramita/libxc/") #"CHECK PATH"
     
     find_package(MPI REQUIRED)
     include_directories(${MPI_Fortran_INCLUDE_PATH})
     
     include_directories(${MATH_ROOT}/include)
     include_directories(${LIBXC_ROOT}/include)
  
   set(MATH_LIB "-L/${MATH_ROOT}/lib/intel64/ -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -ldl -lm -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64")
   set(LIBXC_LIB "-L${LIBXC_ROOT}lib -lxcf90 -lxc")
   
   set (LIBS "${MPI_Fortran_LIBRARIES}" "${MATH_LIB}" "${LIBXC_LIB}")

#**************************************************************************
# GNU 
#**************************************************************************
elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
     message(STATUS "FORTRAN COMPILER: " ${CMAKE_Fortran_COMPILER_ID})

    
     find_package(MPI REQUIRED)
     include_directories(${MPI_Fortran_INCLUDE_PATH})
     
     include_directories(${FFTW_ROOT}include)
     include_directories(${LIBXC_ROOT}include)
     
     set(MATH_LIB "-L${MATH_ROOT} -llapack -lrefblas")
     set(FFTW_LIB "-L${FFTW_ROOT} -lfftw3")
     set(LIBXC_LIB "-L${LIBXC_ROOT}lib -lxcf90 -lxc")
    
     set (LIBS "${MPI_Fortran_LIBRARIES}" "${MATH_LIB}" "${FFTW_LIB}" "${LIBXC_LIB}")


else()
   message(FATAL_ERROR "could not find any compiler")
endif()
#*******************************************************************************************
#NOTHING TO CHANGE 
#*******************************************************************************************
add_subdirectory(yaml-cpp)
set(CMAKE_CXX_STANDARD 11)
add_library(fortran_yaml_cpp yaml_types.f90 yaml.f90 yaml.cpp)
target_link_libraries(fortran_yaml_cpp yaml-cpp)

add_executable(${CMAKE_PROJECT_NAME}.x print_logo_amdkiit.f90 yamlread.f90 system_data_types.f90 max_parameter_pp.f90 kinds.f90 atom.f90 atomic_basis.f90 cg.f90 constants.f90 density.f90 exp_igr.f90 fft_interface.f90 fft_support.f90 gradient.f90 gvectors.f90 interpolgrid.f90 kinds.f90 lsfbtr.f90 main.f90 math.f90 max_parameter_pp.f90 mympi.f90 nuclear_grad.f90 pbc.f90 potential.f90 pseudopotential.f90 read_pp_file.f90 readpsp.f90 readstring.f90 set_cell.f90 system_data_types.f90 total_energy.f90 wfn_initialize.f90 xc.f90 bfgs_mod.f90 function_val.f90 backtracking.f90 readupf.f90 spline_inter.f90 spline.f90 spherical_func.f90 bessel_func.f90 init_module.f90 do_opt.f90 version.f90 eprint.f90 set_environment.f90)

target_link_libraries(${CMAKE_PROJECT_NAME}.x fortran_yaml_cpp ${LIBS})

# Define the version number macro
target_compile_definitions(${CMAKE_PROJECT_NAME}.x PRIVATE VERSION="${CMAKE_PROJECT_VERSION}")

# Set the Fortran standard to 2008
set_target_properties(${CMAKE_PROJECT_NAME}.x PROPERTIES
    Fortran_STANDARD 2008
    Fortran_STANDARD_REQUIRED ON
)

