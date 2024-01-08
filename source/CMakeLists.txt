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
set(USE_INTEL_MPI ON)
set(USE_MKL ON)

if(USE_INTEL_MPI)

    set(CMAKE_C_COMPILER "mpiicc")
    set(CMAKE_CXX_COMPILER "mpiicpc")
    set(CMAKE_Fortran_COMPILER "mpiifort")

    if(USE_MKL)
      set(LIBXC_ROOT "/home/paramitag/amkdiit_intel/libxc/")  # change path for INTEL MKL COMPILER
    else()
      set(MATH_ROOT "/home/paramitag/amkdiit_intel/lapack-3.11/")  # change path for INTEL COMPILER
      set(FFTW_ROOT "/home/paramitag/amkdiit_intel/fftw/") 
      set(LIBXC_ROOT "/home/paramitag/amkdiit_intel/libxc/") 
   endif()

else()

    set(CMAKE_C_COMPILER "mpicc")
    set(CMAKE_CXX_COMPILER "mpicxx")
    set(CMAKE_Fortran_COMPILER "mpif90")
                                                              
     set(MATH_ROOT "/home/paramitag/amdkiit_gnu/lapack-3.11/") # change path for GNU COMPILER 
     set(FFTW_ROOT "/home/paramitag/amdkiit_gnu/fftw/") 
     set(LIBXC_ROOT "/home/paramitag/amdkiit_gnu/libxc/") 
endif()


project(amdkiit VERSION 1.0.1 LANGUAGES CXX C Fortran) 
find_package(MPI REQUIRED)
message(STATUS " ")
message(STATUS " *********************** ")

message(STATUS " PROJECT VERSION: ${CMAKE_PROJECT_VERSION} ")

string(TIMESTAMP CURRENT_DATE "%d-%m-%Y") # %H:%M:%S")
message(STATUS " DATE: ${CURRENT_DATE} ")

string(TIMESTAMP CURRENT_TIME "%H:%M:%S")
message(STATUS " TIME: ${CURRENT_TIME} ")

configure_file(version.h.in ../version.f90)

#**************************************************************************
# Intel Fortran
#**************************************************************************
if( CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" )

   message(STATUS " FORTRAN COMPILER: " ${CMAKE_Fortran_COMPILER_ID})
   message(STATUS " ")

   if(USE_MKL)
     set(MATH_ROOT $ENV{MKLROOT})
     include_directories(${MATH_ROOT}/include)
     set(MATH_LIB "-L${MATH_ROOT}/lib/intel64/ -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64")
   
   else()
     set(MATH_LIB "-L${MATH_ROOT} -llapack -lrefblas")

     include_directories(${FFTW_ROOT}include)
     set(FFTW_LIB "-L${FFTW_ROOT}lib -lfftw3")
   endif()

   include_directories(${LIBXC_ROOT}/include)
   set(LIBXC_LIB "-L${LIBXC_ROOT}lib -lxcf90 -lxc")
   
   set (LIBS "${MPI_Fortran_LIBRARIES} ${MATH_LIB} ${FFTW_LIB} ${LIBXC_LIB}")
   message(STATUS " *********************** ")     
   message(STATUS " MATH LIBS:" ${MATH_LIB})
   message(STATUS " ")
   message(STATUS " FFTW LIBS:" ${FFTW_LIB})
   message(STATUS " ")
   message(STATUS " LIBXC LIBS:" ${LIBXC_LIB})   
   message(STATUS " *********************** ")
   message(STATUS " ")
 
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -pg -O0 -g -traceback -fpe0 -check all -heap-arrays")

#**************************************************************************
# GNU 
#**************************************************************************

elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    
   message(STATUS " FORTRAN COMPILER: " ${CMAKE_Fortran_COMPILER_ID})
   message(STATUS " ")

   include_directories(${FFTW_ROOT}include)
   include_directories(${LIBXC_ROOT}include)

   set(MATH_LIB "-L${MATH_ROOT} -llapack -lrefblas")
   set(FFTW_LIB "-L${FFTW_ROOT}lib -lfftw3")
   set(LIBXC_LIB "-L${LIBXC_ROOT}lib -lxcf90 -lxc")

   set (LIBS "${MPI_Fortran_LIBRARIES}" "${MATH_LIB}" "${FFTW_LIB}" "${LIBXC_LIB}")
   message(STATUS " *********************** ")
   message(STATUS " MATH LIBS:" ${MATH_LIB})
   message(STATUS " ")
   message(STATUS " FFTW LIBS:" ${FFTW_LIB})
   message(STATUS " ")
   message(STATUS " LIBXC LIBS:" ${LIBXC_LIB})
   message(STATUS " *********************** ")
   message(STATUS " ")

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -pg -O0 -g -fbacktrace -ffpe-trap=zero -fcheck=all")
else()
   message(FATAL_ERROR "could not find any compiler")
endif()
#*******************************************************************************************
#NOTHING TO CHANGE 
#*******************************************************************************************
set(CMAKE_BUILD_TYPE Debug)

include(FortranCInterface)
FortranCInterface_VERIFY(CXX)

add_subdirectory(yaml-cpp EXCLUDE_FROM_ALL)

set(CMAKE_CXX_STANDARD 11)
add_library(fortran_yaml_cpp yaml_types.f90 yaml.f90 yaml.cpp)
target_link_libraries(fortran_yaml_cpp yaml-cpp)
target_include_directories(fortran_yaml_cpp PUBLIC yaml-cpp/include)


add_executable(${CMAKE_PROJECT_NAME}.x lbfgs.f print_logo_amdkiit.f90 yamlread.f90 system_data_types.f90 max_parameter_pp.f90 kinds.f90 atom_info.f90 atomic_basis.f90 cg.f90 constants.f90 density.f90 exp_igr.f90 fft_interface.f90 fft_support.f90 gradient.f90 gvectors.f90 interpolgrid.f90 kinds.f90 sbt.f90 main.f90 math.f90 max_parameter_pp.f90 mympi.f90 nuclear_grad.f90 apply_pbc.f90 potential.f90 pseudopotential.f90 read_pp_file.f90 readpsp.f90 readstring.f90 set_cell.f90 system_data_types.f90 total_energy.f90 wfn_initialize.f90 xc.f90 bfgs_mod.f90 function_val.f90 backtracking.f90 readupf.f90 spline_inter.f90 spline.f90 spherical_func.f90 bessel_func.f90 init_module.f90 do_opt.f90 version.f90 eprint.f90 set_environment.f90 random.f90 md_tools.f90 update_velver.f90 molecular_dynamics.f90)
set_target_properties(${CMAKE_PROJECT_NAME}.x PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(${CMAKE_PROJECT_NAME}.x fortran_yaml_cpp)
target_link_libraries(${CMAKE_PROJECT_NAME}.x ${FFTW_LIB})
target_link_libraries(${CMAKE_PROJECT_NAME}.x ${MATH_LIB})
target_link_libraries(${CMAKE_PROJECT_NAME}.x ${LIBXC_LIB})

target_include_directories(${CMAKE_PROJECT_NAME}.x PRIVATE ${FFTW_ROOT}include)
target_include_directories(${CMAKE_PROJECT_NAME}.x PRIVATE ${LIBXC_ROOT}include)
# Define the version number macro
target_compile_definitions(${CMAKE_PROJECT_NAME}.x PRIVATE VERSION="${CMAKE_PROJECT_VERSION}")
