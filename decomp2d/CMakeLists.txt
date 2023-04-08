file(GLOB files_decomp decomp_2d.f90
                       glassman.f90
		       io.f90)
include_directories(${CMAKE_SOURCE_DIR}/decomp2d)
if(${FFT_Choice} MATCHES "generic")
  file(GLOB files_fft fft_generic.f90)
elseif(${FFT_Choice} MATCHES "fftw3")
  file(GLOB files_fft fft_fftw3.f90)
elseif(${FFT_Choice} MATCHES "mkl")
  file(GLOB files_fft fft_mkl.f90)
endif(${FFT_Choice} MATCHES "generic")

set(SRCFILES ${files_decomp} ${files_fft})

add_library(decomp2d STATIC ${SRCFILES})
if (MPI_FOUND)
  target_link_libraries(decomp2d PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND) 
if(ADIOS2_FOUND)
  target_link_libraries(decomp2d PRIVATE adios2::fortran_mpi)
endif()
if(${FFT_Choice} MATCHES "fftw3")
  # target_link_libraries(decomp2d PRIVATE fftw3)
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(FFTW REQUIRED IMPORTED_TARGET fftw3)

  target_link_libraries(decomp2d PUBLIC ${FFTW_LIBRARIES})
  target_include_directories(decomp2d PUBLIC ${FFTW_INCLUDE_DIRS})

endif()


install(TARGETS decomp2d
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} 
)
