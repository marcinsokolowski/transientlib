cmake_minimum_required(VERSION 3.10.0)
project(transientlib VERSION 1.0.1)
include(GNUInstallDirs)

file(GLOB transientlib_sources "src/*.cpp")
file(GLOB transientlib_sources_c "src/*.c")
file(GLOB transientlib_headers "src/*.h")
add_library(transientlib SHARED ${transientlib_sources} ${transientlib_sources_c})

add_definitions(-D_UNIX -D_ENABLE_PROFILER_ -D_ENABLE_TRACING_ -D_ENABLE_VISUALIZATION_ -D_MONTE_CARLO_ -D_NO_ROOT_ -I${CMAKE_SOURCE_DIR}/src)

set_target_properties(transientlib PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
    PUBLIC_HEADER "${transientlib_headers}")
 
configure_file(transientlib.pc.in transientlib.pc @ONLY)


install(TARGETS transientlib
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(FILES ${CMAKE_BINARY_DIR}/transientlib.pc
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)


find_library(CFITSIO_LIB cfitsio HINTS ENV LD_LIBRARY_PATH)
# find_library(FFTW3_LIB fftw3 HINTS ENV FFTW_LIB REQUIRED)
# find_library(LIBNOVA_LIB nova PATHS ENV LD_LIBRARY_PATH)


# required to properly link the dynamic library :
target_link_libraries(transientlib ${CFITSIO_LIB} ${LIBNOVA_LIB} ${ROOT_LIBRARIES} ${FFTW3_LIB} -ldl -lpthread)

# add_executable(nan_test apps/nan_test.cpp)
# add_executable(libtest  apps/libtest.cpp)
# target_link_libraries(libtest transientlib ${CFITSIO_LIB} ${LIBNOVA_LIB} ${ROOT_LIBRARIES} ${FFTW3_LIB} -ldl -lpthread)
# add_executable(radec2azh apps/radec2azh.cpp)
# target_link_libraries(radec2azh transientlib ${CFITSIO_LIB} ${LIBNOVA_LIB} ${ROOT_LIBRARIES} ${FFTW3_LIB} -ldl -lpthread)
# add_executable(sid2ux apps/sid2ux.cpp)
# target_link_libraries(sid2ux transientlib ${CFITSIO_LIB} ${LIBNOVA_LIB} ${ROOT_LIBRARIES} ${FFTW3_LIB} -ldl -lpthread)
# add_executable(ux2sid   apps/ux2sid.cpp) 
# target_link_libraries(ux2sid transientlib ${CFITSIO_LIB} ${LIBNOVA_LIB} ${ROOT_LIBRARIES} ${FFTW3_LIB} -ldl -lpthread)
# add_executable(ux2sid_file   apps/ux2sid_file.cpp)
# target_link_libraries(ux2sid_file transientlib ${CFITSIO_LIB} ${LIBNOVA_LIB} ${ROOT_LIBRARIES} ${FFTW3_LIB} -ldl -lpthread)
# add_executable(eda_spectrometer   apps/main_fft_file.cpp)
# target_link_libraries(eda_spectrometer transientlib ${CFITSIO_LIB} ${LIBNOVA_LIB} ${ROOT_LIBRARIES} ${FFTW3_LIB} -ldl -lpthread)

# larger programs :
# add_executable(avg_images  apps/avg_images.cpp)
# target_link_libraries(avg_images transientlib ${CFITSIO_LIB} ${LIBNOVA_LIB} ${ROOT_LIBRARIES} ${FFTW3_LIB} -ldl -lpthread)

# add_executable(calcfits_bg  apps/calcfits_bg.cpp)
# target_link_libraries(calcfits_bg transientlib ${CFITSIO_LIB} ${LIBNOVA_LIB} ${ROOT_LIBRARIES} ${FFTW3_LIB} -ldl -lpthread)

add_executable(sattest  apps/sattest/main.cpp)
target_link_libraries(sattest transientlib ${CFITSIO_LIB} ${LIBNOVA_LIB} ${ROOT_LIBRARIES} ${FFTW3_LIB} -ldl -lpthread -lcurses)

# INSTALLATION:
install(TARGETS sattest RUNTIME DESTINATION bin)
