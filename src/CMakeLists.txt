cmake_minimum_required(VERSION 2.6)

project (VizHTM)

message( "test: CMAKE_PREFIX_PATH $ENV{CMAKE_PREFIX_PATH}" )
message( "test: COIN3D_CONFIG_COMMAND $ENV{COIN3D_CONFIG_COMMAND}" )

#    find_package(SoQt REQUIRED)
#    target_link_libraries(my_target_name SoQt::SoQt)

find_package(Coin REQUIRED)
message( "test: COIN_NAME  ${COIN_NAME}" )
message( "test: COIN_VERSION  ${COIN_VERSION}" )

find_package(SoQt REQUIRED)
message( "test: SOQT_VERSION  ${SOQT_VERSION}" )

set( GEOMPACK_INCLUDE_DIR_FLAG "-I/home/mrilee/src/geompack/geompack-ree/" )
set( GEOMPACK_LIBRARY_DIR "/home/mrilee/src/geompack/geompack-ree/" )
set( GEOMPACK_LIBRARY_LIB "geompack" )

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set( STARE_CPPFLAGS  "-std=c++11 -g -Wall -D_BOOL_EXISTS -D__macosx -UDIAGNOSE -fPIC -I/home/mrilee/workspace/STARE-CMAKE/include/" )
  set( STARE_CFLAGS    "-g -Wall -D_BOOL_EXISTS -D__macosx -UDIAGNOSE -fPIC -I/home/mrilee/workspace/STARE-CMAKE/include/" )
  set( STARE_LINK_DIR "/home/mrilee/workspace/STARE-CMAKE/build/default/src/" )  
  set( STARE_LINK_LIB "STARE" )
elseif(UNIX)
  set( STARE_CPPFLAGS  "-std=c++11 -g -Wall -D_BOOL_EXISTS -D__unix -UDIAGNOSE -fPIC -I/home/mrilee/workspace/STARE-CMAKE/include/" )
  set( STARE_CFLAGS    "-g -Wall -D_BOOL_EXISTS -D__unix -UDIAGNOSE -fPIC -I/home/mrilee/workspace/STARE-CMAKE/include/" )
  set( STARE_LINK_DIR "/home/mrilee/workspace/STARE-CMAKE/build/default/src/" )  
  set( STARE_LINK_LIB "STARE" )
else()
  message(FATAL_ERROR "Platform ${CMAKE_SYSTEM_NAME} not handled.")
endif()

# HDF4 for gcc
set( HDF4_DIR                 /usr/local/hdf4 )
set( HDF4_INCLUDE_DIRECTORIES ${HDF4_DIR}/include ) 
set( HDF4_CFLAGS  "-I/usr/include/tirpc -O3 -fomit-frame-pointer" )  
set( HDF4_LDFLAGS "-L/usr/local/lib:/usr/lib64 -L/usr/lib64 -L/usr/local/lib -L/usr/local/hdf4/lib " )
set( HDF4_LIBS    "${HDF4_DIR}/lib/libmfhdf.a ${HDF4_DIR}/lib/libdf.a -lsz -ljpeg -lz -ltirpc " )
set( HDF4_LIBRARIES ${HDF4_DIR}/lib/libmfhdf.a ${HDF4_DIR}/lib/libdf.a /usr/local/lib/libsz.so /usr/lib64/libjpeg.so /usr/lib64/libz.so /usr/lib64/libtirpc.so )

set( HDFEOS2_DIR     /usr/local/hdfeos2 )
set( HDFEOS2_CFLAGS  "-I${HDFEOS2_DIR}/include" )
set( HDFEOS2_LDFLAGS "-L${HDFEOS2_DIR}/lib" )
set( HDFEOS2_LIBS "-lhdfeos -lGctp " )
set( HDFEOS2_LIBRARIES "${HDFEOS2_DIR}/lib/libhdfeos.a" "${HDFEOS2_DIR}/lib/libGctp.a" )
set( HDFEOS2_INCLUDE_DIRECTORIES "${HDFEOS2_DIR}/include" )

set( HDF_INCLUDE_DIRECTORIES "${HDF4_INCLUDE_DIRECTORIES} ${HDFEOS2_INCLUDE_DIRECTORIES}" ) 
set( HDF_LINK_DIRECTORIES ${HDFEOS2_DIR}/lib ${HDF4_DIR}/lib )
set( HDF_LIBRARIES ${HDFEOS2_LIBRARIES} ${HDF4_LIBRARIES} )

set( CMAKE_CXX_FLAGS " ${STARE_CPPFLAGS} ${CMAKE_CXX_FLAGS} ${GEOMPACK_INCLUDE_DIR_FLAG}"  )
set( CMAKE_C_FLAGS   " ${STARE_CFLAGS}   ${CMAKE_C_FLAGS}   ${GEOMPACK_INCLUDE_DIR_FLAG}"  )
message( "CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}" )
message( "CMAKE_C_FLAGS   = ${CMAKE_C_FLAGS}"   )

set (
 VIZ_SUPPORT_SRC
 misc.C
 OffScreenViz.C
 VizHTM.C
 )

# set (
#  STARE_SrcFiles
#  ${STARE_SrcFiles}
#  VarStr.C )
#  
# TODO Trouble with VarStr.hpp# TODO Trouble with VarStr.hpp
# Var.Str.C is a kluge because I don't know how to get CMake to handle the hpp files correctly.
# INTERFACE May be what works: http://stackoverflow.com/questions/5957134/how-to-setup-cmake-to-generate-header-only-projects

link_directories( ${GEOMPACK_LIBRARY_DIR} ${STARE_LINK_DIR} ${HDF_LINK_DIRECTORIES} )

add_library( VizHTM_LIB SHARED ${VIZ_SUPPORT_SRC} )
target_include_directories( VizHTM_LIB PUBLIC /usr/include /usr/include/qt5 /usr/include/qt5/QtCore /usr/include/qt5/QtGui /usr/local/Coin3D/include ${HDF_INCLUDE_DIRECTORIES} )

add_subdirectory( tests )

# add_executable( testSTARE testSTARE.C VizHTM_main_support.C ${VIZ_SUPPORT_SRC} )

add_executable( testSTARE testSTARE.C VizHTM_main_support.C ) 
target_link_libraries(testSTARE SoQt::SoQt testSTARE_tests ${STARE_LINK_LIB} VizHTM_LIB ${GEOMPACK_LIBRARY_LIB} ${HDF_LIBRARIES} shp gfortran erfa )
target_include_directories( testSTARE PUBLIC ${CMAKE_SOURCE_DIR}/src /usr/include /usr/include/qt5 /usr/include/qt5/QtCore /usr/include/qt5/QtGui /usr/local/Coin3D/include ${HDF_INCLUDE_DIRECTORIES} )

add_executable( VizHTM_main VizHTM_main.C ) 
target_link_libraries(VizHTM_main SoQt::SoQt ${STARE_LINK_LIB} VizHTM_LIB ${GEOMPACK_LIBRARY_LIB} ${HDF_LIBRARIES} shp gfortran erfa )
target_include_directories( VizHTM_main PUBLIC ${CMAKE_SOURCE_DIR}/src /usr/include /usr/include/qt5 /usr/include/qt5/QtCore /usr/include/qt5/QtGui /usr/local/Coin3D/include ${HDF_INCLUDE_DIRECTORIES} )

# target_include_directories( VizHTM_LIB PUBLIC /usr/include /usr/include/qt5 /usr/include/qt5/QtCore /usr/include/qt5/QtGui /usr/local/Coin3D/include )




