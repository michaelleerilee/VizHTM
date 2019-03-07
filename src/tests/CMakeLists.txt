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

set( CMAKE_CXX_FLAGS " ${STARE_CPPFLAGS} ${CMAKE_CXX_FLAGS} ${GEOMPACK_INCLUDE_DIR_FLAG}"  )
set( CMAKE_C_FLAGS   " ${STARE_CFLAGS}   ${CMAKE_C_FLAGS}   ${GEOMPACK_INCLUDE_DIR_FLAG}"  )
message( "CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}" )
message( "CMAKE_C_FLAGS   = ${CMAKE_C_FLAGS}"   )

set (
  TESTS_SRC
  ${TESTS_SRC}
  BoundingBox.C
  Edges.C
  PolePosition.C
)
#  
# TODO Trouble with VarStr.hpp# TODO Trouble with VarStr.hpp
# Var.Str.C is a kluge because I don't know how to get CMake to handle the hpp files correctly.
# INTERFACE May be what works: http://stackoverflow.com/questions/5957134/how-to-setup-cmake-to-generate-header-only-projects

message( "check: CMAKE_SOURCE_DIR = ${CMAKE_SOURCE_DIR}/src" ) 

link_directories( ${GEOMPACK_LIBRARY_DIR} ${STARE_LINK_DIR} )

add_library( testSTARE_tests SHARED ${TESTS_SRC} )
target_include_directories( testSTARE_tests PUBLIC ${CMAKE_SOURCE_DIR}/src /usr/include /usr/include/qt5 /usr/include/qt5/QtCore /usr/include/qt5/QtGui /usr/local/Coin3D/include )



# add_executable( testSTARE testSTARE.C VizHTM_main_support.C ${VIZ_SUPPORT_SRC} )
# add_executable( testSTARE testSTARE.C VizHTM_main_support.C ) 
# target_link_libraries(testSTARE SoQt::SoQt ${STARE_LINK_LIB} VizHTM_LIB ${GEOMPACK_LIBRARY_LIB} shp gfortran )







