cmake_minimum_required( VERSION 2.8.4 )
project( imgtabdet )
set( tool_EXE imgtabdet )
include( CheckIncludeFileCXX )
include( FindPkgConfig )
pkg_check_modules( MagickCore REQUIRED MagickCore>=6.8 )
pkg_check_modules( opencv REQUIRED opencv>=2.4 )

set( CMAKE_REQUIRED_INCLUDES "${CMAKE_REQUIRED_INCLUDES};${opencv_INCLUDE_DIRS}" )

check_include_file_cxx( opencv2/line_descriptor.hpp ENABLE_OPENCV_LSD )
if( ENABLE_OPENCV_LSD )
  add_definitions( -DENABLE_OPENCV_LSD )
  add_executable( ${tool_EXE} imgtabdet.cpp iomat.cpp imgio_magick.cpp log.cpp mem.cpp imgproc_cv.cpp )
else( ENABLE_OPENCV_LSD )
  add_executable( ${tool_EXE} imgtabdet.cpp iomat.cpp imgio_magick.cpp log.cpp mem.cpp )
endif( ENABLE_OPENCV_LSD )

add_executable( lsd LSD/lsd_cmd.c LSD/lsd.c )

string( REPLACE ";" " " CFLAGS_STR "-Wall -W ${MagickCore_CFLAGS} ${opencv_CFLAGS}" )
set_target_properties( ${tool_EXE} PROPERTIES COMPILE_FLAGS "${CFLAGS_STR}" )

include_directories( SYSTEM ${MagickCore_INCLUDEDIR} ) # To suppress system header warnings

target_link_libraries( ${tool_EXE} ${MagickCore_LDFLAGS} ${opencv_LDFLAGS} z m )

install( TARGETS ${tool_EXE} lsd DESTINATION bin )

add_custom_target( realclean cd ${CMAKE_HOME_DIRECTORY} COMMAND rm -fr ${tool_EXE} ${tool_EXE}.exe CMakeFiles CMakeCache.txt cmake_install.cmake install_manifest.txt Makefile )
