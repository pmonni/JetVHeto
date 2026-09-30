set_directory_properties (PROPERTIES CLEAN_NO_CUSTOM TRUE)
set (version "1.2")
set (filename chaplin-${version}.tar)
set (url https://chaplin.hepforge.org/code/${filename})
set (md5 d2e944f2011f8c1d8a1d00a487ad9d6f)
set (download_dir ${CMAKE_CURRENT_SOURCE_DIR})
set (source_file ${CMAKE_CURRENT_SOURCE_DIR}/${filename})
set (source_dir ${CMAKE_CURRENT_SOURCE_DIR}/chaplin-${version})

include (ExternalProject)

set (configure_command ./configure
  --prefix=${EXTERNAL_BUILD_DIR}
  --with-pic
  --disable-static
)

# this defines a target chaplin which does everything to build and install CHAPLIN
ExternalProject_Add (chaplin
  URL ${url}
  URL_MD5 ${md5}
  DOWNLOAD_DIR ${download_dir}
  SOURCE_DIR ${source_dir}
  CONFIGURE_COMMAND ${configure_command}
  BUILD_IN_SOURCE 1
  BUILD_COMMAND make
  INSTALL_COMMAND make install
)
