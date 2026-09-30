set_directory_properties (PROPERTIES CLEAN_NO_CUSTOM TRUE)
set (version 6.2.1)
set (filename LHAPDF-${version}.tar.gz)
set (url https://lhapdf.hepforge.org/downloads/${filename})
set (md5 9e05567d538fdb4862d4781cd076d7db)
set (download_dir ${CMAKE_CURRENT_SOURCE_DIR})
set (source_file ${CMAKE_CURRENT_SOURCE_DIR}/${filename})
set (source_dir ${CMAKE_CURRENT_SOURCE_DIR}/LHAPDF-${version})

include (ExternalProject)

set (configure_command ./configure
  --prefix=${EXTERNAL_BUILD_DIR}
)

# this defines a target lhapdf which does everything to build and install LHAPDF
ExternalProject_Add (lhapdf
  URL ${url}
  URL_MD5 ${md5}
  DOWNLOAD_DIR ${download_dir}
  SOURCE_DIR ${source_dir}
  CONFIGURE_COMMAND ${configure_command}
  BUILD_IN_SOURCE 1
  BUILD_COMMAND make
  INSTALL_COMMAND make install
)
