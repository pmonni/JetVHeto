set_directory_properties (PROPERTIES CLEAN_NO_CUSTOM TRUE)
set (version "1.2.0")
set (filename hoppet-${version}.tgz)
set (url https://hoppet.hepforge.org/downloads/${filename})
set (md5 6fe7c9d6da12283ba825b913f79a331e)
set (download_dir ${CMAKE_CURRENT_SOURCE_DIR})
set (source_file ${CMAKE_CURRENT_SOURCE_DIR}/${filename})
set (source_dir ${CMAKE_CURRENT_SOURCE_DIR}/hoppet-${version})

include (ExternalProject)

# this defines a target hoppet which does everything to build and install HOPPET
ExternalProject_Add (hoppet
  URL ${url}
  URL_MD5 ${md5}
  DOWNLOAD_DIR ${download_dir}
  SOURCE_DIR ${source_dir}
  CONFIGURE_COMMAND ./configure --prefix=${EXTERNAL_BUILD_DIR}
  BUILD_IN_SOURCE 1
  BUILD_COMMAND make
  INSTALL_COMMAND make install
)
