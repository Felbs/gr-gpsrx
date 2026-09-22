find_package(PkgConfig)

PKG_CHECK_MODULES(PC_GR_GPSRX gnuradio-gpsrx)

FIND_PATH(
    GR_GPSRX_INCLUDE_DIRS
    NAMES gnuradio/gpsrx/api.h
    HINTS $ENV{GPSRX_DIR}/include
        ${PC_GPSRX_INCLUDEDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/include
          /usr/local/include
          /usr/include
)

FIND_LIBRARY(
    GR_GPSRX_LIBRARIES
    NAMES gnuradio-gpsrx
    HINTS $ENV{GPSRX_DIR}/lib
        ${PC_GPSRX_LIBDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/lib
          ${CMAKE_INSTALL_PREFIX}/lib64
          /usr/local/lib
          /usr/local/lib64
          /usr/lib
          /usr/lib64
          )

include("${CMAKE_CURRENT_LIST_DIR}/gnuradio-gpsrxTarget.cmake")

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GR_GPSRX DEFAULT_MSG GR_GPSRX_LIBRARIES GR_GPSRX_INCLUDE_DIRS)
MARK_AS_ADVANCED(GR_GPSRX_LIBRARIES GR_GPSRX_INCLUDE_DIRS)
