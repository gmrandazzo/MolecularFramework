# Install script for directory: /home/marco/Nextcloud/Software/MolecularFramework/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/usr/local/lib/libmoldesc.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/lib/libmoldesc.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/usr/local/lib/libmoldesc.so"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/lib/libmoldesc.so")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/usr/local/lib" TYPE SHARED_LIBRARY FILES "/home/marco/Nextcloud/Software/MolecularFramework/build/src/libmoldesc.so")
  if(EXISTS "$ENV{DESTDIR}/usr/local/lib/libmoldesc.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/lib/libmoldesc.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/usr/local/lib/libmoldesc.so"
         OLD_RPATH "/usr/local/lib:/home/marco/Nextcloud/Software/MolecularFramework/build/src:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/usr/local/lib/libmoldesc.so")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/include/moldesc/version.h;/usr/local/include/moldesc/molecule.h;/usr/local/include/moldesc/misc.h;/usr/local/include/moldesc/periodic_table.h;/usr/local/include/moldesc/massanalysis.h;/usr/local/include/moldesc/atomanalysis.h;/usr/local/include/moldesc/protanalysis.h;/usr/local/include/moldesc/atomdesc.h;/usr/local/include/moldesc/forcefield.h;/usr/local/include/moldesc/fielddesc.h;/usr/local/include/moldesc/geomdesc.h;/usr/local/include/moldesc/shapedesc.h;/usr/local/include/moldesc/miscdesc.h;/usr/local/include/moldesc/wfingerprint.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/usr/local/include/moldesc" TYPE FILE FILES
    "/home/marco/Nextcloud/Software/MolecularFramework/src/version.h"
    "/home/marco/Nextcloud/Software/MolecularFramework/src/molecule.h"
    "/home/marco/Nextcloud/Software/MolecularFramework/src/misc.h"
    "/home/marco/Nextcloud/Software/MolecularFramework/src/periodic_table.h"
    "/home/marco/Nextcloud/Software/MolecularFramework/src/massanalysis.h"
    "/home/marco/Nextcloud/Software/MolecularFramework/src/atomanalysis.h"
    "/home/marco/Nextcloud/Software/MolecularFramework/src/protanalysis.h"
    "/home/marco/Nextcloud/Software/MolecularFramework/src/atomdesc.h"
    "/home/marco/Nextcloud/Software/MolecularFramework/src/forcefield.h"
    "/home/marco/Nextcloud/Software/MolecularFramework/src/fielddesc.h"
    "/home/marco/Nextcloud/Software/MolecularFramework/src/geomdesc.h"
    "/home/marco/Nextcloud/Software/MolecularFramework/src/shapedesc.h"
    "/home/marco/Nextcloud/Software/MolecularFramework/src/miscdesc.h"
    "/home/marco/Nextcloud/Software/MolecularFramework/src/wfingerprint.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/marco/Nextcloud/Software/MolecularFramework/build/src/tests/cmake_install.cmake")
  include("/home/marco/Nextcloud/Software/MolecularFramework/build/src/bin_src/cmake_install.cmake")

endif()

