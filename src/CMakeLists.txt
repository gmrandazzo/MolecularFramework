set(MOLECULARFRAMEWORK_C_SRCS version.c
                   misc.c
                   periodic_table.c
                   massanalysis.c
                   atomanalysis.c
                   protanalysis.c
                   atomdesc.c
                   forcefield.c
                   fielddesc.c
                   geomdesc.c
                   shapedesc.c
                   miscdesc.c
                   molecule.c
                   wfingerprint.c
                   mol3Daligner.c
                   mol3Dfields.c)

set(MOLECULARFRAMEWORK_C_H version.h
                molecule.h
                misc.h
                periodic_table.h
                massanalysis.h
                atomanalysis.h
                protanalysis.h
                atomdesc.h
                forcefield.h
                fielddesc.h
                geomdesc.h
                shapedesc.h
                miscdesc.h
                wfingerprint.h)

set(CPPWRAP_CPP_SRC cppwrap.cpp
                    ringsearch.cpp
                    graphanalysis.cpp)

set(CPPWRAP_CPP_H cppwrap.h
                  ringsearch.h
                  graphanalysis.h)

if(WIN32)
  message(STATUS "Set mingw options")
  set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS_RELEASE)
  set(CMAKE_EXECUTABLE_SUFFIX ".exe")
  set(CMAKE_STATIC_LIBRARY_SUFFIX ".dll")
endif()

add_library(libcppwrap SHARED ${CPPWRAP_CPP_SRC} ${CPPWRAP_CPP_H})
#target_link_libraries(libcppwrap m ${RDKIT_LIBRARIES})
target_link_libraries(libcppwrap m)

add_library(libmoldesc SHARED ${MOLECULARFRAMEWORK_C_SRCS} ${MOLDESC_C_H})
set_target_properties(libmoldesc PROPERTIES OUTPUT_NAME "moldesc")
set_property(TARGET libmoldesc PROPERTY LINKER_LANGUAGE C)
target_link_libraries(libmoldesc m ${LIBSCIENTIFIC_LIBRARY} libcppwrap)

install(TARGETS libmoldesc LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
                      ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
                      RUNTIME DESTINATION bin)

install(FILES ${MOLECULARFRAMEWORK_C_H} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/moldesc)

add_subdirectory(tests)
add_subdirectory(bin_src)
