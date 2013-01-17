################################################################################
# TODO:  Generate license.html
################################################################################

find_package(Doxygen REQUIRED)

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/programmer/.config/Doxyfile.in
    ${CMAKE_CURRENT_SOURCE_DIR}/programmer/Doxyfile
    )

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/programmer/.config/programmer.html.in
    ${CMAKE_CURRENT_SOURCE_DIR}/programmer/programmer.html
    )

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/programmer/.config/mainpage.h.in
    ${CMAKE_CURRENT_SOURCE_DIR}/programmer/mainpage.h
    )

add_custom_target(
    doc
    ALL
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/programmer/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/programmer
    COMMENT "Generating documentation with Doxygen" VERBATIM
)

add_custom_target(
    doc-pdf
    ALL
    COMMAND make
    COMMAND cp refman.pdf ../programmer.pdf
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/programmer/latex
    COMMENT "Generating pdf documentation" VERBATIM
)
