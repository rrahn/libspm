cmake_minimum_required (VERSION 3.20)

# Find doxygen.
find_package (Doxygen QUIET)

if (${DOXYGEN_FOUND})
    message (STATUS "Found Doxygen: ${DOXYGEN_EXECUTABLE}")

    # Configure doxygen options.
    set (LIBJST_DOC_SOURCE_DIR ${CMAKE_SOURCE_DIR}/doc)
    set (LIBJST_DOC_BINARY_DIR ${PROJECT_BINARY_DIR}/doc)

    set (APP_TEMPLATE_DOXYGEN_EXCLUDE_SYMBOLS "")
    set (APP_TEMPLATE_DOXYGEN_PREDEFINED_NDEBUG "")
    set (APP_TEMPLATE_DOXYGEN_ENABLED_SECTIONS "DEV")
    set (APP_TEMPLATE_DOXYGEN_EXTRACT_PRIVATE "YES")
    set (JSTMAP_DOXYGEN_LOGFILE "${LIBJST_DOC_BINARY_DIR}/doxygen.log")

    configure_file (${LIBJST_DOC_SOURCE_DIR}/sphinx_doxygen_cfg ${LIBJST_DOC_BINARY_DIR}/build/Doxyfile)
    configure_file (${LIBJST_DOC_SOURCE_DIR}/conf.py.in ${LIBJST_DOC_BINARY_DIR}/source/conf.py)
    configure_file (${LIBJST_DOC_SOURCE_DIR}/index.rst.in ${LIBJST_DOC_BINARY_DIR}/source/index.rst)
    FILE(MAKE_DIRECTORY ${LIBJST_DOC_BINARY_DIR}/source/_static)
    FILE(MAKE_DIRECTORY ${LIBJST_DOC_BINARY_DIR}/source/_template)

    add_custom_target (docs_base
                       COMMAND ${DOXYGEN_EXECUTABLE}
                       WORKING_DIRECTORY ${LIBJST_DOC_BINARY_DIR}/build
                       COMMENT "Generating API documentation with Doxygen."
                       VERBATIM)

    message (STATUS "${FontBold}You can run `make docs_base` to build api documentation.${FontReset}")

    find_package (Sphinx HINTS "${CMAKE_SOURCE_DIR}/doc" NO_DEFAULT_PATH)
    if (SPHINX_FOUND)
        message (STATUS "Found Sphinx: ${SPHINX_EXECUTABLE}")

        # Configure sphinx options.
        add_custom_target (docs_sphinx
                           COMMAND ${SPHINX_EXECUTABLE} -b html ${LIBJST_DOC_BINARY_DIR}/source ${LIBJST_DOC_BINARY_DIR}/build
                           WORKING_DIRECTORY ${LIBJST_DOC_BINARY_DIR}
                           COMMENT "Generating documentation with Sphinx."
                           VERBATIM
                           DEPENDS docs_base)
        message (STATUS "${FontBold}You can run `make docs_sphinx` to build html documentation.${FontReset}")
    else ()
        message (STATUS "Sphinx not found.")
    endif ()
else ()
    message (STATUS "Doxygen not found.")
endif ()

### Enable testing

enable_testing ()

message (STATUS "Add doc test.")
add_test (NAME doc_test
          COMMAND bash -c "cat ${JSTMAP_DOXYGEN_LOGFILE}; test ! -s ${JSTMAP_DOXYGEN_LOGFILE}"
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
