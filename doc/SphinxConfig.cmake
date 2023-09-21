##############################################################################
# @file  FindSphinx.cmake
# @brief Find Sphinx documentation build tools.
#
# @par Input variables:
# <table border="0">
#   <tr>
#     @tp @b Sphinx_DIR @endtp
#     <td>Installation directory of Sphinx tools. Can also be set as environment variable.</td>
#   </tr>
#   <tr>
#     @tp @b SPHINX_DIR @endtp
#     <td>Alternative environment variable for @c Sphinx_DIR.</td>
#   </tr>
#   <tr>
#     @tp @b Sphinx_FIND_COMPONENTS @endtp
#     <td>Sphinx build tools to look for, i.e., 'apidoc' and/or 'build'.</td>
#   </tr>
# </table>
#
# @par Output variables:
# <table border="0">
#   <tr>
#     @tp @b Sphinx_FOUND @endtp
#     <td>Whether all or only the requested Sphinx build tools were found.</td>
#   </tr>
#   <tr>
#     @tp @b SPHINX_FOUND @endtp
#     <td>Alias for @c Sphinx_FOUND.<td>
#   </tr>
#   <tr>
#     @tp @b SPHINX_EXECUTABLE @endtp
#     <td>Non-cached alias for @c Sphinx-build_EXECUTABLE.</td>
#   </tr>
#   <tr>
#     @tp @b Sphinx_PYTHON_EXECUTABLE @endtp
#     <td>Python executable used to run sphinx-build. This is either the
#         by default found Python interpreter or a specific version as
#         specified by the shebang (#!) of the sphinx-build script.</td>
#   </tr>
#   <tr>
#     @tp @b Sphinx_PYTHON_OPTIONS @endtp
#     <td>A list of Python options extracted from the shebang (#!) of the
#         sphinx-build script. The -E option is added by this module
#         if the Python executable is not the system default to avoid
#         problems with a differing setting of the @c PYTHONHOME.</td>
#   </tr>
#   <tr>
#     @tp @b Sphinx-build_EXECUTABLE @endtp
#     <td>Absolute path of the found sphinx-build tool.</td>
#   </tr>
#   <tr>
#     @tp @b Sphinx-apidoc_EXECUTABLE @endtp
#     <td>Absolute path of the found sphinx-apidoc tool.</td>
#   </tr>
#   <tr>
#     @tp @b Sphinx_VERSION_STRING @endtp
#     <td>Sphinx version found e.g. 1.1.2.</td>
#   </tr>
#   <tr>
#     @tp @b Sphinx_VERSION_MAJOR @endtp
#     <td>Sphinx major version found e.g. 1.</td>
#   </tr>
#   <tr>
#     @tp @b Sphinx_VERSION_MINOR @endtp
#     <td>Sphinx minor version found e.g. 1.</td>
#   </tr>
#   <tr>
#     @tp @b Sphinx_VERSION_PATCH @endtp
#     <td>Sphinx patch version found e.g. 2.</td>
#   </tr>
# </table>
#
# @ingroup CMakeFindModules
##############################################################################

#=============================================================================
# Copyright 2011-2012 University of Pennsylvania
# Copyright 2013-2016 Andreas Schuh <andreas.schuh.84@gmail.com>
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
include(FindPackageHandleStandardArgs)

# We are likely to find Sphinx near the Python interpreter
find_package(Python COMPONENTS Interpreter Development)
if(Python_Interpreter_FOUND)
    get_filename_component(_PYTHON_DIR "${Python_EXECUTABLE}" DIRECTORY)
    set(
        _PYTHON_PATHS
        "${_PYTHON_DIR}"
        "${_PYTHON_DIR}/bin"
        "${_PYTHON_DIR}/Scripts")
endif()

find_program(
    SPHINX_EXECUTABLE
    NAMES sphinx-build
    HINTS ${_PYTHON_PATHS})
mark_as_advanced(SPHINX_EXECUTABLE)

find_package_handle_standard_args(Sphinx DEFAULT_MSG SPHINX_EXECUTABLE)
