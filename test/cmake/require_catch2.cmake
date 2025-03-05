
cmake_minimum_required (VERSION 3.20)

# Exposes the catch2 target
# CMake 3.24: https://cmake.org/cmake/help/latest/module/FetchContent.html#variable:FETCHCONTENT_TRY_FIND_PACKAGE_MODE
macro (require_catch2)
    enable_testing ()

    set (catch2_version "3.4.0")

    find_package (Catch2 ${catch2_version} EXACT QUIET)

    if (NOT Catch2_FOUND)
        message (STATUS "Fetching Catch2 ${catch2_version}")

        Include(FetchContent)

        FetchContent_Declare(
          Catch2
          GIT_REPOSITORY https://github.com/catchorg/Catch2.git
          GIT_TAG        "v${catch2_version}" # or a later release
        )

        FetchContent_MakeAvailable(Catch2)
    else ()
        message (STATUS "Found Catch2 ${catch2_version}")
    endif ()

endmacro ()
