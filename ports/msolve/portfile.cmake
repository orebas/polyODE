vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO algebraic-solving/msolve
    REF "v${VERSION}" # Uses the version from vcpkg.json
    SHA512 "1fc8c1cc61d8d67ce0e0b881422b6c2a62b4a9b309105d3e4eaedb9974e7a66699eacc4baf1cdbcc5262a60728b0d42868ef4d30ae3b6ec6a912900d8562bba6"
    HEAD_REF master # or main, or a specific development branch if needed for HEAD builds
)

# Options for ./configure
# Define CONFIGURE_OPTIONS as a list of options
set(CONFIGURE_OPTIONS_LIST "--enable-shared") # Initial option
list(APPEND CONFIGURE_OPTIONS_LIST "--with-flint=${CURRENT_INSTALLED_DIR}") # Append as a new element

# Forcefully add FLINT's include path to CFLAGS and CPPFLAGS
# This will be seen by both the configure script and the make compilation.
set(ENV{CFLAGS} "$ENV{CFLAGS} -I${CURRENT_INSTALLED_DIR}/include")
set(ENV{CPPFLAGS} "$ENV{CPPFLAGS} -I${CURRENT_INSTALLED_DIR}/include")
set(ENV{LDFLAGS} "$ENV{LDFLAGS} -L${CURRENT_INSTALLED_DIR}/lib")

if(VCPKG_TARGET_IS_WINDOWS)
    # msolve is GCC/Clang only. On Windows, it's typically built with MinGW.
    # vcpkg might handle this by using a MinGW toolchain for this port if one is available and configured.
    # If not, this port might fail on Windows with MSVC.
    # Forcing a specific CFLAGS for MinGW might be needed here if default vcpkg flags are insufficient.
    # set(ENV{CFLAGS} "$ENV{CFLAGS} -O3") # Example if needed, vcpkg usually handles optimization flags.
    message(WARNING "msolve port on Windows might require a MinGW toolchain or specific setup. It is primarily a GCC/Clang library.")
endif()


vcpkg_configure_make(
    SOURCE_PATH ${SOURCE_PATH}
    AUTOCONFIG # This will run ./autogen.sh first, then ./configure
    OPTIONS ${CONFIGURE_OPTIONS_LIST} # Pass the list of options
)

vcpkg_install_make()

vcpkg_fixup_pkgconfig()

# Install license file
file(INSTALL "${SOURCE_PATH}/COPYING" DESTINATION "${CURRENT_PACKAGES_DIR}/share/${PORT}" RENAME copyright)

# Remove unnecessary files from include directory if any (e.g. .git files if they sneak in)
# file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/include/.git")

# Handle CMake integration if msolve provided its own CMake files (it doesn't seem to, relies on pkg-config)
# vcpkg_cmake_config_fixup(PACKAGE_NAME msolve CONFIG_PATH "lib/cmake/msolve") # If it had CMake configs 