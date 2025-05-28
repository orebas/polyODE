vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO algebraic-solving/msolve
    REF "v${VERSION}" # Uses the version from vcpkg.json
    SHA512 "50b92227e051211b75256fa9ec35ef9069d635f16c027db0d5ff57aeddebb26ec2ab8f56208fb59a18b53a2c109da24971c3d112844cb92faba2f2509e6991c1"
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

message(STATUS "[DEBUG] Contents of ${CURRENT_PACKAGES_DIR}/bin after vcpkg_install_make:")
file(GLOB DEBUG_BIN_FILES "${CURRENT_PACKAGES_DIR}/bin/*")
foreach(FILE ${DEBUG_BIN_FILES})
    message(STATUS "  ${FILE}")
endforeach()
message(STATUS "[DEBUG] Contents of ${CURRENT_PACKAGES_DIR}/tools after vcpkg_install_make:")
file(GLOB DEBUG_TOOLS_FILES "${CURRENT_PACKAGES_DIR}/tools/*")
foreach(FILE ${DEBUG_TOOLS_FILES})
    message(STATUS "  ${FILE}")
endforeach()
message(STATUS "[DEBUG] Contents of ${CURRENT_PACKAGES_DIR} after vcpkg_install_make:")
file(GLOB_RECURSE DEBUG_ALL_FILES LIST_DIRECTORIES true "${CURRENT_PACKAGES_DIR}/*")
foreach(FILE ${DEBUG_ALL_FILES})
    message(STATUS "  ${FILE}")
endforeach()

# The msolve binary is installed to tools/msolve/bin/msolve by the Makefile
# Copy it to the tools directory structure expected by vcpkg
file(MAKE_DIRECTORY "${CURRENT_PACKAGES_DIR}/tools/msolve")
file(COPY "${CURRENT_PACKAGES_DIR}/tools/msolve/bin/msolve" 
     DESTINATION "${CURRENT_PACKAGES_DIR}/tools/msolve"
     FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

# Clean up the extra bin directory
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/tools/msolve/bin")
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/tools/msolve/debug")

# If msolve also installs shared libraries to lib that are needed by the executable,
# vcpkg's automatic fixup should handle making them available.
# We might need to also copy them to tools/${PORT}/lib if they are not found otherwise at runtime.
# For now, let's assume the msolve executable is statically linked or system libraries are sufficient.

vcpkg_fixup_pkgconfig()

# Install license file
file(INSTALL "${SOURCE_PATH}/COPYING" DESTINATION "${CURRENT_PACKAGES_DIR}/share/${PORT}" RENAME copyright)

# Remove unnecessary files from include directory if any (e.g. .git files if they sneak in)
# file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/include/.git")

# Handle CMake integration if msolve provided its own CMake files (it doesn't seem to, relies on pkg-config)
# vcpkg_cmake_config_fixup(PACKAGE_NAME msolve CONFIG_PATH "lib/cmake/msolve") # If it had CMake configs 