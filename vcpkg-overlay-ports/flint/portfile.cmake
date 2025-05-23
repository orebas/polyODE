set(VERSION 3.2.1) # Define version for clarity

vcpkg_download_distfile(ARCHIVE # This variable will be set to the path of the downloaded file
    URLS "https://github.com/flintlib/flint/releases/download/v${VERSION}/flint-${VERSION}.tar.gz"
    FILENAME "flint-${VERSION}.tar.gz"
    SHA256 ca7be46d77972277eb6fe0c4f767548432f56bb534aa17d6dba2d7cce15cd23f
    SHA512 3d23e78fb57f0ce95d30fef474063b7a92565e474200d9aba65d8e2fcebe80f557acaabb936845a2bd138e2b275c1258595a6e0538f474050e772ccfc1d72e5b
)

# Extracts the archive. SOURCE_PATH will be set to the root of the extracted sources
# (e.g., where FLINT's CMakeLists.txt is located).
vcpkg_extract_source_archive(SOURCE_PATH
    ARCHIVE "${ARCHIVE}"
    # PATCHES # Add any patches if needed
)

# FLINT depends on GMP and MPFR.
# Vcpkg will handle these dependencies automatically as they are specified in vcpkg-overlay-ports/flint/vcpkg.json

# Options for FLINT's ./configure script
# Vcpkg sets environment variables like CFLAGS, LDFLAGS, PKG_CONFIG_PATH
# to help find dependencies (GMP, MPFR).
# --with-gmp and --with-mpfr can often be omitted if vcpkg's environment is sufficient,
# but explicitly passing them using VCPKG_CROSSCOMPILING information can be more robust.
# We'll let vcpkg's environment try to handle it first.

# Explicitly pass standard autotools directory layout options.
# vcpkg_configure_make will also pass --prefix=${CURRENT_PACKAGES_DIR}
# These should tell FLINT's configure script to generate Makefiles that install
# into standard subdirectories of the prefix.
set(CONFIGURE_OPTIONS
    "--disable-tests"
    "--includedir=\${prefix}/include"  # Headers go into $PREFIX/include
    "--libdir=\${prefix}/lib"          # Libraries go into $PREFIX/lib
    # For pkg-specific headers like flint/*.h, the convention is often that the software's
    # own install rules create the 'flint' subdirectory under $(includedir).
    # So, if includedir is $PREFIX/include, headers should go to $PREFIX/include/flint.
)

# Let vcpkg_configure_make handle shared/static options based on VCPKG_LIBRARY_LINKAGE.
# We only need to add extra options not automatically handled.
vcpkg_configure_make(
    SOURCE_PATH "${SOURCE_PATH}"
    OPTIONS ${CONFIGURE_OPTIONS}
)

vcpkg_install_make()

# After make install, manually copy the generated headers from build tree
# The FLINT Makefile installs all headers but misses flint.h and flint-config.h which are generated

# Copy generated headers from build tree to the proper location
file(COPY "${CURRENT_BUILDTREES_DIR}/${TARGET_TRIPLET}-rel/src/flint.h" 
     DESTINATION "${CURRENT_PACKAGES_DIR}/include/flint/"
     FILE_PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)

file(COPY "${CURRENT_BUILDTREES_DIR}/${TARGET_TRIPLET}-rel/src/flint-config.h" 
     DESTINATION "${CURRENT_PACKAGES_DIR}/include/flint/"
     FILE_PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)

# Also copy to debug directory if it exists
if(EXISTS "${CURRENT_BUILDTREES_DIR}/${TARGET_TRIPLET}-dbg/src/flint.h")
    file(COPY "${CURRENT_BUILDTREES_DIR}/${TARGET_TRIPLET}-dbg/src/flint.h" 
         DESTINATION "${CURRENT_PACKAGES_DIR}/debug/include/flint/"
         FILE_PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)
endif()

if(EXISTS "${CURRENT_BUILDTREES_DIR}/${TARGET_TRIPLET}-dbg/src/flint-config.h")
    file(COPY "${CURRENT_BUILDTREES_DIR}/${TARGET_TRIPLET}-dbg/src/flint-config.h" 
         DESTINATION "${CURRENT_PACKAGES_DIR}/debug/include/flint/"
         FILE_PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)
endif()

# FLINT installs flint.pc, so fix it up for pkg-config users
# Fix the incorrect includedir path in flint.pc
vcpkg_replace_string("${CURRENT_PACKAGES_DIR}/lib/pkgconfig/flint.pc" "includedir=/include" "includedir=\${prefix}/include")
if(EXISTS "${CURRENT_PACKAGES_DIR}/debug/lib/pkgconfig/flint.pc")
    vcpkg_replace_string("${CURRENT_PACKAGES_DIR}/debug/lib/pkgconfig/flint.pc" "includedir=/include" "includedir=\${prefix}/include")
endif()
vcpkg_fixup_pkgconfig()

# Standard vcpkg cleanup.
# The debug/include will be removed if it's empty or only contains empty directories.
# If our moves above worked, debug/include/flint/ should exist and have files,
# so debug/include should not be entirely removed.
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/share") # Always safe

# Install the license file
file(INSTALL "${SOURCE_PATH}/COPYING.LESSER" DESTINATION "${CURRENT_PACKAGES_DIR}/share/${PORT}" RENAME copyright)
# Install usage instructions
configure_file("${CMAKE_CURRENT_LIST_DIR}/usage" "${CURRENT_PACKAGES_DIR}/share/${PORT}/usage" COPYONLY)

# Post-build validation (optional)
# vcpkg_test_cmake(PACKAGE_NAME FLINT) 