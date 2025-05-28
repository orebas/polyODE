cmake_minimum_required(VERSION 3.13)

project(PolyODE LANGUAGES CXX)

# Set newer CMake policies
cmake_policy(SET CMP0135 NEW) # Handle FetchContent DOWNLOAD_EXTRACT_TIMESTAMP

# Generate compile_commands.json for clang-tidy and other tools
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- Sanitizers (ASan/UBSan) ---
option(ENABLE_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer" ON)
if(ENABLE_SANITIZERS)
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        message(STATUS "Sanitizers enabled.")
        set(SANITIZER_FLAGS "-fsanitize=address,undefined")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SANITIZER_FLAGS}")
        # Add to linker flags as well, important for ASan
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${SANITIZER_FLAGS}")
        set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${SANITIZER_FLAGS}")
    else()
        message(WARNING "Sanitizers requested but only supported on GCC/Clang.")
    endif()
endif()

# --- Static Analysis (Clang-Tidy) ---
find_program(CLANGTIDY_EXE clang-tidy)
if(CLANGTIDY_EXE)
    message(STATUS "clang-tidy found: ${CLANGTIDY_EXE}. Use 'check-tidy' target.")
    # Remove the automatic execution during build:
    # set(CMAKE_CXX_CLANG_TIDY "${CLANGTIDY_EXE}")

    # Define a custom target to run clang-tidy on demand using the script
    add_custom_target(check-tidy
        # Use the shell script to run clang-tidy
        # Pass PATH using cmake -E env to ensure clang-tidy is found by the script if needed
        COMMAND ${CMAKE_COMMAND} -E env PATH="$ENV{PATH}" ${CMAKE_CURRENT_SOURCE_DIR}/run_clang_tidy.sh
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} # Run script from project root
        COMMENT "Running clang-tidy static analysis via script (report in clang-tidy-report.txt)"
        VERBATIM
    )
    # No explicit dependency needed here as the script handles finding compile_commands.json itself.
    # The script generates the report file in the directory it's run from (project root).

    # Diagnostic target removed.

else()
    message(WARNING "clang-tidy not found. Static analysis target 'check-tidy' will not be available.")
    # Define a dummy target if clang-tidy is not found
    add_custom_target(check-tidy
        COMMAND ${CMAKE_COMMAND} -E echo "clang-tidy not found, skipping static analysis."
    )
    # Dummy diagnostic target removed.
endif()

# Include directories
# Add the project's include directory
include_directories(include)

# Find dependencies using find_package. vcpkg integration relies on this.
find_package(Boost REQUIRED) # vcpkg handles finding this via toolchain
find_package(Ceres REQUIRED) # vcpkg handles finding this via toolchain
find_package(PkgConfig REQUIRED) # Needed to find msolve via its .pc file

# Set PKG_CONFIG_PATH to include vcpkg installed packages
if(DEFINED VCPKG_INSTALLED_DIR)
    set(ENV{PKG_CONFIG_PATH} "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
endif()

pkg_check_modules(MSOLVE REQUIRED IMPORTED_TARGET msolve) # vcpkg ensures msolve.pc is findable

# Find the msolve executable
find_program(MSOLVE_EXECUTABLE msolve
    HINTS "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/tools/msolve" # Where the port actually installs it
          "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/tools/msolve/bin" # Original hint for vcpkg
          "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin" # Fallback if installed to bin
    PATHS ENV PATH # Check system PATH as well
    DOC "Path to msolve executable"
)

if(NOT MSOLVE_EXECUTABLE)
    message(FATAL_ERROR "MSolve executable (msolve) not found. Ensure it is installed correctly by vcpkg or is in PATH.")
else()
    message(STATUS "Found msolve executable: ${MSOLVE_EXECUTABLE}")
endif()

# Find FLINT library (used by univariate_solver_flint) - using pkg-config with vcpkg FLINT
pkg_check_modules(FLINT REQUIRED IMPORTED_TARGET flint)

# GTest will be found here now instead of FetchContent
# find_package(GTest REQUIRED) # Moved lower, near testing section

# Remove global include directories for dependencies managed by vcpkg.
# target_link_libraries with imported targets (e.g., Boost::boost, Ceres::ceres)
# should propagate include directories.
# include_directories(${Boost_INCLUDE_DIRS}) # Removed
# include_directories(${CERES_INCLUDE_DIRS}) # Removed

# --- Library Target ---
# Define the library using sources from src/
add_library(polyode_lib STATIC
    src/polynomial.cpp
    src/parameter_estimation.cpp
    src/parameter_estimator.cpp
    src/identifiability_analyzer.cpp
    src/approximation/observable_approximator.cpp # Base class (might be header-only, but good practice)
    src/approximation/aa_approximator.cpp       # AAA implementation
    src/phc_solver.cpp                          # PHC solver implementation
    src/ceres_algebraic_solver.cpp              # Added CeresAlgebraicSolver
    src/MSolveSolver.cpp                        # MSolve solver implementation
    src/example_systems.cpp                     # Added to compile and link with the library
    src/univariate_solver_flint.cpp             # Added univariate solver implementation
    src/solution_polisher.cpp                   # Added SolutionPolisher
    tests/model_test_framework.cpp              # Test framework for systematic model testing
    tests/model_registrations.cpp               # Model registrations for test framework
)
# Specify that the library includes headers from the include/ directory.
# PUBLIC makes these includes available to targets linking against polyode_lib.
target_include_directories(polyode_lib PUBLIC include)
# Add tests directory for test_utils.hpp needed by parameter_estimation.hpp
# PRIVATE because users of the library shouldn't need test utils.
target_include_directories(polyode_lib PUBLIC tests)
# Add build directory include for generated headers like msolve_config.hpp
target_include_directories(polyode_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/include)

# Generate config header with msolve executable path
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/msolve_config.hpp.in
               ${CMAKE_CURRENT_BINARY_DIR}/include/msolve_config.hpp @ONLY)

# Define __GMP_DONT_USE_CXX_STREAM_OPS for polyode_lib and any other
# targets that will compile files including FLINT headers.
# set(GMP_CXX_GUARD __GMP_DONT_USE_CXX_STREAM_OPS=1) # REMOVED
# target_compile_definitions(polyode_lib PRIVATE ${GMP_CXX_GUARD}) # REMOVED

# Link Boost if the library itself needs compiled Boost components (unlikely for header-only like Odeint)
# Note: Boost::boost is a common interface target for header-only usage.
# If specific compiled components are needed, list them in find_package and link e.g., Boost::system
# target_link_libraries(polyode_lib PUBLIC Ceres::ceres PkgConfig::MSOLVE FLINT::flint) # OLD from previous attempt
target_link_libraries(polyode_lib PUBLIC Ceres::ceres PkgConfig::MSOLVE PkgConfig::FLINT) # REVERTED to PkgConfig::FLINT

# --- Example Targets ---
# Example 1: Lotka-Volterra
add_executable(lotka_volterra_example examples/lotka_volterra.cpp)
target_link_libraries(lotka_volterra_example PRIVATE polyode_lib Boost::boost)
if(UNIX AND NOT APPLE)
    target_link_libraries(lotka_volterra_example PRIVATE m)
endif()

# Example 2: Holling
add_executable(holling_example examples/holling_test.cpp)
target_link_libraries(holling_example PRIVATE polyode_lib Boost::boost)
if(UNIX AND NOT APPLE)
    target_link_libraries(holling_example PRIVATE m)
endif()

# Example 3: Jacobian Test
add_executable(jacobian_test examples/jacobian_test.cpp)
target_link_libraries(jacobian_test PRIVATE polyode_lib Boost::boost Ceres::ceres)
if(UNIX AND NOT APPLE)
    target_link_libraries(jacobian_test PRIVATE m)
endif()

# Example 4: Basic Parameter Estimation Test
add_executable(basic_estimation examples/basic_estimation.cpp)
target_link_libraries(basic_estimation PRIVATE polyode_lib Boost::boost Ceres::ceres)
if(UNIX AND NOT APPLE)
    target_link_libraries(basic_estimation PRIVATE m)
endif()

# Example 5: Estimate IC and Param Test
add_executable(estimate_ic_param examples/estimate_ic_param.cpp)
target_link_libraries(estimate_ic_param PRIVATE polyode_lib Boost::boost Ceres::ceres)
if(UNIX AND NOT APPLE)
    target_link_libraries(estimate_ic_param PRIVATE m)
endif()

# Example 6: Blinded Parameter Estimation Test
add_executable(blinded_estimation examples/blinded_estimation.cpp)
target_link_libraries(blinded_estimation PRIVATE polyode_lib Boost::boost Ceres::ceres)
if(UNIX AND NOT APPLE)
    target_link_libraries(blinded_estimation PRIVATE m)
endif()

# Example 7: Difficult Estimation Scenario
add_executable(difficult_estimation examples/difficult_estimation.cpp)
target_link_libraries(difficult_estimation PRIVATE polyode_lib Boost::boost Ceres::ceres)
if(UNIX AND NOT APPLE)
    target_link_libraries(difficult_estimation PRIVATE m)
endif()

# Example 8: Complex System Identifiability Analysis
add_executable(complex_system_identifiability examples/complex_system_identifiability.cpp)
target_link_libraries(complex_system_identifiability PRIVATE polyode_lib Boost::boost Ceres::ceres)
if(UNIX AND NOT APPLE)
    target_link_libraries(complex_system_identifiability PRIVATE m)
endif()

# Add copy of PHC script to build directory (before defining the tests)
# This ensures the script is available at runtime from the build dir
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/scripts/phc_dict_to_json.py
               ${CMAKE_CURRENT_BINARY_DIR}/scripts/phc_dict_to_json.py
               COPYONLY)

# Make sure the script is executable
file(CHMOD ${CMAKE_CURRENT_BINARY_DIR}/scripts/phc_dict_to_json.py
     PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

# --- Testing Framework (Google Test) ---
# Remove FetchContent block for googletest
# include(FetchContent)
# FetchContent_Declare(
#     googletest
#     # Specify the repository and tag (version)
#     URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
# )
# # For Windows: Prevent overriding the parent project's compiler/linker settings
# # This is often needed when integrating gtest via FetchContent
# set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
# FetchContent_MakeAvailable(googletest)
#
# # Disable warnings-as-errors for Google Test targets (no longer needed)
# # set_target_properties(gtest gmock PROPERTIES COMPILE_WARNING_AS_ERROR FALSE)

# Find GTest package provided by vcpkg
find_package(GTest REQUIRED)

# Enable CTest support
enable_testing()

# Include GoogleTest CMake module (provides gtest_discover_tests)
# This should work with the find_package version as well
include(GoogleTest)

# --- Test Targets ---
# Helper function to simplify adding test executables and discovering tests
function(add_polyode_test test_name)
    add_executable(${test_name} tests/${test_name}.cpp)
    target_link_libraries(${test_name} PRIVATE polyode_lib GTest::gtest_main)
    if(UNIX AND NOT APPLE)
        target_link_libraries(${test_name} PRIVATE m)
    endif()
    # Use modern gtest_discover_tests instead of gtest_add_tests
    gtest_discover_tests(${test_name})
endfunction()

# Add each test using the helper function
add_polyode_test(analytic_test)
add_polyode_test(holling_odesystem_test)
add_polyode_test(monomial_test)
add_polyode_test(polynomial_test)
add_polyode_test(rational_function_test)
add_polyode_test(variable_test)
add_polyode_test(ode_system_test)
add_polyode_test(exponential_decay_test)
# add_polyode_test(parameter_estimation_test) # Removed old test
add_polyode_test(identifiability_analyzer_test)
add_polyode_test(approximation_test)
add_polyode_test(univariate_solver_flint_test) # Added univariate solver test

# Only link to FLINT for the univariate solver test
# target_link_libraries(univariate_solver_flint_test PRIVATE FLINT::flint) # OLD from previous attempt
target_link_libraries(univariate_solver_flint_test PRIVATE PkgConfig::FLINT) # REVERTED to PkgConfig::FLINT
# Also add the compile definition to the test target if it includes FLINT directly
# target_compile_definitions(univariate_solver_flint_test PRIVATE ${GMP_CXX_GUARD}) # REMOVED

# Add the new parameter estimator test using the helper
add_polyode_test(parameter_estimator_test)
add_polyode_test(parameter_estimator_scenarios_test)
add_polyode_test(systematic_model_tests)
# Link specific libraries needed by this test (e.g., for PHCSolver)
# Assuming vcpkg provides nlohmann_json::nlohmann_json target
find_package(nlohmann_json REQUIRED)
# Link PHCSolver dependencies if needed (nlohmann_json)
# Check if polyode_lib already links nlohmann_json PUBLICLY.
# If not, link it directly to the test target.
# Let's assume it's not linked publicly by the lib for now.
target_link_libraries(parameter_estimator_test PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(parameter_estimator_scenarios_test PRIVATE nlohmann_json::nlohmann_json)
# Link Ceres if ParameterEstimator directly uses Ceres types (it doesn't anymore)
# target_link_libraries(parameter_estimator_test PRIVATE Ceres::ceres)

# --- Compiler Warnings ---
# Apply common warning flags to the library target
# INTERFACE ensures these options are used when compiling the library
# and also propagated to targets linking against it (like examples and tests)
if(MSVC)
    target_compile_options(polyode_lib INTERFACE /W4) # Example for MSVC
else()
    target_compile_options(polyode_lib INTERFACE -Wall -Wextra -pedantic) # GCC/Clang
endif()

# --- Doxygen Documentation ---
find_package(Doxygen)
if(DOXYGEN_FOUND)
    # Set Doxygen variables for the template
    set(DOXYGEN_PROJECT_NAME ${PROJECT_NAME})
    # Use CMAKE_PROJECT_VERSION if set, otherwise fallback
    if(NOT DEFINED CMAKE_PROJECT_VERSION)
        set(PROJECT_VERSION "0.1.0") # Default version if not set
    else()
        set(PROJECT_VERSION ${CMAKE_PROJECT_VERSION})
    endif()

    # Configure the Doxyfile using the template
    set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
    set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)
    message(STATUS "Doxygen build enabled. Configured Doxyfile: ${DOXYFILE_OUT}")

    # Add a target to generate documentation
    # Use the new doxygen_add_docs function (CMake 3.9+)
    # Note: This function might implicitly depend on the input source files,
    # but explicitly adding dependencies can be safer if needed.
    doxygen_add_docs(
        doc # Target name (e.g., build with 'make doc')
        ${DOXYFILE_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} # Doxygen runs from build dir
        COMMENT "Generating API documentation with Doxygen to docs/html"
    )
    # If you want to ensure 'doc' depends on your library sources (optional):
    # get_target_property(polyode_lib_sources polyode_lib SOURCES)
    # add_dependencies(doc ${polyode_lib_sources})
else()
    message(WARNING "Doxygen not found. Documentation target 'doc' will not be available.")
    # Define a dummy target if Doxygen is not found
    add_custom_target(doc
        COMMAND ${CMAKE_COMMAND} -E echo "Doxygen not found, skipping documentation generation."
    )
endif()

# --- Installation (Optional - Example structure) ---
# install(TARGETS polyode_lib DESTINATION lib)
# install(DIRECTORY include/ DESTINATION include) # Installs the entire include directory
# install(TARGETS lotka_volterra_example holling_example DESTINATION bin) # Install example executables


# --- Old Targets/CTest commands (Removed for cleanup) ---
