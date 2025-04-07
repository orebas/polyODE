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

# Find Boost (needed by the library and potentially examples/tests)
find_package(Boost REQUIRED)
# Find Ceres (needed for jacobian_test example)
find_package(Ceres REQUIRED)

# Add Boost include directories globally for simplicity, or target_include_directories per target
include_directories(${Boost_INCLUDE_DIRS})
# Add Ceres include directory
include_directories(${CERES_INCLUDE_DIRS})

# --- Library Target ---
# Define the library using sources from src/
add_library(polyode_lib STATIC
    src/polynomial.cpp
    src/parameter_estimation.cpp
)
# Specify that the library includes headers from the include/ directory.
# PUBLIC makes these includes available to targets linking against polyode_lib.
target_include_directories(polyode_lib PUBLIC include)
# Add tests directory for test_utils.hpp needed by parameter_estimation.hpp
# PRIVATE because users of the library shouldn't need test utils.
target_include_directories(polyode_lib PUBLIC tests)

# Link Boost if the library itself needs compiled Boost components (unlikely for header-only like Odeint)
# target_link_libraries(polyode_lib PUBLIC Boost::odeint) # Example


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


# --- Testing Framework (Google Test) ---
# Use FetchContent to download and configure Google Test
include(FetchContent)
FetchContent_Declare(
    googletest
    # Specify the repository and tag (version)
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
# This is often needed when integrating gtest via FetchContent
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# Disable warnings-as-errors for Google Test targets
set_target_properties(gtest gmock PROPERTIES COMPILE_WARNING_AS_ERROR FALSE)

# Enable CTest support
enable_testing()

# Include GoogleTest CMake module (provides gtest_discover_tests)
include(GoogleTest)

# --- Test Targets ---
# Helper function to simplify adding test executables and discovering tests
function(add_polyode_test test_name)
    add_executable(${test_name} tests/${test_name}.cpp)
    # Link the test executable against our library and the gtest_main library
    # gtest_main provides a main() function and links gtest and gmock
    target_link_libraries(${test_name} PRIVATE polyode_lib gtest_main)
    # Link math library if needed (based on original CMakeLists)
    if(UNIX AND NOT APPLE)
        target_link_libraries(${test_name} PRIVATE m)
    endif()
    # Use gtest_add_tests instead to explicitly parse source files
    # Correct order: AUTO comes before SOURCES
    gtest_add_tests(TARGET ${test_name} AUTO SOURCES tests/${test_name}.cpp)
endfunction()

# Add each test using the helper function
add_polyode_test(analytic_test)
add_polyode_test(holling_odesystem_test)
# add_polyode_test(rational_ode_test) # Removed as redundant
add_polyode_test(monomial_test)
add_polyode_test(polynomial_test)
add_polyode_test(rational_function_test)
add_polyode_test(variable_test)
add_polyode_test(ode_system_test)
add_polyode_test(exponential_decay_test)


# --- Compiler Warnings ---
# Apply common warning flags to the library target
# INTERFACE ensures these options are used when compiling the library
# and also propagated to targets linking against it (like examples and tests)
if(MSVC)
    target_compile_options(polyode_lib INTERFACE /W4) # Example for MSVC
else()
    target_compile_options(polyode_lib INTERFACE -Wall -Wextra -pedantic) # GCC/Clang
endif()


# --- Installation (Optional - Example structure) ---
# install(TARGETS polyode_lib DESTINATION lib)
# install(DIRECTORY include/ DESTINATION include) # Installs the entire include directory
# install(TARGETS lotka_volterra_example holling_example DESTINATION bin) # Install example executables


# --- Old Targets/CTest commands (Removed for cleanup) ---
