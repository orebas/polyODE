# PolyODE: Polynomial ODE System Toolkit

This repository contains C++ tools for defining, simulating, and analyzing systems of ordinary differential equations (ODEs) described by polynomials or rational functions.

## Prerequisites

Before building, ensure you have the following installed:

1.  **CMake**: Version 3.15 or higher. ([Download CMake](https://cmake.org/download/))
2.  **A C++17 Compiler**:
    *   **Windows**: Visual Studio 2019 or later (with C++ workload).
    *   **macOS**: Xcode Command Line Tools (install with `xcode-select --install`).
    *   **Linux**: GCC 7 or later, or Clang 6 or later. Install via your package manager (e.g., `sudo apt update && sudo apt install build-essential g++` on Debian/Ubuntu).
3.  **Git**: Required for cloning this repository and for vcpkg. ([Download Git](https://git-scm.com/downloads))

## Building with vcpkg (Recommended)

This project uses [vcpkg](https://github.com/microsoft/vcpkg) for dependency management (Boost, Ceres Solver, Google Test).

1.  **Clone the Repository**:
    ```bash
    git clone <your-repository-url> # Replace with your repo URL
    cd polyODE # Or your repository's directory name
    ```

2.  **Install vcpkg**:
    *   Clone the vcpkg repository (it's recommended to install it outside your project directory, e.g., in your home directory or `C:\tools`).
      ```bash
      git clone https://github.com/microsoft/vcpkg.git
      cd vcpkg
      ```
    *   Bootstrap vcpkg:
      *   Windows: `.\bootstrap-vcpkg.bat`
      *   Linux/macOS: `./bootstrap-vcpkg.sh`
    *   *(Optional but Recommended)* Integrate with user-wide CMake (makes finding vcpkg easier):
      *   Windows: `.\vcpkg integrate install`
      *   Linux/macOS: `./vcpkg integrate install`

3.  **Configure with CMake**:
    *   Create a build directory:
      ```bash
      cd ../polyODE # Go back to your project directory
      mkdir build
      cd build
      ```
    *   Run CMake, telling it where to find vcpkg's toolchain file. **Replace `<path/to/vcpkg>` with the actual path where you cloned vcpkg.**
      ```bash
      cmake .. -DCMAKE_TOOLCHAIN_FILE=<path/to/vcpkg>/scripts/buildsystems/vcpkg.cmake
      ```
      *   *Note*: If you ran `vcpkg integrate install`, CMake *might* find vcpkg automatically, and you *might* be able to omit the `-DCMAKE_TOOLCHAIN_FILE=...` part, but explicitly specifying it is more reliable.
      *   *Optional CMake Options*:
          *   `-DCMAKE_BUILD_TYPE=Release` (or `Debug`)
          *   `-DENABLE_SANITIZERS=OFF` (to disable ASan/UBSan, ON by default)

4.  **Build**:
    *   Use CMake's build command (works cross-platform):
      ```bash
      cmake --build .
      ```
    *   Alternatively, use the native build tools (e.g., `make` on Linux/macOS, or open the generated `.sln` file in Visual Studio on Windows).

5.  **Run Examples/Tests**:
    *   Executables will be located in the `build` directory (or subdirectories like `build/Debug` depending on configuration).
    *   Run an example: `./lotka_volterra_example` (or `.\Debug\lotka_volterra_example.exe` on Windows).
    *   Run tests using CTest: `ctest`

This setup ensures that all required libraries (Boost, Ceres, GTest) are automatically downloaded, built, and found by CMake during the configuration step.

## Testing

The project uses Google Test for unit and integration testing.

### Running Tests

1.  Configure the build using CMake (from a `build` directory):
    ```bash
    cd build
    cmake ..
    ```
2.  Build the tests:
    ```bash
    # Using Ninja (if configured)
    ninja
    # Or using Make
    make -j
    ```
3.  Run all tests using CTest:
    ```bash
    # From the build directory
    ctest
    # Or via the build tool
    ninja test 
    # make test
    ```
4.  Run a specific test executable:
    ```bash
    # From the build directory
    ./<test_executable_name> --gtest_filter=TestSuiteName.TestCaseName
    # e.g., ./polynomial_test --gtest_filter=PolynomialTest.Arithmetic
    ```

### Adding New Tests

1.  Create a new test source file in the `tests/` directory (e.g., `tests/my_feature_test.cpp`). Ensure the filename follows the pattern `<target_name>_test.cpp`.
2.  Include necessary headers and `<gtest/gtest.h>`.
3.  Write test cases using Google Test macros (e.g., `TEST`, `TEST_F`, `EXPECT_EQ`, `ASSERT_TRUE`).
4.  Add **one line** to the main `CMakeLists.txt` file in the "Test Targets" section:
    ```cmake
    # Add each test using the helper function
    # ... existing add_polyode_test calls ...
    add_polyode_test(my_feature_test) # Add this line (use the target name without _test.cpp)
    ```
5.  **If your test needs extra libraries** not already linked by `polyode_lib` (like Ceres, nlohmann_json, etc.):
    *   Ensure the dependency is found (e.g., `find_package(MyLib REQUIRED)`).
    *   Add a line *after* the `add_polyode_test` call to link it:
        ```cmake
        add_polyode_test(my_feature_test)
        # Example: Link nlohmann_json if needed by my_feature_test.cpp
        target_link_libraries(my_feature_test PRIVATE nlohmann_json::nlohmann_json)
        ```
6.  Re-run CMake configuration and build. 