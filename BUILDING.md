# Building PolyODE

This document provides detailed build instructions for developers and contributors.

## Overview

PolyODE uses:
- **CMake** for build system configuration
- **vcpkg** for dependency management with custom overlay ports
- **C++17** standard

## Dependencies

### Required Libraries (managed by vcpkg)

- **FLINT 3.2.1** - Fast Library for Number Theory (custom overlay port)
- **msolve 0.7.5** - Multivariate polynomial system solver (custom overlay port)
- **Ceres Solver** - Non-linear least squares optimization
- **Boost** (components: odeint, filesystem, compute, headers)
- **Google Test** - Testing framework
- **Eigen3** - Linear algebra (dependency of Ceres)
- **GMP/MPFR** - Multiple precision arithmetic (dependencies of FLINT)

### System Requirements

- C++17 compatible compiler (GCC 7+, Clang 6+, MSVC 2019+)
- CMake 3.15 or higher
- Git
- vcpkg package manager

## Quick Build

```bash
# Install vcpkg if not already available
git clone https://github.com/microsoft/vcpkg.git ~/vcpkg
cd ~/vcpkg && ./bootstrap-vcpkg.sh

# Build the project
cd /path/to/polyODE
./build.sh
```

## Manual Build Process

### 1. vcpkg Setup

Install vcpkg if you haven't already:

```bash
git clone https://github.com/microsoft/vcpkg.git ~/vcpkg
cd ~/vcpkg
./bootstrap-vcpkg.sh
```

### 2. Configure Build

```bash
mkdir build && cd build
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_TOOLCHAIN_FILE=~/vcpkg/scripts/buildsystems/vcpkg.cmake \
    -DENABLE_SANITIZERS=OFF
```

### 3. Build

```bash
cmake --build . --parallel
```

## Build Options

### CMake Options

- `CMAKE_BUILD_TYPE`: `Debug`, `Release`, `RelWithDebInfo`, `MinSizeRel`
- `ENABLE_SANITIZERS`: `ON`/`OFF` - Enable AddressSanitizer and UBSanitizer
- `CMAKE_TOOLCHAIN_FILE`: Path to vcpkg toolchain file

### Build Script Options

The `build.sh` script supports several options:

```bash
./build.sh --help                    # Show help
./build.sh                          # Default Release build
./build.sh --debug                  # Debug build
./build.sh --sanitizers             # Enable sanitizers
./build.sh --clean                  # Clean build
./build.sh --verbose                # Verbose output
./build.sh --debug --sanitizers     # Debug with sanitizers
```

## Custom Overlay Ports

This project includes custom vcpkg overlay ports that are essential for building:

### FLINT 3.2.1 (`vcpkg-overlay-ports/flint/`)

**Purpose**: Provides FLINT 3.2.1 with fixes for vcpkg compatibility.

**Key Fixes**:
- Copies generated headers (`flint.h`, `flint-config.h`) from build tree
- Fixes pkg-config file paths (`includedir=${prefix}/include`)
- Proper header directory structure

**Files**:
- `vcpkg.json` - Package metadata
- `portfile.cmake` - Build instructions
- `usage` - Usage instructions

### msolve 0.7.5 (`ports/msolve/`)

**Purpose**: Multivariate polynomial system solver with FLINT 3.x compatibility.

**Key Fixes**:
- Sets proper environment variables for FLINT discovery
- Configures autotools build system for vcpkg
- Links against FLINT, GMP, MPFR

**Files**:
- `vcpkg.json` - Package metadata  
- `portfile.cmake` - Build instructions

## Troubleshooting

### Common Build Issues

#### 1. FLINT Headers Not Found
```
fatal error: flint/flint.h: No such file or directory
```

**Solution**: Ensure you're using the custom FLINT overlay port included in this repo. The standard vcpkg FLINT port may not work.

#### 2. msolve Build Fails
```
error: building msolve:x64-linux failed with: BUILD_FAILED
```

**Solution**: 
- Make sure FLINT built successfully first
- Check that overlay ports are being used (`vcpkg-configuration.json`)
- Clean vcpkg cache: `vcpkg remove --outdated`

#### 3. vcpkg Manifest Errors
```
Error: Could not locate a manifest file
```

**Solution**: 
- Make sure you're in the project root directory
- Verify `vcpkg.json` exists
- Update vcpkg: `cd ~/vcpkg && git pull && ./bootstrap-vcpkg.sh`

#### 4. Pkg-config Path Issues
```
Imported target "PkgConfig::FLINT" includes non-existent path "/include"
```

**Solution**: This should be fixed by the custom FLINT overlay port. If it persists, manually edit the pkg-config file or rebuild FLINT.

### Debug Build Issues

Debug builds may take significantly longer due to:
- Unoptimized code generation
- Additional debugging symbols
- Sanitizer overhead (if enabled)

### Memory Usage

Building FLINT and msolve can be memory-intensive. If you encounter out-of-memory errors:

- Reduce parallel build jobs: `cmake --build . --parallel 2`
- Close other applications
- Consider building in Release mode first

## Development Workflow

### Adding New Dependencies

1. **Standard vcpkg package**:
   ```bash
   # Add to vcpkg.json
   "dependencies": ["new-package"]
   
   # Update CMakeLists.txt
   find_package(NewPackage REQUIRED)
   target_link_libraries(polyode_lib PUBLIC NewPackage::NewPackage)
   ```

2. **Custom overlay port**:
   - Create port in `ports/` or `vcpkg-overlay-ports/`
   - Add to `vcpkg-configuration.json` overlay-ports list
   - Follow vcpkg port guidelines

### Testing Changes

```bash
# Build and test
./build.sh --clean --debug
cd build && ctest --verbose

# Test specific component
./build/polynomial_test
./build/parameter_estimation_test
```

### Code Quality

```bash
# Enable sanitizers for development
./build.sh --debug --sanitizers

# Generate compile commands for IDE
# (automatically done by build script)
```

## Platform-Specific Notes

### Linux
- Requires build-essential or equivalent
- May need additional packages: `libgmp-dev libmpfr-dev`

### macOS  
- Requires Xcode Command Line Tools
- Homebrew packages: `cmake git`

### Windows
- Requires Visual Studio 2019+ with C++ workload
- Use `bootstrap-vcpkg.bat` instead of `.sh`
- Build script not supported (use manual process)

## Performance Notes

### Build Time Optimization
- Use parallel builds: `--parallel` (default in build script)
- Pre-built vcpkg binary cache speeds up subsequent builds
- Release builds are faster to compile than Debug

### Runtime Performance
- Release builds are significantly faster than Debug
- Sanitizers add substantial runtime overhead
- Consider `RelWithDebInfo` for debugging optimized code

## Continuous Integration

For CI systems, consider:

```yaml
# GitHub Actions example
- name: Setup vcpkg
  run: |
    git clone https://github.com/microsoft/vcpkg.git
    ./vcpkg/bootstrap-vcpkg.sh
    
- name: Build
  run: |
    cmake -B build -DCMAKE_TOOLCHAIN_FILE=vcpkg/scripts/buildsystems/vcpkg.cmake
    cmake --build build --parallel
    
- name: Test  
  run: cd build && ctest --output-on-failure
```