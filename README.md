# PolyODE: Polynomial ODE System Toolkit

A C++ library for analyzing polynomial and rational function systems of ordinary differential equations (ODEs), with support for parameter estimation, identifiability analysis, and algebraic solving.

## Quick Start

```bash
# Clone the repository
git clone <your-repo-url>
cd polyODE

# Build (requires vcpkg - see setup below if you don't have it)
./build.sh

# Run examples
./build/lotka_volterra
./build/parameter_estimation_example
```

## Prerequisites

1. **C++17 compiler** (GCC 7+, Clang 6+, or MSVC 2019+)
2. **CMake 3.15+**
3. **Git**
4. **vcpkg** (for dependency management)

## Dependencies

This project depends on several libraries that are automatically managed via vcpkg:
- **FLINT 3.2.1** - Fast Library for Number Theory (custom overlay)
- **msolve 0.7.5** - Multivariate polynomial system solver (custom overlay) 
- **Ceres Solver** - Non-linear optimization
- **Boost** - Various utilities (odeint, filesystem, etc.)
- **Google Test** - Testing framework

## Building

### Option 1: Using the Build Script (Recommended)

```bash
# Make sure vcpkg is installed (see vcpkg setup below)
./build.sh
```

The build script assumes vcpkg is installed at `~/vcpkg`. If it's elsewhere, edit `build.sh` to point to your vcpkg installation.

### Option 2: Manual CMake Build

```bash
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release \
         -DCMAKE_TOOLCHAIN_FILE=/path/to/vcpkg/scripts/buildsystems/vcpkg.cmake
cmake --build .
```

## vcpkg Setup

If you don't have vcpkg installed:

```bash
# Install vcpkg (recommended location: ~/vcpkg)
git clone https://github.com/microsoft/vcpkg.git ~/vcpkg
cd ~/vcpkg
./bootstrap-vcpkg.sh  # or .\bootstrap-vcpkg.bat on Windows

# Optional: integrate with CMake globally
./vcpkg integrate install
```

## Custom Overlay Ports

This project includes custom vcpkg overlay ports for libraries not in the main vcpkg registry:

- `vcpkg-overlay-ports/flint/` - FLINT 3.2.1 with fixes for proper header installation
- `ports/msolve/` - msolve 0.7.5 with FLINT 3.x compatibility

These are automatically used when building via vcpkg.

## Project Structure

```
polyODE/
├── include/           # Header files
├── src/              # Source files  
├── examples/         # Example programs
├── tests/            # Unit tests
├── ports/            # Custom vcpkg ports (msolve)
├── vcpkg-overlay-ports/  # More custom vcpkg ports (FLINT)
├── scripts/          # Utility scripts
├── build.sh          # Build script
├── vcpkg.json        # vcpkg manifest
└── CMakeLists.txt    # Main CMake configuration
```

## Testing

```bash
# Build and run all tests
cd build
ctest

# Run specific test
./polynomial_test
./parameter_estimation_test

# Run with verbose output
ctest --verbose
```

## Examples

The `examples/` directory contains sample programs demonstrating key functionality:

- `basic_estimation.cpp` - Basic parameter estimation
- `lotka_volterra.cpp` - Classic predator-prey model
- `identifiability_test.cpp` - Parameter identifiability analysis
- `holling_test.cpp` - Holling-type functional response

## Troubleshooting

### Common Issues

1. **FLINT headers not found**: Make sure you're using the custom overlay ports included in this repo
2. **msolve build fails**: Ensure FLINT built successfully first
3. **vcpkg manifest errors**: Make sure vcpkg is up to date (`git pull` in vcpkg directory)

### Build Script Fails

If `./build.sh` fails, try:

```bash
# Clean and rebuild
rm -rf build
./build.sh

# Or build manually with more verbose output
mkdir build && cd build
cmake .. -DCMAKE_TOOLCHAIN_FILE=~/vcpkg/scripts/buildsystems/vcpkg.cmake -DCMAKE_VERBOSE_MAKEFILE=ON
make VERBOSE=1
```

### vcpkg Issues

```bash
# Update vcpkg
cd ~/vcpkg
git pull
./bootstrap-vcpkg.sh

# Clear vcpkg cache if needed
./vcpkg remove --outdated
```

## Contributing

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass: `ctest`
5. Submit a pull request

## License

[Your license here]

## Citations

If you use this software in academic work, please cite:

```
[Your citation format]
```

## Development Notes

### Adding New Dependencies

If you need to add new dependencies:

1. Add to `vcpkg.json` if available in vcpkg registry
2. Create custom overlay port in `ports/` or `vcpkg-overlay-ports/` if not available
3. Update `CMakeLists.txt` to find and link the package
4. Update this README's dependency list

### Custom Overlay Ports

The FLINT and msolve overlay ports include specific fixes:
- FLINT: Copies generated headers (`flint.h`, `flint-config.h`) and fixes pkg-config paths
- msolve: Sets proper environment variables for FLINT discovery

These fixes are necessary for FLINT 3.x compatibility and should be preserved when updating.