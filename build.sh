#!/bin/bash

# PolyODE Build Script
# Automatically builds the project with vcpkg dependencies

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}PolyODE Build Script${NC}"
echo "===================="

# Check for vcpkg installation
VCPKG_ROOT=""
if [ -n "$VCPKG_ROOT" ]; then
    echo -e "${GREEN}Using vcpkg from VCPKG_ROOT: $VCPKG_ROOT${NC}"
elif [ -d "$HOME/vcpkg" ]; then
    VCPKG_ROOT="$HOME/vcpkg"
    echo -e "${GREEN}Found vcpkg at: $VCPKG_ROOT${NC}"
elif [ -d "/usr/local/share/vcpkg" ]; then
    VCPKG_ROOT="/usr/local/share/vcpkg"
    echo -e "${GREEN}Found vcpkg at: $VCPKG_ROOT${NC}"
else
    echo -e "${RED}ERROR: vcpkg not found!${NC}"
    echo ""
    echo "Please install vcpkg first:"
    echo "  git clone https://github.com/microsoft/vcpkg.git ~/vcpkg"
    echo "  cd ~/vcpkg"
    echo "  ./bootstrap-vcpkg.sh"
    echo ""
    echo "Or set VCPKG_ROOT environment variable to your vcpkg installation."
    exit 1
fi

# Verify vcpkg toolchain file exists
TOOLCHAIN_FILE="$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake"
if [ ! -f "$TOOLCHAIN_FILE" ]; then
    echo -e "${RED}ERROR: vcpkg toolchain file not found at: $TOOLCHAIN_FILE${NC}"
    echo "Make sure vcpkg is properly bootstrapped."
    exit 1
fi

# Parse command line arguments
BUILD_TYPE="Release"
ENABLE_SANITIZERS="OFF"
CLEAN_BUILD=false
VERBOSE=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --debug)
            BUILD_TYPE="Debug"
            shift
            ;;
        --sanitizers)
            ENABLE_SANITIZERS="ON"
            shift
            ;;
        --clean)
            CLEAN_BUILD=true
            shift
            ;;
        --verbose)
            VERBOSE=true
            shift
            ;;
        --help)
            echo "Usage: $0 [options]"
            echo "Options:"
            echo "  --debug        Build in Debug mode (default: Release)"
            echo "  --sanitizers   Enable AddressSanitizer/UBSan (default: OFF)"
            echo "  --clean        Clean build directory first"
            echo "  --verbose      Verbose build output"
            echo "  --help         Show this help message"
            exit 0
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            echo "Use --help for usage information."
            exit 1
            ;;
    esac
done

echo ""
echo "Build Configuration:"
echo "  Build Type: $BUILD_TYPE"
echo "  Sanitizers: $ENABLE_SANITIZERS"
echo "  vcpkg Root: $VCPKG_ROOT"
echo ""

# Clean build directory if requested or if it doesn't exist
if [ "$CLEAN_BUILD" = true ] || [ ! -d "build" ]; then
    echo -e "${YELLOW}Cleaning build directory...${NC}"
    rm -rf build
fi

# Create build directory
mkdir -p build
cd build

# Configure with CMake
echo -e "${BLUE}Configuring with CMake...${NC}"
cmake .. \
    -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
    -DENABLE_SANITIZERS="$ENABLE_SANITIZERS" \
    -DCMAKE_TOOLCHAIN_FILE="$TOOLCHAIN_FILE" \
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON

# Build
echo -e "${BLUE}Building...${NC}"
if [ "$VERBOSE" = true ]; then
    cmake --build . --verbose
else
    cmake --build . --parallel
fi

echo ""
echo -e "${GREEN}✓ Build completed successfully!${NC}"
echo ""
echo "Built executables are in the build/ directory."
echo "Run tests with: cd build && ctest"
echo "Run examples: ./build/<example_name>"