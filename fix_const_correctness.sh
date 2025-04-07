#!/bin/bash
# Script to run clang-tidy to fix ONLY misc-const-correctness issues

# Determine the build directory and compile_commands.json location
if [ -f "compile_commands.json" ]; then
  # We're in the build directory
  COMPILE_DB="$(pwd)/compile_commands.json"
  BUILD_DIR="$(pwd)"
elif [ -f "build/compile_commands.json" ]; then
  # We're in the project root and compile_commands.json is in build/
  COMPILE_DB="$(pwd)/build/compile_commands.json"
  BUILD_DIR="$(pwd)/build"
elif [ -d "build" ]; then
  # We're in the project root but compile_commands.json isn't generated yet
  echo "Compile commands database not found. Generating it now..."
  cmake -S . -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
  COMPILE_DB="$(pwd)/build/compile_commands.json"
  BUILD_DIR="$(pwd)/build"
else
  echo "Error: Neither compile_commands.json nor build directory found."
  echo "Please run this script from either the build directory or the project root."
  exit 1
fi

echo "Using compile commands database: $COMPILE_DB"

# Find all source files
SOURCE_FILES=$(find src examples tests -name "*.cpp" -type f 2>/dev/null)
if [ -z "$SOURCE_FILES" ]; then
  # If no files found, we might be in the build directory
  SOURCE_FILES=$(find ../src ../examples ../tests -name "*.cpp" -type f 2>/dev/null)
fi

if [ -z "$SOURCE_FILES" ]; then
  echo "Error: No source files found."
  exit 1
fi

echo "Running clang-tidy to fix ONLY misc-const-correctness issues..."
echo "Results will be saved to const-fixes-report.txt"

# Clear the report file
echo "Clang-Tidy Const-Correctness Fixes Report" > const-fixes-report.txt
echo "=======================================" >> const-fixes-report.txt
echo "" >> const-fixes-report.txt

# Define only the misc-const-correctness check
CONST_CHECK="-checks=-*,misc-const-correctness,performance-avoid-endl"

# Process each file and apply fixes
for file in $SOURCE_FILES; do
  echo "Fixing const-correctness in $file..."
  echo "File: $file" >> const-fixes-report.txt
  echo "-------------------" >> const-fixes-report.txt
  
  # Run clang-tidy with -fix to automatically apply const-correctness fixes
  clang-tidy -p "$BUILD_DIR" --header-filter='.*' $CONST_CHECK -fix "$file" >> const-fixes-report.txt 2>&1
  
  echo "" >> const-fixes-report.txt
  echo "" >> const-fixes-report.txt
done

echo "Const-correctness fixes complete. See const-fixes-report.txt for details." 
