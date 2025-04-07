#!/bin/bash
# Simple script to run clang-tidy on all source files

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

echo "Running clang-tidy on all source files..."
echo "Results will be saved to clang-tidy-report.txt"

# Clear the report file
echo "Clang-Tidy Report" > clang-tidy-report.txt
echo "==================" >> clang-tidy-report.txt
echo "" >> clang-tidy-report.txt

# Process each file
for file in $SOURCE_FILES; do
  echo "Checking $file..."
  echo "File: $file" >> clang-tidy-report.txt
  echo "-------------------" >> clang-tidy-report.txt
  clang-tidy -p "$BUILD_DIR" --header-filter='.*' "$file" >> clang-tidy-report.txt 2>&1
  echo "" >> clang-tidy-report.txt
  echo "" >> clang-tidy-report.txt
done

echo "Clang-tidy analysis complete. See clang-tidy-report.txt for details." 