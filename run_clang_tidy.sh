#!/bin/bash
# Simple script to run clang-tidy on all source files

# Ensure the compile database exists
if [ ! -f "build/compile_commands.json" ]; then
  echo "Error: compile_commands.json not found. Run cmake with -DCMAKE_EXPORT_COMPILE_COMMANDS=ON first."
  exit 1
fi

# Find all source files
SOURCE_FILES=$(find src examples tests -name "*.cpp" -type f)

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
  clang-tidy -p build "$file" >> clang-tidy-report.txt 2>&1
  echo "" >> clang-tidy-report.txt
  echo "" >> clang-tidy-report.txt
done

echo "Clang-tidy analysis complete. See clang-tidy-report.txt for details." 