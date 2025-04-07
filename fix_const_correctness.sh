#!/bin/bash
# Script to run clang-tidy to fix ONLY misc-const-correctness issues

# Ensure the compile database exists
if [ ! -f "build/compile_commands.json" ]; then
  echo "Error: compile_commands.json not found. Run cmake with -DCMAKE_EXPORT_COMPILE_COMMANDS=ON first."
  exit 1
fi

# Find all source files
SOURCE_FILES=$(find src examples tests -name "*.cpp" -type f)

echo "Running clang-tidy to fix ONLY misc-const-correctness issues..."
echo "Results will be saved to const-fixes-report.txt"

# Clear the report file
echo "Clang-Tidy Const-Correctness Fixes Report" > const-fixes-report.txt
echo "=======================================" >> const-fixes-report.txt
echo "" >> const-fixes-report.txt

# Define only the misc-const-correctness check
CONST_CHECK="-checks=-*,misc-const-correctness"

# Process each file and apply fixes
for file in $SOURCE_FILES; do
  echo "Fixing const-correctness in $file..."
  echo "File: $file" >> const-fixes-report.txt
  echo "-------------------" >> const-fixes-report.txt
  
  # Run clang-tidy with -fix to automatically apply const-correctness fixes
  clang-tidy -p build $CONST_CHECK -fix "$file" >> const-fixes-report.txt 2>&1
  
  echo "" >> const-fixes-report.txt
  echo "" >> const-fixes-report.txt
done

echo "Const-correctness fixes complete. See const-fixes-report.txt for details." 