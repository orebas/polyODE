# run_clang_tidy.cmake
# This script runs clang-tidy on a list of source files

# Convert the list of source files to a CMake list
string(REPLACE " " ";" SOURCE_FILES "${ALL_CXX_SOURCES}")

# Create or clear the report file
file(WRITE "${SOURCE_DIR}/clang-tidy-report.txt" "Clang-Tidy Report\n==================\n\n")

# Process each source file
foreach(SOURCE_FILE ${SOURCE_FILES})
    # Construct the full path to the source file
    set(FULL_PATH "${SOURCE_DIR}/${SOURCE_FILE}")
    
    # Print a status message
    message(STATUS "Running clang-tidy on ${SOURCE_FILE}")
    
    # Run clang-tidy and capture the output
    execute_process(
        COMMAND ${CLANG_TIDY_EXE} -p ${CMAKE_BINARY_DIR} ${FULL_PATH}
        OUTPUT_VARIABLE TIDY_OUTPUT
        ERROR_VARIABLE TIDY_ERROR
        RESULT_VARIABLE TIDY_RESULT
    )
    
    # Append the result to the report file
    file(APPEND "${SOURCE_DIR}/clang-tidy-report.txt" "File: ${SOURCE_FILE}\n")
    file(APPEND "${SOURCE_DIR}/clang-tidy-report.txt" "Result: ${TIDY_RESULT}\n")
    
    if(TIDY_OUTPUT)
        file(APPEND "${SOURCE_DIR}/clang-tidy-report.txt" "Output:\n${TIDY_OUTPUT}\n\n")
    endif()
    
    if(TIDY_ERROR)
        file(APPEND "${SOURCE_DIR}/clang-tidy-report.txt" "Errors/Warnings:\n${TIDY_ERROR}\n\n")
    endif()
    
    file(APPEND "${SOURCE_DIR}/clang-tidy-report.txt" "-------------------------------------------\n\n")
endforeach()

message(STATUS "Clang-tidy analysis complete. See ${SOURCE_DIR}/clang-tidy-report.txt for details.") 