# This file appears to be unused by the main CMake configuration.
# Tests are added via the add_polyode_test() function in the root CMakeLists.txt.
# Commenting out existing content to avoid confusion.

# # Find GTest
# find_package(GTest REQUIRED)
# 
# # Enable testing using CTest
# enable_testing()
# 
# # Add test executables
# # Note: Link PRIVATE GTest::gtest_main to automatically get main()
# 
# add_executable(monomial_test monomial_test.cpp)
# target_link_libraries(monomial_test PRIVATE GTest::gtest_main poly_ode)
# 
# add_executable(polynomial_test polynomial_test.cpp)
# target_link_libraries(polynomial_test PRIVATE GTest::gtest_main poly_ode)
# 
# add_executable(rational_function_test rational_function_test.cpp)
# target_link_libraries(rational_function_test PRIVATE GTest::gtest_main poly_ode)
# 
# add_executable(variable_test variable_test.cpp)
# target_link_libraries(variable_test PRIVATE GTest::gtest_main poly_ode)
# 
# add_executable(ode_system_test ode_system_test.cpp)
# target_link_libraries(ode_system_test PRIVATE GTest::gtest_main poly_ode)
# 
# add_executable(identifiability_analyzer_test identifiability_analyzer_test.cpp)
# target_link_libraries(identifiability_analyzer_test PRIVATE GTest::gtest_main poly_ode)
# 
# add_executable(approximation_test approximation_test.cpp)
# target_link_libraries(approximation_test PRIVATE GTest::gtest_main poly_ode)
# 
# add_executable(exponential_decay_test exponential_decay_test.cpp)
# target_link_libraries(exponential_decay_test PRIVATE GTest::gtest_main poly_ode)
# 
# add_executable(analytic_test analytic_test.cpp)
# target_link_libraries(analytic_test PRIVATE GTest::gtest_main poly_ode)
# 
# add_executable(holling_odesystem_test holling_odesystem_test.cpp)
# target_link_libraries(holling_odesystem_test PRIVATE GTest::gtest_main poly_ode)
# 
# # New test for the parameter estimator setup and symbolic derivatives
# add_executable(parameter_estimator_test parameter_estimator_test.cpp)
# target_link_libraries(parameter_estimator_test PRIVATE GTest::gtest_main poly_ode)
# 
# # --- Old target removed --- 
# # add_executable(parameter_estimation_test OLD_parameter_estimation_test.cpp)
# # target_link_libraries(parameter_estimation_test PRIVATE GTest::gtest_main poly_ode)
# 
# 
# # Discover tests for all defined test executables
# include(GoogleTest)
# gtest_discover_tests(monomial_test)
# gtest_discover_tests(polynomial_test)
# gtest_discover_tests(rational_function_test)
# gtest_discover_tests(variable_test)
# gtest_discover_tests(ode_system_test)
# gtest_discover_tests(identifiability_analyzer_test)
# gtest_discover_tests(approximation_test)
# gtest_discover_tests(exponential_decay_test)
# gtest_discover_tests(analytic_test)
# gtest_discover_tests(holling_odesystem_test)
# gtest_discover_tests(parameter_estimator_test) # Discover tests for the NEW executable
# # gtest_discover_tests(parameter_estimation_test) # Discovery for OLD executable removed 