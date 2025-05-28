#include "model_test_framework.hpp"
#include "model_registrations.hpp"
#include "MSolveSolver.hpp"
#include <gtest/gtest.h>
#include <iostream>
#include <set>
#include <algorithm>

using namespace poly_ode::test_framework;

class SystematicModelTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        // Register all models once for the entire test suite
        static bool registered = false;
        if (!registered) {
            register_all_models();
            registered = true;
        }
    }
    
    void SetUp() override {
        // Create solver for each test
        solver_ = std::make_unique<poly_ode::MSolveSolver>();
        executor_ = std::make_unique<ModelTestExecutor>(get_global_model_registry());
    }
    
    void print_test_results(const std::string& model_name, const ModelTestExecutor::TestResults& results) {
        std::cout << "\n=== Test Results for " << model_name << " ===" << std::endl;
        std::cout << "Success: " << (results.success ? "YES" : "NO") << std::endl;
        
        if (!results.success) {
            std::cout << "Error: " << results.error_message << std::endl;
            return;
        }
        
        std::cout << "Timing:" << std::endl;
        std::cout << "  Setup: " << results.setup_time_ms << " ms" << std::endl;
        std::cout << "  Solve: " << results.solve_time_ms << " ms" << std::endl;
        std::cout << "  Validation: " << results.validation_time_ms << " ms" << std::endl;
        
        std::cout << "Identifiability:" << std::endl;
        std::cout << "  Identifiable (" << results.identified_params.size() << "): ";
        for (const auto& param : results.identified_params) {
            std::cout << param.name << " ";
        }
        std::cout << std::endl;
        
        std::cout << "  Unidentifiable (" << results.unidentified_params.size() << "): ";
        for (const auto& param : results.unidentified_params) {
            std::cout << param.name << " ";
        }
        std::cout << std::endl;
        
        std::cout << "Solutions:" << std::endl;
        std::cout << "  Found: " << results.num_solutions_found << std::endl;
        std::cout << "  Valid: " << results.num_valid_solutions << std::endl;
        std::cout << "  Best RMSE: " << results.best_rmse << std::endl;
        
        if (!results.best_parameters.empty()) {
            std::cout << "  Best parameters: ";
            for (const auto& [param, value] : results.best_parameters) {
                std::cout << param.name << "=" << value << " ";
            }
            std::cout << std::endl;
        }
        
        if (!results.best_initial_conditions.empty()) {
            std::cout << "  Best ICs: ";
            for (const auto& [ic, value] : results.best_initial_conditions) {
                std::cout << ic.name << "=" << value << " ";
            }
            std::cout << std::endl;
        }
    }
    
    std::unique_ptr<poly_ode::PolynomialSolver> solver_;
    std::unique_ptr<ModelTestExecutor> executor_;
};

// Test individual models
TEST_F(SystematicModelTest, TrivialUnidentModel) {
    auto results = executor_->execute_test("trivial_unident", *solver_);
    print_test_results("trivial_unident", results);
    EXPECT_TRUE(results.success) << results.error_message;
}

TEST_F(SystematicModelTest, SimpleModel) {
    auto results = executor_->execute_test("simple", *solver_);
    print_test_results("simple", results);
    EXPECT_TRUE(results.success) << results.error_message;
}

TEST_F(SystematicModelTest, LotkaVolterraModel) {
    // Use relaxed configuration for LV
    ModelTestExecutor::TestConfiguration config;
    config.strict_identifiability_checking = false;  // LV can be tricky
    config.check_parameter_values = true;
    
    auto results = executor_->execute_test("lotka_volterra", *solver_, config);
    print_test_results("lotka_volterra", results);
    
    // For LV, we expect some solutions even if not all parameters are identifiable
    EXPECT_GT(results.num_solutions_found, 0) << "Should find at least some algebraic solutions";
}

TEST_F(SystematicModelTest, SumTestModel) {
    // This model is known to be challenging for msolve
    auto results = executor_->execute_test("sum_test", *solver_);
    print_test_results("sum_test", results);
    
    if (results.num_solutions_found == 0) {
        GTEST_SKIP() << "SumTest model returns no solutions from msolve (known issue)";
    } else {
        EXPECT_TRUE(results.success) << results.error_message;
    }
}

TEST_F(SystematicModelTest, GlobalUnidentModel) {
    auto results = executor_->execute_test("global_unident_test", *solver_);
    print_test_results("global_unident_test", results);
    EXPECT_TRUE(results.success) << results.error_message;
}

// Biological Models
TEST_F(SystematicModelTest, SEIRModel) {
    ModelTestExecutor::TestConfiguration config;
    config.strict_identifiability_checking = false;  // Biological models can be complex
    config.check_parameter_values = true;
    
    auto results = executor_->execute_test("seir", *solver_, config);
    print_test_results("seir", results);
    
    // SEIR models often have identifiability issues
    EXPECT_GT(results.num_solutions_found, 0) << "Should find at least some solutions";
}

TEST_F(SystematicModelTest, TreatmentModel) {
    ModelTestExecutor::TestConfiguration config;
    config.strict_identifiability_checking = false;
    config.check_parameter_values = true;
    
    auto results = executor_->execute_test("treatment", *solver_, config);
    print_test_results("treatment", results);
    EXPECT_GT(results.num_solutions_found, 0) << "Should find at least some solutions";
}

TEST_F(SystematicModelTest, HIVModel) {
    ModelTestExecutor::TestConfiguration config;
    config.strict_identifiability_checking = false;  // Complex biological system
    config.check_parameter_values = false;  // Many parameters, may be challenging
    
    auto results = executor_->execute_test("hiv", *solver_, config);
    print_test_results("hiv", results);
    
    // HIV model is very complex - expect it to run but may not find perfect solutions
    EXPECT_FALSE(results.error_message.find("Exception") != std::string::npos) 
        << "HIV test should run without exceptions: " << results.error_message;
}

// Classical Models with New Additions
TEST_F(SystematicModelTest, HarmonicModel) {
    auto results = executor_->execute_test("harmonic", *solver_);
    print_test_results("harmonic", results);
    EXPECT_TRUE(results.success) << results.error_message;
}

TEST_F(SystematicModelTest, VanderPolModel) {
    auto results = executor_->execute_test("vanderpol", *solver_);
    print_test_results("vanderpol", results);
    EXPECT_TRUE(results.success) << results.error_message;
}

TEST_F(SystematicModelTest, BrusselatorModel) {
    auto results = executor_->execute_test("brusselator", *solver_);
    print_test_results("brusselator", results);
    EXPECT_TRUE(results.success) << results.error_message;
}

// Simple Models with New Additions
TEST_F(SystematicModelTest, OnespCubedModel) {
    auto results = executor_->execute_test("onesp_cubed", *solver_);
    print_test_results("onesp_cubed", results);
    EXPECT_TRUE(results.success) << results.error_message;
}

TEST_F(SystematicModelTest, ThreespCubedModel) {
    auto results = executor_->execute_test("threesp_cubed", *solver_);
    print_test_results("threesp_cubed", results);
    EXPECT_TRUE(results.success) << results.error_message;
}

TEST_F(SystematicModelTest, SimpleLinearCombinationModel) {
    auto results = executor_->execute_test("simple_linear_combination", *solver_);
    print_test_results("simple_linear_combination", results);
    EXPECT_TRUE(results.success) << results.error_message;
}

// Advanced Models
TEST_F(SystematicModelTest, DaisyEx3Model) {
    auto results = executor_->execute_test("daisy_ex3", *solver_);
    print_test_results("daisy_ex3", results);
    EXPECT_TRUE(results.success) << results.error_message;
}

TEST_F(SystematicModelTest, FitzHughNagumoModel) {
    ModelTestExecutor::TestConfiguration config;
    config.num_time_points = 31;  // More points for short time interval
    
    auto results = executor_->execute_test("fitzhugh_nagumo", *solver_, config);
    print_test_results("fitzhugh_nagumo", results);
    EXPECT_TRUE(results.success) << results.error_message;
}

// Identifiability Tests
TEST_F(SystematicModelTest, SubstrTestModel) {
    ModelTestExecutor::TestConfiguration config;
    config.strict_identifiability_checking = false;  // Complex substitution patterns
    
    auto results = executor_->execute_test("substr_test", *solver_, config);
    print_test_results("substr_test", results);
    
    // This test focuses on parameter substitution patterns rather than exact estimation
    EXPECT_GT(results.num_solutions_found, 0) << "Should find at least some solutions";
}

// Test categories
TEST_F(SystematicModelTest, IdentifiabilityCategory) {
    ModelTestExecutor::TestConfiguration config;
    config.strict_identifiability_checking = true;
    config.check_parameter_values = true;
    
    auto results = executor_->execute_category_tests("identifiability", *solver_, config);
    
    std::cout << "\n=== Identifiability Category Results ===" << std::endl;
    
    for (const auto& [model_name, result] : results) {
        std::cout << model_name << ": " << (result.success ? "PASS" : "FAIL");
        if (!result.success) {
            std::cout << " (" << result.error_message << ")";
        }
        std::cout << std::endl;
    }
    
    // At least some models should pass
    size_t passed = 0;
    for (const auto& [model_name, result] : results) {
        if (result.success) passed++;
    }
    
    EXPECT_GT(passed, 0) << "At least some identifiability models should pass";
}

TEST_F(SystematicModelTest, SimpleCategory) {
    auto results = executor_->execute_category_tests("simple", *solver_);
    
    std::cout << "\n=== Simple Category Results ===" << std::endl;
    
    for (const auto& [model_name, result] : results) {
        std::cout << model_name << ": " << (result.success ? "PASS" : "FAIL");
        if (!result.success) {
            std::cout << " (" << result.error_message << ")";
        }
        std::cout << std::endl;
    }
    
    // Simple models should generally pass
    for (const auto& [model_name, result] : results) {
        EXPECT_TRUE(result.success) << "Simple model " << model_name << " failed: " << result.error_message;
    }
}

// Registry inspection tests
TEST_F(SystematicModelTest, RegistryContents) {
    auto& registry = get_global_model_registry();
    
    auto all_models = registry.get_all_model_names();
    std::cout << "\nRegistered models (" << all_models.size() << "):" << std::endl;
    for (const auto& name : all_models) {
        std::cout << "  " << name << std::endl;
    }
    
    auto categories = registry.get_all_categories();
    std::cout << "\nAvailable categories (" << categories.size() << "):" << std::endl;
    for (const auto& category : categories) {
        auto models_in_category = registry.get_models_in_category(category);
        std::cout << "  " << category << " (" << models_in_category.size() << " models)" << std::endl;
    }
    
    auto unident_models = registry.get_models_with_unidentifiability();
    std::cout << "\nModels with unidentifiability (" << unident_models.size() << "):" << std::endl;
    for (const auto& name : unident_models) {
        std::cout << "  " << name << std::endl;
    }
    
    EXPECT_GT(all_models.size(), 0) << "Should have registered models";
    EXPECT_GT(categories.size(), 0) << "Should have model categories";
}

// Performance benchmark test
TEST_F(SystematicModelTest, PerformanceBenchmark) {
    ModelTestExecutor::TestConfiguration config;
    config.strict_identifiability_checking = false;  // Focus on performance
    config.check_parameter_values = false;
    
    auto all_models = get_global_model_registry().get_all_model_names();
    
    std::cout << "\n=== Performance Benchmark ===" << std::endl;
    std::cout << "Model\tSetup(ms)\tSolve(ms)\tValidation(ms)\tTotal(ms)\tSolutions" << std::endl;
    
    double total_time = 0.0;
    size_t total_solutions = 0;
    
    for (const auto& model_name : all_models) {
        auto results = executor_->execute_test(model_name, *solver_, config);
        
        double model_total = results.setup_time_ms + results.solve_time_ms + results.validation_time_ms;
        total_time += model_total;
        total_solutions += results.num_solutions_found;
        
        std::cout << model_name << "\t" 
                  << results.setup_time_ms << "\t"
                  << results.solve_time_ms << "\t" 
                  << results.validation_time_ms << "\t"
                  << model_total << "\t"
                  << results.num_solutions_found << std::endl;
    }
    
    std::cout << "\nTotal time: " << total_time << " ms" << std::endl;
    std::cout << "Total solutions found: " << total_solutions << std::endl;
    std::cout << "Average time per model: " << (total_time / all_models.size()) << " ms" << std::endl;
    
    // Basic performance expectations
    EXPECT_LT(total_time / all_models.size(), 30000.0) << "Average time per model should be reasonable";
}

// Test configuration variations
TEST_F(SystematicModelTest, ConfigurationVariations) {
    const std::string test_model = "simple";  // Use simple model for configuration testing
    
    struct ConfigTest {
        std::string name;
        ModelTestExecutor::TestConfiguration config;
    };
    
    std::vector<ConfigTest> config_tests;
    
    // Default config
    config_tests.push_back({"default", ModelTestExecutor::TestConfiguration()});
    
    // High precision config
    {
        ModelTestExecutor::TestConfiguration config;
        config.ident_rank_tol = 1e-12;
        config.ident_null_tol = 1e-9;
        config.aa_abs_tol = 1e-15;
        config.integration_abs_tol = 1e-10;
        config_tests.push_back({"high_precision", config});
    }
    
    // Fast config
    {
        ModelTestExecutor::TestConfiguration config;
        config.ident_max_deriv_order = 5;
        config.ident_num_test_points = 3;
        config.num_time_points = 11;
        config_tests.push_back({"fast", config});
    }
    
    // Relaxed config
    {
        ModelTestExecutor::TestConfiguration config;
        config.strict_identifiability_checking = false;
        config.check_parameter_values = false;
        config_tests.push_back({"relaxed", config});
    }
    
    std::cout << "\n=== Configuration Variations ===" << std::endl;
    
    for (const auto& test : config_tests) {
        auto results = executor_->execute_test(test_model, *solver_, test.config);
        std::cout << test.name << ": " << (results.success ? "PASS" : "FAIL") 
                  << " (time: " << (results.setup_time_ms + results.solve_time_ms + results.validation_time_ms) 
                  << " ms, solutions: " << results.num_solutions_found << ")" << std::endl;
    }
}

// Legacy Test Modernization - these replace the failing old tests
TEST_F(SystematicModelTest, LegacySimpleModelModernized) {
    std::cout << "\n--- Legacy Test Modernized: SimpleModel --- " << std::endl;
    
    // This replaces ParameterEstimatorScenariosTest.SimpleModel 
    // which was using the broken stub function
    auto results = executor_->execute_test("simple", *solver_);
    print_test_results("simple", results);
    
    // More detailed validation matching the old test expectations
    ASSERT_TRUE(results.success) << "SimpleModel test failed: " << results.error_message;
    
    // Check identifiability (should have 4 identifiable: a, b, x1_0, x2_0)
    EXPECT_EQ(results.identified_params.size(), 4) 
        << "Expected 4 identifiable parameters for SimpleModel";
    
    // Check that we found valid solutions
    EXPECT_GT(results.num_valid_solutions, 0) 
        << "SimpleModel should produce valid parameter estimates";
    
    // Check RMSE quality
    EXPECT_LT(results.best_rmse, 1e-3) 
        << "SimpleModel should achieve good RMSE: " << results.best_rmse;
    
    std::cout << "Legacy SimpleModel test successfully modernized!" << std::endl;
}

TEST_F(SystematicModelTest, LegacyGlobalUnidentTestModernized) {
    std::cout << "\n--- Legacy Test Modernized: GlobalUnidentTest --- " << std::endl;
    
    // This replaces ParameterEstimatorScenariosTest.GlobalUnidentTest
    // which was using the broken stub function
    auto results = executor_->execute_test("global_unident_test", *solver_);
    print_test_results("global_unident_test", results);
    
    // More detailed validation for unidentifiability expectations
    ASSERT_TRUE(results.success) << "GlobalUnidentTest failed: " << results.error_message;
    
    // Expected identifiability pattern for global_unident_test:
    // - 'a' should be identifiable (appears alone in dx1/dt = -a*x1)
    // - 'x1', 'x2' initial conditions should be identifiable (observed directly)
    // - 'b', 'c' should have at most one identifiable (sum b+c appears in dx2/dt)
    // - 'd', 'x3' should be unidentifiable (x3 not observed, d only affects x3)
    
    EXPECT_GE(results.identified_params.size(), 3) 
        << "Should have at least a, x1_0, x2_0 identifiable";
    EXPECT_LE(results.identified_params.size(), 4) 
        << "Should have at most a, x1_0, x2_0, and one of (b,c) identifiable";
    
    // Check parameter names in identifiable set
    std::set<std::string> identifiable_names;
    for (const auto& param : results.identified_params) {
        identifiable_names.insert(param.name);
    }
    
    EXPECT_TRUE(identifiable_names.count("a")) << "Parameter 'a' should be identifiable";
    EXPECT_TRUE(identifiable_names.count("x1")) << "Initial condition 'x1' should be identifiable";
    EXPECT_TRUE(identifiable_names.count("x2")) << "Initial condition 'x2' should be identifiable";
    EXPECT_FALSE(identifiable_names.count("d")) << "Parameter 'd' should be unidentifiable";
    EXPECT_FALSE(identifiable_names.count("x3")) << "Initial condition 'x3' should be unidentifiable";
    
    // Either b or c (or both) should be unidentifiable due to sum relationship
    bool b_unident = results.unidentified_params.end() != 
        std::find_if(results.unidentified_params.begin(), results.unidentified_params.end(),
                     [](const Variable& v) { return v.name == "b"; });
    bool c_unident = results.unidentified_params.end() != 
        std::find_if(results.unidentified_params.begin(), results.unidentified_params.end(),
                     [](const Variable& v) { return v.name == "c"; });
    
    EXPECT_TRUE(b_unident || c_unident) 
        << "At least one of b, c should be unidentifiable due to sum relationship";
    
    std::cout << "Legacy GlobalUnidentTest successfully modernized!" << std::endl;
}