#ifndef MODEL_TEST_FRAMEWORK_HPP
#define MODEL_TEST_FRAMEWORK_HPP

#include "test_utils.hpp"
#include <functional>
#include <map>
#include <string>
#include <vector>
#include <optional>

namespace poly_ode {
namespace test_framework {

/**
 * @brief Metadata for a test model including expected properties
 */
struct ModelTestMetadata {
    std::string name;
    std::string category;  // "identifiability", "simple", "classical", "biological", "advanced"
    std::string description;
    
    // Expected identifiability properties
    int expected_identifiable_count = -1;  // -1 means don't check
    int expected_unidentifiable_count = -1;
    std::vector<std::string> expected_identifiable_params;
    std::vector<std::string> expected_unidentifiable_params;
    
    // Expected derivative orders
    std::map<std::string, int> expected_derivative_orders;
    
    // Simulation parameters
    std::optional<std::pair<double, double>> recommended_time_interval;
    double recommended_dt = 0.01;
    double noise_level = 0.0;
    
    // Validation tolerances
    double parameter_tolerance = 1e-4;
    double ic_tolerance = 1e-4;
    double rmse_threshold = 1e-3;
    
    // Special properties
    bool has_rational_functions = false;
    bool is_stiff = false;
    bool requires_special_handling = false;
    std::string special_notes;
};

/**
 * @brief Factory function type for creating test models
 */
using ModelBuilderFunction = std::function<poly_ode::test_utils::OdeSystemTestBuilder()>;

/**
 * @brief Registry for test models organized by category
 */
class ModelTestRegistry {
public:
    struct RegisteredModel {
        ModelBuilderFunction builder_func;
        ModelTestMetadata metadata;
    };
    
    // Register a new test model
    void register_model(const std::string& name, 
                       ModelBuilderFunction builder_func,
                       const ModelTestMetadata& metadata);
    
    // Get model by name
    std::optional<RegisteredModel> get_model(const std::string& name) const;
    
    // Get all models in a category
    std::vector<std::string> get_models_in_category(const std::string& category) const;
    
    // Get all available models
    std::vector<std::string> get_all_model_names() const;
    
    // Get all categories
    std::vector<std::string> get_all_categories() const;
    
    // Get models with specific properties
    std::vector<std::string> get_models_with_unidentifiability() const;
    std::vector<std::string> get_models_with_rational_functions() const;
    
private:
    std::map<std::string, RegisteredModel> models_;
};

/**
 * @brief Standardized test execution framework
 */
class ModelTestExecutor {
public:
    struct TestConfiguration {
        // Identifiability analysis config
        int ident_max_deriv_order;
        int ident_num_test_points;
        double ident_rank_tol;
        double ident_null_tol;
        
        // AAA approximation config
        double aa_abs_tol;
        unsigned int aa_max_order_hint;
        
        // Integration config
        double integration_abs_tol;
        double integration_rel_tol;
        double integration_dt_hint;
        
        // Solution polishing config
        double real_tolerance;
        
        // Data generation config
        int num_time_points;
        bool use_recommended_time_interval;
        
        // Validation config
        bool strict_identifiability_checking;
        bool check_parameter_values;
        bool check_rmse_threshold;
        
        // Constructor with default values
        TestConfiguration() 
            : ident_max_deriv_order(10)
            , ident_num_test_points(5)
            , ident_rank_tol(1e-9)
            , ident_null_tol(1e-6)
            , aa_abs_tol(1e-12)
            , aa_max_order_hint(10)
            , integration_abs_tol(1e-7)
            , integration_rel_tol(1e-7)
            , integration_dt_hint(0.001)
            , real_tolerance(1e-9)
            , num_time_points(21)
            , use_recommended_time_interval(true)
            , strict_identifiability_checking(true)
            , check_parameter_values(true)
            , check_rmse_threshold(true)
        {}
    };
    
    struct TestResults {
        bool success = false;
        std::string error_message;
        
        // Identifiability results
        std::vector<Variable> identified_params;
        std::vector<Variable> unidentified_params;
        std::map<Observable, int> derivative_orders;
        
        // Estimation results
        size_t num_solutions_found = 0;
        size_t num_valid_solutions = 0;
        std::vector<poly_ode::EstimationResult> validated_solutions;
        
        // Best solution metrics
        double best_rmse = std::numeric_limits<double>::infinity();
        std::map<Variable, double> best_parameters;
        std::map<Variable, double> best_initial_conditions;
        
        // Timing information
        double setup_time_ms = 0.0;
        double solve_time_ms = 0.0;
        double validation_time_ms = 0.0;
    };
    
    // Execute a full test for a given model
    TestResults execute_test(const std::string& model_name,
                           PolynomialSolver& solver,
                           const TestConfiguration& config = TestConfiguration());
    
    // Execute test with custom metadata overrides
    TestResults execute_test_with_overrides(const std::string& model_name,
                                          PolynomialSolver& solver,
                                          const ModelTestMetadata& metadata_overrides,
                                          const TestConfiguration& config = TestConfiguration());
    
    // Batch execute tests for multiple models
    std::map<std::string, TestResults> execute_category_tests(const std::string& category,
                                                            PolynomialSolver& solver,
                                                            const TestConfiguration& config = TestConfiguration());
    
private:
    ModelTestRegistry& registry_;
    
public:
    explicit ModelTestExecutor(ModelTestRegistry& registry) : registry_(registry) {}
    
private:
    // Helper methods
    poly_ode::ExperimentalData generate_test_data(const poly_ode::test_utils::OdeSystemTestBuilder& builder,
                                                  const ModelTestMetadata& metadata,
                                                  const TestConfiguration& config);
    
    bool validate_identifiability_results(const poly_ode::EstimationSetupData& setup_data,
                                         const ModelTestMetadata& metadata,
                                         TestResults& results);
    
    bool validate_estimation_results(const std::vector<poly_ode::EstimationResult>& estimation_results,
                                   const poly_ode::test_utils::OdeSystemTestBuilder& builder,
                                   const ModelTestMetadata& metadata,
                                   TestResults& results);
};

/**
 * @brief Utility functions for creating common model metadata patterns
 */
namespace metadata_builders {

// Create metadata for identifiability test models
ModelTestMetadata create_identifiability_metadata(const std::string& name,
                                                 const std::string& description,
                                                 int expected_identifiable_count,
                                                 const std::vector<std::string>& identifiable_params = {},
                                                 const std::vector<std::string>& unidentifiable_params = {});

// Create metadata for classical system models
ModelTestMetadata create_classical_metadata(const std::string& name,
                                           const std::string& description,
                                           const std::pair<double, double>& time_interval,
                                           bool is_stiff = false);

// Create metadata for biological system models
ModelTestMetadata create_biological_metadata(const std::string& name,
                                            const std::string& description,
                                            const std::pair<double, double>& time_interval,
                                            bool has_rational_functions = false);

// Create metadata for simple test models
ModelTestMetadata create_simple_metadata(const std::string& name,
                                        const std::string& description,
                                        int expected_identifiable_count = -1);

} // namespace metadata_builders

// Global registry instance
extern ModelTestRegistry& get_global_model_registry();

} // namespace test_framework
} // namespace poly_ode

#endif // MODEL_TEST_FRAMEWORK_HPP