#include "model_test_framework.hpp"
#include "parameter_estimator.hpp"
#include "approximation/aa_approximator.hpp"
#include <chrono>
#include <algorithm>
#include <set>

namespace poly_ode {
namespace test_framework {

// ModelTestRegistry implementation
void ModelTestRegistry::register_model(const std::string& name, 
                                      ModelBuilderFunction builder_func,
                                      const ModelTestMetadata& metadata) {
    models_[name] = RegisteredModel{std::move(builder_func), metadata};
}

std::optional<ModelTestRegistry::RegisteredModel> 
ModelTestRegistry::get_model(const std::string& name) const {
    auto it = models_.find(name);
    if (it != models_.end()) {
        return it->second;
    }
    return std::nullopt;
}

std::vector<std::string> 
ModelTestRegistry::get_models_in_category(const std::string& category) const {
    std::vector<std::string> result;
    for (const auto& [name, model] : models_) {
        if (model.metadata.category == category) {
            result.push_back(name);
        }
    }
    return result;
}

std::vector<std::string> ModelTestRegistry::get_all_model_names() const {
    std::vector<std::string> result;
    for (const auto& [name, model] : models_) {
        result.push_back(name);
    }
    return result;
}

std::vector<std::string> ModelTestRegistry::get_all_categories() const {
    std::set<std::string> categories;
    for (const auto& [name, model] : models_) {
        categories.insert(model.metadata.category);
    }
    return std::vector<std::string>(categories.begin(), categories.end());
}

std::vector<std::string> ModelTestRegistry::get_models_with_unidentifiability() const {
    std::vector<std::string> result;
    for (const auto& [name, model] : models_) {
        if (model.metadata.expected_unidentifiable_count > 0 || 
            !model.metadata.expected_unidentifiable_params.empty()) {
            result.push_back(name);
        }
    }
    return result;
}

std::vector<std::string> ModelTestRegistry::get_models_with_rational_functions() const {
    std::vector<std::string> result;
    for (const auto& [name, model] : models_) {
        if (model.metadata.has_rational_functions) {
            result.push_back(name);
        }
    }
    return result;
}

// ModelTestExecutor implementation
ModelTestExecutor::TestResults 
ModelTestExecutor::execute_test(const std::string& model_name,
                               PolynomialSolver& solver,
                               const TestConfiguration& config) {
    auto model_opt = registry_.get_model(model_name);
    if (!model_opt) {
        TestResults results;
        results.error_message = "Model '" + model_name + "' not found in registry";
        return results;
    }
    
    return execute_test_with_overrides(model_name, solver, model_opt->metadata, config);
}

ModelTestExecutor::TestResults 
ModelTestExecutor::execute_test_with_overrides(const std::string& model_name,
                                             PolynomialSolver& solver,
                                             const ModelTestMetadata& metadata,
                                             const TestConfiguration& config) {
    TestResults results;
    auto start_time = std::chrono::high_resolution_clock::now();
    
    try {
        // Get the model
        auto model_opt = registry_.get_model(model_name);
        if (!model_opt) {
            results.error_message = "Model '" + model_name + "' not found in registry";
            return results;
        }
        
        // Build the system
        auto builder = model_opt->builder_func();
        auto system = builder.get_system();
        
        // Generate test data
        auto data = generate_test_data(builder, metadata, config);
        
        // Determine parameters to analyze (all parameters and initial conditions)
        std::vector<Variable> params_to_analyze;
        for (const auto& param : system.parameters) {
            params_to_analyze.push_back(param);
        }
        for (const auto& state : system.state_variables) {
            params_to_analyze.push_back(state);  // Initial conditions
        }
        
        auto setup_end = std::chrono::high_resolution_clock::now();
        results.setup_time_ms = std::chrono::duration<double, std::milli>(setup_end - start_time).count();
        
        // Run identifiability analysis and setup
        auto solve_start = std::chrono::high_resolution_clock::now();
        
        EstimationSetupData setup_data = setup_estimation(
            system, 
            params_to_analyze,
            config.ident_max_deriv_order,
            config.ident_num_test_points,
            config.ident_rank_tol,
            config.ident_null_tol
        );
        
        // Validate identifiability results
        if (config.strict_identifiability_checking) {
            if (!validate_identifiability_results(setup_data, metadata, results)) {
                return results;
            }
        }
        
        // Store identifiability results
        results.identified_params = setup_data.identifiable_parameters;
        results.unidentified_params.clear();
        for (const auto& [param, value] : setup_data.non_identifiable_parameters) {
            results.unidentified_params.push_back(param);
        }
        results.derivative_orders = setup_data.required_derivative_orders;
        
        // Fit approximators and evaluate derivatives
        double t_eval = data.times[data.times.size() / 2];  // Middle time point
        std::map<Variable, double> approx_obs_values;
        
        for (const auto& [obs, order] : setup_data.required_derivative_orders) {
            AAApproximator<double> approximator(config.aa_abs_tol, 100, order + 1);
            approximator.fit(data.times, data.measurements.at(obs));
            
            for (int deriv_order = 0; deriv_order <= order; ++deriv_order) {
                Variable obs_deriv_var(obs.name, deriv_order);
                approx_obs_values[obs_deriv_var] = approximator.derivative(t_eval, deriv_order);
            }
        }
        
        // Create estimator and solve
        ParameterEstimator estimator(solver, setup_data, approx_obs_values, t_eval);
        
        auto solve_start_algebraic = std::chrono::high_resolution_clock::now();
        PolynomialSolutionSet solutions = estimator.solve();
        auto solve_end = std::chrono::high_resolution_clock::now();
        
        results.solve_time_ms = std::chrono::duration<double, std::milli>(solve_end - solve_start_algebraic).count();
        results.num_solutions_found = solutions.size();
        
        // Validate solutions
        auto validation_start = std::chrono::high_resolution_clock::now();
        
        if (!solutions.empty()) {
            double rmse_threshold = (metadata.rmse_threshold > 0) ? metadata.rmse_threshold : config.integration_abs_tol * 10;
            
            std::vector<EstimationResult> estimation_results = estimator.process_solutions_and_validate(
                solutions,
                system,
                data,
                data.times.front(),
                rmse_threshold,
                config.integration_abs_tol,
                config.integration_rel_tol,
                config.integration_dt_hint,
                config.real_tolerance
            );
            
            results.validated_solutions = estimation_results;
            results.num_valid_solutions = estimation_results.size();
            
            if (config.check_parameter_values) {
                validate_estimation_results(estimation_results, builder, metadata, results);
            }
        }
        
        auto validation_end = std::chrono::high_resolution_clock::now();
        results.validation_time_ms = std::chrono::duration<double, std::milli>(validation_end - validation_start).count();
        
        // Find best solution
        if (!results.validated_solutions.empty()) {
            auto best_it = std::min_element(results.validated_solutions.begin(), 
                                          results.validated_solutions.end(),
                                          [](const EstimationResult& a, const EstimationResult& b) {
                                              return a.error_metric < b.error_metric;
                                          });
            
            results.best_rmse = best_it->error_metric;
            results.best_parameters = best_it->parameters;
            results.best_initial_conditions = best_it->initial_conditions;
            results.success = true;
        }
        
    } catch (const std::exception& e) {
        results.error_message = "Exception during test execution: " + std::string(e.what());
        results.success = false;
    }
    
    return results;
}

std::map<std::string, ModelTestExecutor::TestResults> 
ModelTestExecutor::execute_category_tests(const std::string& category,
                                         PolynomialSolver& solver,
                                         const TestConfiguration& config) {
    std::map<std::string, TestResults> results;
    auto models = registry_.get_models_in_category(category);
    
    for (const auto& model_name : models) {
        results[model_name] = execute_test(model_name, solver, config);
    }
    
    return results;
}

poly_ode::ExperimentalData 
ModelTestExecutor::generate_test_data(const poly_ode::test_utils::OdeSystemTestBuilder& builder,
                                     const ModelTestMetadata& metadata,
                                     const TestConfiguration& config) {
    // Determine time interval
    std::pair<double, double> time_interval{0.0, 5.0};  // Default
    if (config.use_recommended_time_interval && metadata.recommended_time_interval) {
        time_interval = *metadata.recommended_time_interval;
    }
    
    // Generate time points
    std::vector<double> time_points;
    double dt = (time_interval.second - time_interval.first) / (config.num_time_points - 1);
    for (int i = 0; i < config.num_time_points; ++i) {
        time_points.push_back(time_interval.first + i * dt);
    }
    
    // Generate data
    double integration_dt = std::min(metadata.recommended_dt, dt / 10.0);
    return builder.generate_data(time_points, metadata.noise_level, integration_dt);
}

bool ModelTestExecutor::validate_identifiability_results(const poly_ode::EstimationSetupData& setup_data,
                                                        const ModelTestMetadata& metadata,
                                                        TestResults& results) {
    // Check expected identifiable count
    if (metadata.expected_identifiable_count >= 0) {
        if (static_cast<int>(setup_data.identifiable_parameters.size()) != metadata.expected_identifiable_count) {
            results.error_message = "Expected " + std::to_string(metadata.expected_identifiable_count) + 
                                   " identifiable parameters, but got " + std::to_string(setup_data.identifiable_parameters.size());
            return false;
        }
    }
    
    // Check expected unidentifiable count
    if (metadata.expected_unidentifiable_count >= 0) {
        if (static_cast<int>(setup_data.non_identifiable_parameters.size()) != metadata.expected_unidentifiable_count) {
            results.error_message = "Expected " + std::to_string(metadata.expected_unidentifiable_count) + 
                                   " unidentifiable parameters, but got " + std::to_string(setup_data.non_identifiable_parameters.size());
            return false;
        }
    }
    
    // Check specific identifiable parameters
    for (const auto& expected_param : metadata.expected_identifiable_params) {
        bool found = false;
        for (const auto& param : setup_data.identifiable_parameters) {
            if (param.name == expected_param) {
                found = true;
                break;
            }
        }
        if (!found) {
            results.error_message = "Expected parameter '" + expected_param + "' to be identifiable";
            return false;
        }
    }
    
    // Check specific unidentifiable parameters
    for (const auto& expected_param : metadata.expected_unidentifiable_params) {
        bool found = false;
        for (const auto& [param, value] : setup_data.non_identifiable_parameters) {
            if (param.name == expected_param) {
                found = true;
                break;
            }
        }
        if (!found) {
            results.error_message = "Expected parameter '" + expected_param + "' to be unidentifiable";
            return false;
        }
    }
    
    return true;
}

bool ModelTestExecutor::validate_estimation_results(const std::vector<poly_ode::EstimationResult>& estimation_results,
                                                   const poly_ode::test_utils::OdeSystemTestBuilder& builder,
                                                   const ModelTestMetadata& metadata,
                                                   TestResults& results) {
    if (estimation_results.empty()) {
        results.error_message = "No valid estimation results to validate";
        return false;
    }
    
    const auto& true_values = builder.get_true_parameter_values();
    bool found_matching_solution = false;
    
    for (const auto& result : estimation_results) {
        bool current_solution_matches = true;
        
        // Check parameters
        for (const auto& [param, estimated_value] : result.parameters) {
            auto true_it = true_values.find(param);
            if (true_it != true_values.end()) {
                if (std::abs(estimated_value - true_it->second) > metadata.parameter_tolerance) {
                    current_solution_matches = false;
                    break;
                }
            }
        }
        
        // Check initial conditions
        if (current_solution_matches) {
            for (const auto& [ic, estimated_value] : result.initial_conditions) {
                auto true_it = true_values.find(ic);
                if (true_it != true_values.end()) {
                    if (std::abs(estimated_value - true_it->second) > metadata.ic_tolerance) {
                        current_solution_matches = false;
                        break;
                    }
                }
            }
        }
        
        if (current_solution_matches) {
            found_matching_solution = true;
            break;
        }
    }
    
    if (!found_matching_solution) {
        results.error_message = "No estimation result matched the true parameter values within tolerance";
        return false;
    }
    
    return true;
}

// Metadata builders implementation
namespace metadata_builders {

ModelTestMetadata create_identifiability_metadata(const std::string& name,
                                                 const std::string& description,
                                                 int expected_identifiable_count,
                                                 const std::vector<std::string>& identifiable_params,
                                                 const std::vector<std::string>& unidentifiable_params) {
    ModelTestMetadata metadata;
    metadata.name = name;
    metadata.category = "identifiability";
    metadata.description = description;
    metadata.expected_identifiable_count = expected_identifiable_count;
    metadata.expected_identifiable_params = identifiable_params;
    metadata.expected_unidentifiable_params = unidentifiable_params;
    metadata.recommended_time_interval = std::make_pair(0.0, 5.0);
    metadata.parameter_tolerance = 1e-3;  // More relaxed for identifiability tests
    metadata.ic_tolerance = 1e-3;
    return metadata;
}

ModelTestMetadata create_classical_metadata(const std::string& name,
                                           const std::string& description,
                                           const std::pair<double, double>& time_interval,
                                           bool is_stiff) {
    ModelTestMetadata metadata;
    metadata.name = name;
    metadata.category = "classical";
    metadata.description = description;
    metadata.recommended_time_interval = time_interval;
    metadata.is_stiff = is_stiff;
    metadata.expected_identifiable_count = -1;  // Don't check by default
    return metadata;
}

ModelTestMetadata create_biological_metadata(const std::string& name,
                                            const std::string& description,
                                            const std::pair<double, double>& time_interval,
                                            bool has_rational_functions) {
    ModelTestMetadata metadata;
    metadata.name = name;
    metadata.category = "biological";
    metadata.description = description;
    metadata.recommended_time_interval = time_interval;
    metadata.has_rational_functions = has_rational_functions;
    metadata.expected_identifiable_count = -1;  // Don't check by default
    metadata.parameter_tolerance = 1e-2;  // More relaxed for biological systems
    metadata.ic_tolerance = 1e-2;
    return metadata;
}

ModelTestMetadata create_simple_metadata(const std::string& name,
                                        const std::string& description,
                                        int expected_identifiable_count) {
    ModelTestMetadata metadata;
    metadata.name = name;
    metadata.category = "simple";
    metadata.description = description;
    metadata.expected_identifiable_count = expected_identifiable_count;
    metadata.recommended_time_interval = std::make_pair(0.0, 5.0);
    return metadata;
}

} // namespace metadata_builders

// Global registry
static ModelTestRegistry global_registry;

ModelTestRegistry& get_global_model_registry() {
    return global_registry;
}

} // namespace test_framework
} // namespace poly_ode