#include "ode_system.hpp" // Include the new ODESystem class
#include "polynomial.hpp" // Still needed for Variable, RationalFunction etc.

#include <cmath>         // For std::floor
#include <gtest/gtest.h> // Include Google Test header
#include <iomanip>
#include <iostream>
#include <map>
#include <memory> // Include for std::unique_ptr
#include <string>
#include <vector>

// Helper function to print the results map nicely (keep for potential debugging)
template<typename Coeff>
void
print_results(const typename ODESystem<Coeff>::ResultsType &results) {
    if (results.empty() || results.find("time") == results.end() || results.at("time").empty()) {
        std::cout << "No results to print." << std::endl;
        return;
    }

    const auto &time_vec = results.at("time");
    size_t n_points = time_vec.size();

    // Print header
    std::cout << std::setw(12) << "Time";
    for (const auto &pair : results) {
        if (pair.first != "time") { std::cout << std::setw(15) << pair.first; }
    }
    std::cout << std::endl;

    // Print dashed line separator
    int total_width = 12;
    for (const auto &pair : results) {
        if (pair.first != "time") { total_width += 15; }
    }
    std::cout << std::string(total_width, '-') << std::endl;

    // Print data rows
    std::cout << std::fixed << std::setprecision(6);
    for (size_t i = 0; i < n_points; ++i) {
        std::cout << std::setw(12) << time_vec[i];
        for (const auto &pair : results) {
            if (pair.first != "time") {
                if (pair.second.size() > i) {
                    std::cout << std::setw(15) << pair.second[i];
                } else {
                    std::cout << std::setw(15) << "N/A"; // Should not happen if observer works
                }
            }
        }
        std::cout << std::endl;
    }
}


// Define a test fixture for common setup
class HollingODESystemTest : public ::testing::Test {
  protected:
    using Coeff = double;
    using StateType = ODESystemStateType<Coeff>;

    // Define shared variables, parameters, system etc.
    const Variable x_var{ "x" };
    const Variable y_var{ "y" };
    const Variable r_var{ "r", 0, true };
    const Variable K_var{ "K", 0, true };
    const Variable a_var{ "a", 0, true };
    const Variable b_var{ "b", 0, true };
    const Variable c_var{ "c", 0, true };
    const Variable d_var{ "d", 0, true };

    std::vector<Variable> state_vars;
    std::vector<Variable> param_vars;
    std::vector<RationalFunction<Coeff>> rhs_equations;
    std::map<std::string, RationalFunction<Coeff>> observables;
    // Hold the system by unique_ptr
    std::unique_ptr<ODESystem<Coeff>> holling_system;

    // Setup runs before each test in the fixture
    void SetUp() override {
        state_vars = { x_var, y_var };
        param_vars = { r_var, K_var, a_var, b_var, c_var, d_var };

        auto holling_response = (a_var * x_var) / (1.0 + b_var * x_var);
        auto logistic_growth = r_var * x_var - (r_var * x_var * x_var) / K_var;
        auto predation = holling_response * y_var;
        auto dx_dt_rf = logistic_growth - predation;

        auto predator_growth = c_var * predation;
        auto predator_death = d_var * y_var;
        auto dy_dt_rf = predator_growth - predator_death;

        rhs_equations = { dx_dt_rf, dy_dt_rf };

        observables = { { "Sum (x+y)", RationalFunction<Coeff>(x_var + y_var) },
                        { "Product (x*y)", RationalFunction<Coeff>(x_var * y_var) },
                        { "Holling Term", holling_response } };

        // Construct the system using make_unique in SetUp
        holling_system = std::make_unique<ODESystem<Coeff>>(state_vars, rhs_equations, param_vars, observables);
    }

    // Helper to validate simulation results structure
    void ValidateResults(const typename ODESystem<Coeff>::ResultsType &results,
                         double t_start,
                         double t_end,
                         double dt_observe) {
        ASSERT_FALSE(results.empty());
        ASSERT_TRUE(results.count("time"));
        ASSERT_TRUE(results.count("x"));
        ASSERT_TRUE(results.count("y"));
        ASSERT_TRUE(results.count("Sum (x+y)"));
        ASSERT_TRUE(results.count("Product (x*y)"));
        ASSERT_TRUE(results.count("Holling Term"));

        ASSERT_FALSE(results.at("time").empty());
        size_t expected_points = static_cast<size_t>(std::floor((t_end - t_start) / dt_observe)) + 1;
        // Allow for slight variations due to floating point end condition
        ASSERT_GE(results.at("time").size(), expected_points - 1);
        ASSERT_LE(results.at("time").size(), expected_points + 1);

        size_t n_points = results.at("time").size();
        ASSERT_GT(n_points, 0); // Should have at least one point

        for (const auto &pair : results) {
            ASSERT_EQ(pair.second.size(), n_points) << "Vector size mismatch for key: " << pair.first;
        }
    }
};

// Test case for a single simulation run using the fixture
TEST_F(HollingODESystemTest, SingleSimulationRun) {
    // Parameter values
    std::map<Variable, Coeff> parameter_values = { { r_var, 1.0 }, { K_var, 10.0 }, { a_var, 1.0 },
                                                   { b_var, 0.2 }, { c_var, 0.5 },  { d_var, 0.2 } };

    // Initial conditions
    StateType initial_conditions = { 5.0, 2.0 }; // x=5, y=2

    // Simulation time settings
    double t_start = 0.0;
    double t_end = 100.0;
    double dt_observe = 0.5;

    typename ODESystem<Coeff>::ResultsType results;
    ASSERT_NO_THROW({
        // Use -> to access simulate via the pointer
        results = holling_system->simulate(initial_conditions, parameter_values, t_start, t_end, dt_observe);
    }) << "Single simulation threw an exception.";

    // Validate the structure and basic properties of the results
    ValidateResults(results, t_start, t_end, dt_observe);
}

// Test case for a batch simulation run using the fixture
TEST_F(HollingODESystemTest, BatchSimulationRun) {
    std::vector<StateType> ic_list = { { 5.0, 2.0 }, { 1.0, 1.0 }, { 9.0, 4.0 } };
    std::map<Variable, Coeff> base_params = { { r_var, 1.0 }, { K_var, 10.0 }, { a_var, 1.0 },
                                              { b_var, 0.2 }, { c_var, 0.5 },  { d_var, 0.2 } };
    std::map<Variable, Coeff> varied_params = { { r_var, 1.2 }, { K_var, 15.0 }, { a_var, 0.8 },
                                                { b_var, 0.3 }, { c_var, 0.6 },  { d_var, 0.25 } };
    std::vector<std::map<Variable, Coeff>> param_list = { base_params, base_params, varied_params };

    double t_start = 0.0;
    double t_end = 50.0; // Shorter time for batch
    double dt_observe = 1.0;

    std::vector<typename ODESystem<Coeff>::ResultsType> batch_results;
    ASSERT_NO_THROW({
        // Use -> to access simulate_batch via the pointer
        batch_results = holling_system->simulate_batch(ic_list, param_list, t_start, t_end, dt_observe);
    }) << "Batch simulation threw an exception.";

    // Check that we got results for each input condition
    ASSERT_EQ(batch_results.size(), ic_list.size());

    // Validate each result set in the batch
    for (const auto &results : batch_results) { ValidateResults(results, t_start, t_end, dt_observe); }
}

// --- Original main() removed ---
/*
int
main() {
    // ... original main contents ...
}
*/