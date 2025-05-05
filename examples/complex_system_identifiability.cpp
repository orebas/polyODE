/**
 * Complex System Identifiability Analysis Example
 * ------------------------------------------------
 *
 * This example demonstrates the use of the IdentifiabilityAnalyzer on a complex
 * system with multiple states, parameters, and observables. It analyzes which
 * model parameters and initial conditions are locally identifiable and determines
 * the minimum derivative orders required for identifiability.
 *
 * System:
 *   x1_dot = x2
 *   x2_dot = x3
 *   x3_dot = theta*x1 + x2*x2
 *   w_dot = -(k1+k2)*w
 *   z_dot = k*k*z
 *   u_dot = -u
 *
 * Observables:
 *   y1 = x1
 *   y2 = w
 *   y3 = z
 */

#include "identifiability_analyzer.hpp"
#include "observed_ode_system.hpp"
#include "polynomial.hpp"
#include "rational_function_operators.hpp"

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace poly_ode;

// Helper function to print a separator line
void
printSeparator(char c = '-', int width = 80) {
    std::cout << std::string(width, c) << std::endl;
}

// Helper function to print the results in a human-readable format
void
printResults(const IdentifiabilityAnalyzer::AnalysisResults &results) {
    printSeparator('=');
    std::cout << "IDENTIFIABILITY ANALYSIS RESULTS" << std::endl;
    printSeparator('=');

    // Print identifiable parameters
    std::cout << "IDENTIFIABLE PARAMETERS (" << results.identifiable_parameters.size() << "):" << std::endl;
    printSeparator('-');
    if (results.identifiable_parameters.empty()) {
        std::cout << "  None found" << std::endl;
    } else {
        for (const auto &param : results.identifiable_parameters) { std::cout << "  " << param << std::endl; }
    }
    std::cout << std::endl;

    // Print non-identifiable parameters and their fixed values
    std::cout << "NON-IDENTIFIABLE PARAMETERS (" << results.non_identifiable_parameters.size() << "):" << std::endl;
    printSeparator('-');
    if (results.non_identifiable_parameters.empty()) {
        std::cout << "  None found" << std::endl;
    } else {
        std::cout << "  Parameter         Fixed Value" << std::endl;
        printSeparator('-', 40);
        for (const auto &pair : results.non_identifiable_parameters) {
            std::cout << "  " << std::left << std::setw(18) << pair.first.name << std::fixed << std::setprecision(6)
                      << pair.second << std::endl;
        }
    }
    std::cout << std::endl;

    // Print required derivative orders for each observable
    std::cout << "REQUIRED DERIVATIVE ORDERS:" << std::endl;
    printSeparator('-');
    if (results.required_derivative_orders.empty()) {
        std::cout << "  Not determined" << std::endl;
    } else {
        std::cout << "  Observable        Min Order" << std::endl;
        printSeparator('-', 40);
        for (const auto &pair : results.required_derivative_orders) {
            std::cout << "  " << std::left << std::setw(18) << pair.first.name << pair.second << std::endl;
        }
    }
    printSeparator('=');
}

int
main() {
    try {
        std::cout << "Complex System Identifiability Analysis" << std::endl;
        printSeparator();

        // Define the system variables
        Variable x1("x1");
        Variable x2("x2");
        Variable x3("x3");
        Variable w("w");
        Variable z("z");
        Variable u("u");
        Variable theta("theta", 0, true);
        Variable k1("k1", 0, true);
        Variable k2("k2", 0, true);
        Variable k("k", 0, true);

        // Define the observables
        Observable y1_obs("y1");
        Observable y2_obs("y2");
        Observable y3_obs("y3");

        // Create the ODE system
        std::cout << "Setting up the ODE system..." << std::endl;

        ObservedOdeSystem complex_system;
        complex_system.state_variables = { x1, x2, x3, w, z, u };
        complex_system.parameters = { theta, k1, k2, k };

        // Define the differential equations using the new operator overloads
        complex_system.equations = { // x1_dot = x2
                                     x2,

                                     // x2_dot = x3
                                     x3,

                                     // x3_dot = theta*x1 + x2*x2
                                     theta * x1 + x2 * x2,

                                     // w_dot = -(k1+k2)*w
                                     -(k1 + k2) * w,

                                     // z_dot = k*k*z
                                     k * k * z,

                                     // u_dot = -u
                                     -u
        };

        // Define the observable equations
        complex_system.observable_definitions = { { y1_obs, x1 }, { y2_obs, w }, { y3_obs, z } };

        // Print summary of the system
        std::cout << "Complex System Summary:" << std::endl;
        std::cout << "  States: " << complex_system.state_variables.size() << std::endl;
        std::cout << "  Parameters: " << complex_system.parameters.size() << std::endl;
        std::cout << "  Observables: " << complex_system.observable_definitions.size() << std::endl;

        // Create the analyzer including ALL parameters and initial conditions
        std::vector<Variable> params_to_analyze = { theta, k1, k2, k, x1, x2, x3, w, z, u }; // 10 total
        int max_deriv_order = 3;  // Need up to order 3 for complete identifiability
        int num_test_points = 10; // Use multiple test points for robust analysis

        std::cout << "\nSetting up identifiability analysis..." << std::endl;
        std::cout << "  Max derivative order: " << max_deriv_order << std::endl;
        std::cout << "  Number of test points: " << num_test_points << std::endl;
        std::cout << "  Parameters & initial conditions to analyze: " << params_to_analyze.size() << std::endl;

        // Create the analyzer and run the analysis
        IdentifiabilityAnalyzer analyzer(complex_system, params_to_analyze, max_deriv_order);

        std::cout << "\nRunning identifiability analysis..." << std::endl;
        auto results = analyzer.analyze(num_test_points, 1e-12); // Use tight tolerance for rank determination

        // Print the results in a nicely formatted way
        printResults(results);

        // Provide interpretation of results
        std::cout << "\nINTERPRETATION:" << std::endl;
        printSeparator('-');
        std::cout << "• The system has " << results.identifiable_parameters.size() << " identifiable parameters/ICs."
                  << std::endl;
        std::cout << "• The non-identifiable parameters represent structural non-identifiabilities." << std::endl;

        bool found_k_sum = false;
        for (const auto &pair : results.non_identifiable_parameters) {
            if (pair.first == k1 || pair.first == k2) {
                found_k_sum = true;
                const auto &other_k = (pair.first == k1) ? k2 : k1;
                std::cout << "• Only the sum (k1+k2) is identifiable, not individual values of k1 and k2." << std::endl;
                std::cout << "  In this analysis, " << pair.first << " was fixed, making " << other_k
                          << " identifiable." << std::endl;
            }
            if (pair.first == u) {
                std::cout << "• The u parameter/state is completely unidentifiable. It doesn't appear in any "
                             "observable equation or their derivatives."
                          << std::endl;
            }
        }

        std::cout << "• For complete identifiability, we need:" << std::endl;
        for (const auto &pair : results.required_derivative_orders) {
            std::cout << "  - " << pair.first.name << " observable: derivatives up to order " << pair.second
                      << std::endl;
            if (pair.first.name == "y1" && pair.second == 3) {
                std::cout << "    (This means we need 3rd derivatives of y1 to identify the parameter theta)"
                          << std::endl;
            }
            if (pair.first.name == "y2" && pair.second == 1) {
                std::cout << "    (This means we need 1st derivatives of y2 to identify the (k1+k2) combination)"
                          << std::endl;
            }
            if (pair.first.name == "y3" && pair.second == 1) {
                std::cout << "    (This means we need 1st derivatives of y3 to identify the parameter k)" << std::endl;
            }
        }

        return 0;
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Unknown error occurred" << std::endl;
        return 1;
    }
}