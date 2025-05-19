#include "MSolveSolver.hpp" // For class definition and includes like polynomial_solver.hpp, msolve/msolve/msolve.h (for types)

#include "algebraic_system.hpp"
#include "polynomial.hpp" // For Variable, Polynomial, Monomial etc.

#include <algorithm> // For std::remove_if, std::transform
#include <cstdio>    // For open_memstream, fileno, FILE
#include <cstdlib>   // For atof, free
#include <cstring>   // For strdup, strtok_r, strlen
#include <iostream>  // For std::cout, std::cerr
#include <map>
#include <sstream>   // For std::stringstream, std::ostringstream
#include <stdexcept> // For std::runtime_error
#include <vector>

namespace poly_ode {

// Initialize static members
bool MSolveSolver::msolve_initialized_ = false;
int MSolveSolver::instance_count_ = 0;

void
MSolveSolver::initialize_msolve_if_needed() {
    if (instance_count_ == 0) {
        std::cout << "  [MSolveSolver] Calling msolve_global_init()..." << std::endl;
        msolve_global_init();       // Provided by user notes
        msolve_initialized_ = true; // Assuming success
    }
    instance_count_++;
}

void
MSolveSolver::uninitialize_msolve_if_needed() {
    instance_count_--;
    if (instance_count_ == 0 && msolve_initialized_) {
        // If msolve had a msolve_global_cleanup(), it would be called here.
        // Based on the user's notes, only msolve_global_init() is mentioned.
        // jl_atexit_hook(0); is for the Julia example, not the C API directly.
        std::cout << "  [MSolveSolver] Last instance destroyed. msolve resources assumed released or managed by OS."
                  << std::endl;
        msolve_initialized_ = false;
    }
}

MSolveSolver::MSolveSolver() {
    initialize_msolve_if_needed();
    std::cout << "MSolveSolver instance created." << std::endl;
}

MSolveSolver::~MSolveSolver() {
    uninitialize_msolve_if_needed();
    std::cout << "MSolveSolver instance destroyed." << std::endl;
}

std::string
MSolveSolver::name() const {
    return "MSolveSolver";
}

// Helper to convert poly_ode::Polynomial to msolve string part
std::string
polynomial_to_msolve_string(const Polynomial<double> &poly, const std::map<Variable, std::string> &var_to_msolve_name) {
    std::stringstream ss;
    bool first_term = true;
    for (const auto &mono : poly.monomials) {
        if (!first_term && mono.coeff >= 0) { ss << "+"; }
        // msolve expects explicit coefficients, e.g., 1*x, -1*x
        ss << mono.coeff;

        // Iterate directly over mono.vars, which is a std::map and iterates in sorted order by Variable key
        for (const auto &var_exponent_pair : mono.vars) { // var_exponent_pair is std::pair<const Variable, int>
            auto it = var_to_msolve_name.find(var_exponent_pair.first);
            if (it != var_to_msolve_name.end()) {
                ss << "*" << it->second;
                if (var_exponent_pair.second > 1) { ss << "^" << var_exponent_pair.second; }
            }
        }
        first_term = false;
    }
    if (poly.monomials.empty()) {
        ss << "0"; // Represent a zero polynomial
    }
    return ss.str();
}

std::string
MSolveSolver::convert_to_msolve_format(const AlgebraicSystem &system,
                                       std::map<Variable, std::string> &var_to_msolve_name) const {
    std::stringstream ss_vars;
    std::vector<std::string> msolve_var_names;

    // Create a mapping from poly_ode::Variable to simple msolve names (e.g., x0, x1, ... or use original names if
    // simple) Using original names if they are simple enough for msolve, otherwise using indexed names For now, let's
    // try to use original names directly, assuming they are valid C identifiers. Msolve example uses "x,y", so spaces
    // or special chars in poly_ode var names might be an issue. Let's stick to original names for now, and sanitize if
    // it becomes a problem.

    for (size_t i = 0; i < system.unknowns.size(); ++i) {
        const auto &var = system.unknowns[i];
        std::string msolve_name = var.name;
        // Basic sanitization: remove spaces, replace special chars if necessary.
        // For now, just using the name as is. Msolve might be flexible.
        // Example: "x_1" -> "x_1"
        // We need to ensure these names are consistent with what msolve_solution_print or other accessors would use.
        var_to_msolve_name[var] = msolve_name;
        msolve_var_names.push_back(msolve_name);
    }

    for (size_t i = 0; i < msolve_var_names.size(); ++i) {
        ss_vars << msolve_var_names[i] << (i == msolve_var_names.size() - 1 ? "" : ",");
    }

    std::stringstream ss_polys;
    for (size_t i = 0; i < system.polynomials.size(); ++i) {
        ss_polys << polynomial_to_msolve_string(system.polynomials[i], var_to_msolve_name)
                 << (i == system.polynomials.size() - 1 ? "" : ", ");
    }

    // Format: "vars; 0; polys"
    return ss_vars.str() + "; 0; " + ss_polys.str();
}


PolynomialSolutionSet
MSolveSolver::solve(const AlgebraicSystem &system) {
    initialize_msolve_if_needed();
    PolynomialSolutionSet solutions;
    std::map<Variable, std::string> var_to_msolve_name; // Map to store msolve internal var names

    std::string system_str = convert_to_msolve_format(system, var_to_msolve_name);
    if (system_str.empty()) {
        std::cerr << "[MSolveSolver] Error: System string is empty after conversion." << std::endl;
        return {};
    }

    std::cout << "  [MSolveSolver] Generated msolve input:\n" << system_str << std::endl;

    // Get a list of unique variable names in the order they appear in AlgebraicSystem unknowns
    std::vector<std::string> var_names;
    for (const auto &unknown : system.unknowns) {
        var_names.push_back(unknown.name); // Assuming Variable has a .name field
    }

    struct msolve_system_s *S = msolve_system_new(static_cast<unsigned int>(var_names.size())); // Restore argument
    if (!S) {
        std::cerr << "[MSolveSolver] Error: msolve_system_new failed." << std::endl;
        return {};
    }

    if (msolve_system_parse(S, system_str.c_str()) != 0) { // Check for msolve_ok or similar if available
        std::cerr << "[MSolveSolver] Error: msolve_system_parse failed: " << msolve_system_get_error(S) << std::endl;
        msolve_system_free(S);
        return {};
    }

    msolve_solution_array_t *roots = msolve_solve(S);

    if (roots) {
        size_t num_solutions = msolve_solution_array_size(roots);
        std::cout << "    [MSolveSolver] Found " << num_solutions << " solution(s) (real or complex)." << std::endl;

        for (size_t i = 0; i < num_solutions; ++i) {
            // The user example uses msolve_solution_print. We need to extract values.
            // The `interfaces/c/` directory in msolve repo would have examples for this.
            // Let's assume a hypothetical msolve_solution_get_value(roots, i, var_name_str, &real, &imag)
            // Or, msolve_solution_print_to_buffer and parse that buffer.
            // For now, trying to adapt based on the print function structure.

            // Alternative: msolve_solution_fprintf with a string stream, then parse.
            // This is a common pattern if direct value access is not straightforward.
            char *sol_buffer = nullptr;
            size_t sol_buffer_size = 0;
            FILE *mem_stream = open_memstream(&sol_buffer, &sol_buffer_size);
            if (mem_stream) {
                msolve_solution_print(mem_stream, roots, i); // Prints one solution
                fclose(mem_stream);                          // Flushes and finalizes buffer

                if (sol_buffer) {
                    std::cout << "      Raw Solution " << i << ": " << sol_buffer << std::endl;
                    PolynomialSolutionMap current_solution_map =
                      parse_msolve_solution_string(sol_buffer, var_to_msolve_name, system.unknowns);
                    if (!current_solution_map.empty()) { solutions.push_back(current_solution_map); }
                    free(sol_buffer); // Free memory allocated by open_memstream
                }
            } else {
                std::cerr << "    [MSolveSolver] Warning: Could not open memstream to capture solution " << i
                          << std::endl;
            }
        }
        msolve_solution_array_free(roots);
    } else {
        std::cerr << "    [MSolveSolver] msolve_solve returned nullptr (no solutions or error)." << std::endl;
        // Check for system error after solve if roots is null
        const char *error_msg = msolve_system_get_error(S);
        if (error_msg) { std::cerr << "    [MSolveSolver] Error after solve: " << error_msg << std::endl; }
    }

    msolve_system_free(S);

    std::cout << "  [MSolveSolver] Solve finished. Processed " << solutions.size() << " valid solution(s)."
              << std::endl;
    return solutions;
}

// Helper to parse a single solution string from msolve_solution_print output.
// Example output from msolve_solution_print for one solution:
// "x = 1.234, y = 5.678 (real solution)" or "x = 1.2+3.4i, y = ... (complex solution)"
PolynomialSolutionMap
MSolveSolver::parse_msolve_solution_string(const char *solution_str_c,
                                           const std::map<Variable, std::string> &var_to_msolve_name,
                                           const std::vector<Variable> &system_unknowns) const {
    PolynomialSolutionMap sol_map;
    if (!solution_str_c) return sol_map;

    std::string solution_str(solution_str_c);
    std::cout << "    [MSolveSolver-Parser] Parsing solution string: '" << solution_str << "'" << std::endl;

    // This is a placeholder parsing strategy. Msolve's actual output format needs to be handled robustly.
    // The example `msolve_solution_print(stdout, roots, i);` implies a human-readable format.
    // We need to map back from msolve variable names to poly_ode::Variable objects.

    // Create reverse map: msolve_name_str -> poly_ode::Variable
    std::map<std::string, Variable> msolve_name_to_var_obj;
    for (const auto &pair : var_to_msolve_name) { msolve_name_to_var_obj[pair.second] = pair.first; }

    // Tokenize by comma for assignments like "var1 = val1, var2 = val2"
    char *mutable_sol_str = strdup(solution_str_c);
    if (!mutable_sol_str) return sol_map;

    char *saveptr1;
    char *token = strtok_r(mutable_sol_str, ",", &saveptr1);

    while (token != nullptr) {
        std::string assignment(token);
        // Trim whitespace
        assignment.erase(assignment.begin(), std::find_if(assignment.begin(), assignment.end(), [](unsigned char ch) {
                             return !std::isspace(ch);
                         }));
        assignment.erase(
          std::find_if(assignment.rbegin(), assignment.rend(), [](unsigned char ch) { return !std::isspace(ch); })
            .base(),
          assignment.end());

        // Split assignment by "="
        size_t eq_pos = assignment.find('=');
        if (eq_pos != std::string::npos) {
            std::string var_name_msolve = assignment.substr(0, eq_pos);
            std::string val_str = assignment.substr(eq_pos + 1);

            // Trim var_name_msolve and val_str
            var_name_msolve.erase(var_name_msolve.begin(),
                                  std::find_if(var_name_msolve.begin(), var_name_msolve.end(), [](unsigned char ch) {
                                      return !std::isspace(ch);
                                  }));
            var_name_msolve.erase(std::find_if(var_name_msolve.rbegin(),
                                               var_name_msolve.rend(),
                                               [](unsigned char ch) { return !std::isspace(ch); })
                                    .base(),
                                  var_name_msolve.end());

            val_str.erase(val_str.begin(), std::find_if(val_str.begin(), val_str.end(), [](unsigned char ch) {
                              return !std::isspace(ch);
                          }));
            val_str.erase(
              std::find_if(val_str.rbegin(), val_str.rend(), [](unsigned char ch) { return !std::isspace(ch); }).base(),
              val_str.end());

            // Check if this msolve variable name is one we are tracking
            auto poly_ode_var_it = msolve_name_to_var_obj.find(var_name_msolve);
            if (poly_ode_var_it != msolve_name_to_var_obj.end()) {
                Variable current_poly_ode_var = poly_ode_var_it->second;
                double real_part = 0.0, imag_part = 0.0;

                // Attempt to parse complex number like "1.2+3.4i" or "1.2-3.4i" or just "1.2"
                // This is a simplified parser.
                size_t i_pos = val_str.find('i');
                size_t plus_pos = val_str.find('+');
                size_t minus_pos = val_str.find('-'); // Could be sign of real or separator for imag

                if (i_pos != std::string::npos) { // Complex number
                    std::string real_str, imag_str;
                    bool imag_is_signed_explicitly = false;

                    // Try to find the split point between real and imaginary part based on '+' or '-' not at the start
                    size_t split_pos = std::string::npos;
                    if (plus_pos != std::string::npos && plus_pos > 0)
                        split_pos = plus_pos;
                    else if (minus_pos != std::string::npos && minus_pos > 0)
                        split_pos = minus_pos;

                    if (split_pos != std::string::npos) {
                        real_str = val_str.substr(0, split_pos);
                        imag_str = val_str.substr(split_pos, i_pos - split_pos);
                        // If imag_str is just "+" or "-", it means +/- 1.0i, but atof needs a number.
                        if (imag_str == "+" || imag_str == "-")
                            imag_str += "1.0";
                        else if (imag_str.empty() && val_str[i_pos - 1] == '+')
                            imag_str = "1.0"; // e.g. 2+i
                        else if (imag_str.empty() && val_str[i_pos - 1] == '-')
                            imag_str = "-1.0"; // e.g. 2-i

                    } else { // Pure imaginary e.g. "3.4i" or "-3.4i" or "i" or "-i"
                        real_str = "0.0";
                        imag_str = val_str.substr(0, i_pos);
                        if (imag_str.empty() || imag_str == "+")
                            imag_str = "1.0";
                        else if (imag_str == "-")
                            imag_str = "-1.0";
                    }

                    try {
                        real_part = std::stod(real_str);
                        imag_part = std::stod(imag_str); // stod should handle signs correctly
                    } catch (const std::exception &e) {
                        std::cerr << "    [MSolveSolver-Parser] Error parsing complex parts for " << var_name_msolve
                                  << ": '" << val_str << "' -> real_str='" << real_str << "', imag_str='" << imag_str
                                  << "'. Error: " << e.what() << std::endl;
                        continue; // Skip this variable assignment
                    }
                } else { // Real number
                    // Remove any trailing text like "(real solution)"
                    size_t paren_pos = val_str.find('(');
                    if (paren_pos != std::string::npos) {
                        val_str = val_str.substr(0, paren_pos);
                        // Trim trailing whitespace again
                        val_str.erase(std::find_if(val_str.rbegin(),
                                                   val_str.rend(),
                                                   [](unsigned char ch) { return !std::isspace(ch); })
                                        .base(),
                                      val_str.end());
                    }
                    try {
                        real_part = std::stod(val_str);
                    } catch (const std::exception &e) {
                        std::cerr << "    [MSolveSolver-Parser] Error parsing real part for " << var_name_msolve
                                  << ": '" << val_str << "'. Error: " << e.what() << std::endl;
                        continue; // Skip this variable assignment
                    }
                }
                sol_map[current_poly_ode_var] = std::complex<double>(real_part, imag_part);
                std::cout << "      [MSolveSolver-Parser] Parsed: " << current_poly_ode_var << " = " << real_part
                          << " + " << imag_part << "i" << std::endl;
            }
        }
        token = strtok_r(nullptr, ",", &saveptr1);
    }
    free(mutable_sol_str);

    // Check if all system unknowns were found in the solution string
    if (sol_map.size() != system_unknowns.size() && !system_unknowns.empty() && !solution_str.empty()) {
        std::cerr << "    [MSolveSolver-Parser] Warning: Parsed solution map size (" << sol_map.size()
                  << ") does not match system unknowns size (" << system_unknowns.size() << ") for solution string: '"
                  << solution_str_c << "'" << std::endl;
        // This might indicate a parsing error or an incomplete solution from msolve for this specific solution entry.
        // Depending on strictness, one might clear sol_map here to invalidate it.
    }
    return sol_map;
}


} // namespace poly_ode