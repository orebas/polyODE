#include "MSolveSolver.hpp" // For class definition and includes like polynomial_solver.hpp, msolve/msolve/msolve.h (for types)

// Attempt to prioritize vcpkg flint includes for this translation unit
#include <flint/acb_poly.h> // For acb_poly_t (expecting vcpkg path)
#include <flint/arb.h>      // For arb_t (expecting vcpkg path)
#include <flint/flint.h>    // General flint utilities (expecting vcpkg path)
#include <flint/fmpq.h>     // For Flint fmpq_t operations (expecting vcpkg path)
#include <flint/fmpz.h>     // Often needed with fmpq
#include <gmp.h>            // GMP must come before any FLINT header

#include "univariate_solver_flint.hpp" // Includes <flint/arb.h> (hopefully now resolves consistently)

#include "algebraic_system.hpp"
#include "polynomial.hpp" // For Variable, Polynomial, Monomial etc.

#include <algorithm> // For std::remove_if, std::transform
#include <array>     // For std::array (buffer for popen)
#include <cmath>     // For std::gcd, std::round, std::abs
#include <cstdio>    // For open_memstream, fileno, FILE
#include <cstdio>    // For popen, pclose, fgets
#include <cstdlib>   // For atof, free
#include <cstring>   // For strdup, strtok_r, strlen

#include <fstream>  // For file I/O
#include <iostream> // For std::cout, std::cerr
#include <map>
#include <nlohmann/json.hpp> // For JSON parsing - ensure this is available
#include <random>            // For generating unique temp file names
#include <regex>             // For std::regex_replace
#include <sstream>           // For std::stringstream, std::ostringstream
#include <stdexcept>         // For std::runtime_error
#include <unistd.h>          // For access(), to check if msolve binary exists
#include <vector>


#ifdef __has_include
#if __has_include(<unsupported/Eigen/Polynomials>)
#include <unsupported/Eigen/Polynomials>
#elif __has_include(<eigen3/unsupported/Eigen/Polynomials>)
#include <eigen3/unsupported/Eigen/Polynomials>
#else
#error "Eigen Polynomials header not found"
#endif
#else
#include <unsupported/Eigen/Polynomials>
#endif

namespace poly_ode {

// Initialize static members
bool MSolveSolver::msolve_initialized_ = false;
int MSolveSolver::instance_count_ = 0;

// Helper function to generate unique temporary filenames
std::string
generate_temp_filename(const std::string &prefix) {
    // Create a random string to avoid collisions
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_int_distribution<> dis(0, 35); // a-z, 0-9

    std::string random_str;
    for (int i = 0; i < 8; i++) {
        int r = dis(gen);
        if (r < 26)
            random_str += 'a' + r;
        else
            random_str += '0' + (r - 26);
    }

    return "/tmp/" + prefix + "_" + random_str;
}

void
MSolveSolver::initialize_msolve_if_needed() {
    if (instance_count_ == 0) {
        std::cout << "  [MSolveSolver] Checking msolve binary..." << std::endl;

        // Check if msolve binary is accessible
        const char *msolve_paths[] = {
            "msolve",                                                  // In PATH
            "build/vcpkg_installed/x64-linux/tools/msolve/bin/msolve", // vcpkg install location
            "/usr/bin/msolve",                                         // System install
            "/usr/local/bin/msolve"                                    // Local install
        };

        bool found = false;
        for (const char *path : msolve_paths) {
            if (access(path, X_OK) == 0) {
                std::cout << "  [MSolveSolver] Found msolve binary at: " << path << std::endl;
                found = true;
                break;
            }
        }

        if (!found) {
            std::cerr << "  [MSolveSolver] WARNING: msolve binary not found in standard locations." << std::endl;
            std::cerr << "  [MSolveSolver] Please ensure msolve is installed and in your PATH." << std::endl;
        }

        msolve_initialized_ = true;
    }
    instance_count_++;
}

void
MSolveSolver::uninitialize_msolve_if_needed() {
    instance_count_--;
    if (instance_count_ == 0 && msolve_initialized_) {
        std::cout << "  [MSolveSolver] Last instance destroyed." << std::endl;
        msolve_initialized_ = false;
    }
}

// Default constructor
MSolveSolver::MSolveSolver()
  : msolve_executable_path_("msolve")
  , python_script_path_("scripts/msolve_p_to_json.py") {
    // Attempt to find msolve and script in common locations or via environment variables if not absolute
    // For now, uses hardcoded relative paths or assumes they are in PATH / relative to execution
    std::cout << "  [MSolveSolver] Default constructor: msolve path: " << msolve_executable_path_
              << ", script path: " << python_script_path_ << std::endl;
}

// Constructor with paths
MSolveSolver::MSolveSolver(std::string msolve_executable_path, std::string python_script_path)
  : msolve_executable_path_(std::move(msolve_executable_path))
  , python_script_path_(std::move(python_script_path)) {
    std::cout << "  [MSolveSolver] Path constructor: msolve path: " << msolve_executable_path_
              << ", script path: " << python_script_path_ << std::endl;
}

MSolveSolver::~MSolveSolver() {
    uninitialize_msolve_if_needed();
    std::cout << "MSolveSolver instance destroyed." << std::endl;
}

std::string
MSolveSolver::name() const {
    return "MSolveSolver";
}

// Static helper method to convert double to exact rational string via GMP -> FLINT
std::string
MSolveSolver::double_to_rational_string(double x) {
    // 1) GMP: double -> mpq_t
    mpq_t g;
    mpq_init(g);
    mpq_set_d(g, x);
    mpq_canonicalize(g); // Reduce to lowest terms

    // 2) FLINT: mpq_t -> fmpq_t
    fmpq_t q;
    fmpq_init(q);
    fmpq_set_mpq(q, g); // Import from GMP rational

    // Clean up GMP object
    mpq_clear(g);

    // 3) FLINT: fmpq_t -> base-10 string
    char *s = fmpq_get_str(NULL, 10, q); // NULL means flint allocates string
    std::string str(s);
    flint_free(s); // Free the string allocated by flint

    // Clean up FLINT object
    fmpq_clear(q);
    return str;
}

// Helper to convert poly_ode::Polynomial to msolve string format
std::string
polynomial_to_msolve_string(const Polynomial<double> &poly, const std::map<Variable, std::string> &var_to_msolve_name) {
    std::stringstream ss;
    bool first_term = true;
    for (const auto &mono : poly.monomials) {
        if (!first_term && mono.coeff >= 0) { ss << "+"; }

        // Use the robust double_to_rational_string converter
        ss << MSolveSolver::double_to_rational_string(mono.coeff);

        for (const auto &var_exponent_pair : mono.vars) {
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
                                       std::map<Variable, std::string> &var_to_msolve_name_map) const {
    std::stringstream ss_vars;
    std::vector<std::string> msolve_var_names_ordered; // Keep track of order for msolve output
    var_to_msolve_name_map.clear();                    // Ensure the map is clean before populating

    // Create mapping from Variables to msolve names
    for (size_t i = 0; i < system.unknowns.size(); ++i) {
        const auto &var = system.unknowns[i];

        std::string msolve_name_base;
        if (var.deriv_level > 0) {
            // Construct a name like "d<level><baseName>"
            msolve_name_base = "d" + std::to_string(var.deriv_level) + var.name;
        } else {
            msolve_name_base = var.name;
        }
        // Sanitize the base name
        msolve_name_base = std::regex_replace(
          msolve_name_base, std::regex("[^a-zA-Z0-9_]"), ""); // Keep only alphanumeric and underscore
        msolve_name_base = std::regex_replace(msolve_name_base, std::regex("\\^"), "");
        msolve_name_base = std::regex_replace(msolve_name_base, std::regex("/"), "");
        msolve_name_base = std::regex_replace(msolve_name_base, std::regex("\\s+"), "");
        if (msolve_name_base.empty() || std::isdigit(msolve_name_base[0])) {
            msolve_name_base = "v" + msolve_name_base; // Ensure starts with a letter
        }


        std::string msolve_name = msolve_name_base;
        int suffix = 1;
        // Check for collisions only with names *already added* to msolve_var_names_ordered
        while (std::find(msolve_var_names_ordered.begin(), msolve_var_names_ordered.end(), msolve_name) !=
               msolve_var_names_ordered.end()) {
            msolve_name = msolve_name_base + "_" + std::to_string(suffix++);
        }

        var_to_msolve_name_map[var] = msolve_name;
        msolve_var_names_ordered.push_back(msolve_name);

        std::cout << "  Mapped '" << var << "' to msolve name '" << msolve_name << "'" << std::endl;
    }

    // First line: variables as comma-separated list, using the ordered list
    for (size_t i = 0; i < msolve_var_names_ordered.size(); ++i) {
        ss_vars << msolve_var_names_ordered[i] << (i == msolve_var_names_ordered.size() - 1 ? "" : ",");
    }
    ss_vars << "\n";

    // Second line: field characteristic (0 for rational)
    ss_vars << "0\n";

    // Following lines: polynomials with comma separators
    for (size_t i = 0; i < system.polynomials.size(); ++i) {
        ss_vars << polynomial_to_msolve_string(system.polynomials[i], var_to_msolve_name_map)
                << (i == system.polynomials.size() - 1 ? "" : ",\n");
    }
    ss_vars << "\n";

    return ss_vars.str();
}

PolynomialSolutionSet
MSolveSolver::solve(const AlgebraicSystem &system) {
    initialize_msolve_if_needed();
    PolynomialSolutionSet solutions;
    std::map<Variable, std::string> var_to_msolve_name_map; // Declare map here

    std::string msolve_input_str = convert_to_msolve_format(system, var_to_msolve_name_map); // Pass by ref
    if (msolve_input_str.empty()) {
        std::cerr << "[MSolveSolver] Error: System string is empty after conversion." << std::endl;
        return {};
    }

    std::cout << "  [MSolveSolver] Generated msolve input format:\n" << msolve_input_str << std::endl;

    std::string temp_msolve_input_file = generate_temp_filename("msolve_in");
    std::string temp_msolve_output_file = generate_temp_filename("msolve_raw_out");

    std::ofstream msolve_in_stream(temp_msolve_input_file);
    if (!msolve_in_stream) {
        std::cerr << "[MSolveSolver] Error: Could not create temporary input file: " << temp_msolve_input_file
                  << std::endl;
        return {};
    }
    msolve_in_stream << msolve_input_str;
    msolve_in_stream.close();

    std::stringstream msolve_cmd_ss;
    const char *msolve_binary_paths[] = {
        "build/vcpkg_installed/x64-linux/tools/msolve/bin/msolve", "msolve", "/usr/bin/msolve", "/usr/local/bin/msolve"
    };
    bool found_msolve_binary = false;
    for (const char *path : msolve_binary_paths) {
        if (access(path, X_OK) == 0) {
            msolve_cmd_ss << path;
            found_msolve_binary = true;
            break;
        }
    }
    if (!found_msolve_binary) {
        std::cerr << "[MSolveSolver] Error: msolve binary not found." << std::endl;
        remove(temp_msolve_input_file.c_str());
        return {};
    }
    msolve_cmd_ss << " -P 2 -v 0 -f " << temp_msolve_input_file << " -o " << temp_msolve_output_file;
    std::string msolve_cmd = msolve_cmd_ss.str();

    std::cout << "  [MSolveSolver] Executing msolve: " << msolve_cmd << std::endl;
    int msolve_exec_result = std::system(msolve_cmd.c_str());

    if (remove(temp_msolve_input_file.c_str()) != 0 && msolve_exec_result == 0) {
        std::cerr << "[MSolveSolver] Warning: Could not remove temp msolve input file: " << temp_msolve_input_file
                  << std::endl;
    }

    if (msolve_exec_result != 0) {
        std::cerr << "[MSolveSolver] Error: msolve command failed with code: " << msolve_exec_result << std::endl;
        remove(temp_msolve_output_file.c_str());
        return {};
    }

    // ---- BEGIN ADDED DEBUG CODE ----
    // Read and print the content of temp_msolve_output_file for debugging
    std::ifstream msolve_out_debug_stream(temp_msolve_output_file);
    if (msolve_out_debug_stream) {
        std::string msolve_raw_output_content((std::istreambuf_iterator<char>(msolve_out_debug_stream)),
                                              std::istreambuf_iterator<char>());
        msolve_out_debug_stream.close();
        std::cout << "  [MSolveSolver] Raw msolve output file content (" << temp_msolve_output_file << "):\n------\n"
                  << msolve_raw_output_content << "\n------" << std::endl;
    } else {
        std::cerr << "  [MSolveSolver] Warning: Could not open raw msolve output file for debugging: "
                  << temp_msolve_output_file << std::endl;
    }
    // ---- END ADDED DEBUG CODE ----

    // Now, call the Python script to parse msolve's output file and get JSON
    std::stringstream python_cmd_ss;
    std::string effective_script_path;

    // 1. Try the configured python_script_path_ first if it's not empty
    if (!python_script_path_.empty()) {
        std::ifstream script_file_check(python_script_path_);
        if (script_file_check.good()) {
            effective_script_path = python_script_path_;
            std::cout << "  [MSolveSolver] Using configured python_script_path_: " << effective_script_path
                      << std::endl;
        }
        script_file_check.close();
    }

    // 2. If not found or not configured, try standard locations
    if (effective_script_path.empty()) {
        std::string path1 = "scripts/msolve_p_to_json.py";
        std::ifstream script_file_check1(path1);
        if (script_file_check1.good()) {
            effective_script_path = path1;
            std::cout << "  [MSolveSolver] Found Python script at: " << effective_script_path << std::endl;
        }
        script_file_check1.close();

        if (effective_script_path.empty()) {
            std::string path2 = "../scripts/msolve_p_to_json.py";
            std::ifstream script_file_check2(path2);
            if (script_file_check2.good()) {
                effective_script_path = path2;
                std::cout << "  [MSolveSolver] Found Python script at: " << effective_script_path << std::endl;
            }
            script_file_check2.close();
        }
    }

    // 3. If still not found, error out
    if (effective_script_path.empty()) {
        std::cerr << "[MSolveSolver] Error: Python script 'msolve_p_to_json.py' not found in standard locations."
                  << std::endl;
        if (!python_script_path_.empty()) {
            std::cerr << "  (Configured python_script_path_ was: '" << python_script_path_
                      << "' but it was not found either)" << std::endl;
        } else {
            std::cerr << "  (python_script_path_ was not configured)" << std::endl;
        }
        remove(temp_msolve_output_file.c_str()); // Clean up msolve output
        return {};
    }

    python_cmd_ss << "python3 " << effective_script_path << " " << temp_msolve_output_file;
    std::string python_cmd = python_cmd_ss.str();
    std::cout << "  [MSolveSolver] Executing Python parser: " << python_cmd << std::endl;

    std::string json_output_str;
    std::array<char, 128> popen_buffer;
    FILE *pipe = popen(python_cmd.c_str(), "r");
    if (!pipe) {
        std::cerr << "[MSolveSolver] Error: popen() failed for Python script!" << std::endl;
        remove(temp_msolve_output_file.c_str());
        return {};
    }
    while (fgets(popen_buffer.data(), popen_buffer.size(), pipe) != nullptr) { json_output_str += popen_buffer.data(); }
    int python_script_status = pclose(pipe);

    if (remove(temp_msolve_output_file.c_str()) != 0) {
        std::cerr << "[MSolveSolver] Warning: Could not remove temp msolve output file: " << temp_msolve_output_file
                  << std::endl;
    }

    std::cout << "  [MSolveSolver] NOTE: Python script processed msolve output file: " << temp_msolve_output_file
              << std::endl;
    std::cout
      << "  [MSolveSolver] To manually debug python script, re-run test and then use this file path with python3."
      << std::endl;

    if (python_script_status != 0) {
        std::cerr << "[MSolveSolver] Error: Python script execution failed or returned error. Exit status: "
                  << python_script_status << std::endl;
        std::cerr << "[MSolveSolver] Python script output (if any):\n" << json_output_str << std::endl;
        return {};
    }

    std::cout << "  [MSolveSolver] JSON output from Python script: \n" << json_output_str << std::endl;

    try {
        if (json_output_str.empty()) {
            std::cerr << "[MSolveSolver] Error: JSON output from Python script is empty." << std::endl;
            return {};
        }
        nlohmann::json parsed_json = nlohmann::json::parse(json_output_str);
        std::cout << "  [MSolveSolver] Successfully parsed JSON. Structure is expected to be rational parametrization."
                  << std::endl;
        std::cout << "  [MSolveSolver] Full parsed JSON for debug: \n"
                  << parsed_json.dump(2) << std::endl; // Print pretty JSON

        // --- Start of new JSON parsing logic for rational parametrization ---

        if (!parsed_json.is_array() || parsed_json.empty()) {
            std::cerr
              << "[MSolveSolver] Error: Parsed JSON is not a non-empty array as expected for rational parametrization."
              << std::endl;
            return {};
        }

        int dim_flag = parsed_json[0].get<int>();
        std::cout << "  [MSolveSolver] Dimension flag from JSON: " << dim_flag << std::endl;

        if (dim_flag != 0) {
            std::cerr << "[MSolveSolver] Error: System dimension is " << dim_flag
                      << ". MSolveSolver currently only supports 0-dimensional systems (finite solutions)."
                      << std::endl;
            // Per PLAN.md, defer handling non-zero-dim cases.
            return {}; // Return empty solution set
        }

        // If dim_flag == 0, proceed to extract parametrization components from parsed_json[1]
        // The structure of parsed_json[1] is [<characteristic>, <num_vars_in_param_rep>, <degree_system>, <vars_list>,
        // <linear_form_coeffs>, <parametrization_data>]

        if (parsed_json.size() < 2 || !parsed_json[1].is_array() || parsed_json[1].size() < 6) {
            std::cerr << "[MSolveSolver] Error: Parsed JSON[1] is not a valid array or is too short for rational "
                         "parametrization data."
                      << std::endl;
            return {};
        }

        const auto &param_components = parsed_json[1];

        // Extract Variable Names: vars_list is at param_components[3]
        if (!param_components[3].is_array()) {
            std::cerr << "[MSolveSolver] Error: Variable list (param_components[3]) is not an array." << std::endl;
            return {};
        }

        std::vector<std::string> stored_poly_vars;
        for (const auto &var_item : param_components[3]) {
            if (!var_item.is_string()) {
                std::cerr << "[MSolveSolver] Error: Variable name in list is not a string." << std::endl;
                return {}; // Or handle more gracefully
            }
            stored_poly_vars.push_back(var_item.get<std::string>());
        }

        std::cout << "  [MSolveSolver] Extracted msolve variable names: ";
        for (size_t i = 0; i < stored_poly_vars.size(); ++i) {
            std::cout << stored_poly_vars[i] << (i == stored_poly_vars.size() - 1 ? "" : ", ");
        }
        std::cout << std::endl;


        // P_ is param_components[5]
        // P_ = [flag_or_count, CoreParametrization]
        // CoreParametrization = [lw_struct, lwp_struct, params_list_struct]
        // lw_struct = [degree, coeffs_list]

        if (!param_components[5].is_array() || param_components[5].size() < 2) {
            std::cerr
              << "[MSolveSolver] Error: Outer parametrization data (param_components[5]) is not an array or too short."
              << std::endl;
            return {};
        }
        const auto &OuterP_data = param_components[5];

        // Assuming OuterP_data[0] is a flag/count, typically 1 for our cases.
        // CoreParametrization is OuterP_data[1]
        if (!OuterP_data[1].is_array() ||
            OuterP_data[1].size() <
              2) { // CoreParametrization must have at least lw, lwp. params_list can be empty/null.
            std::cerr << "[MSolveSolver] Error: Core parametrization data (OuterP_data[1]) is not an array or too "
                         "short for lw, lwp."
                      << std::endl;
            return {};
        }
        const auto &CoreParametrization = OuterP_data[1];

        // Extract Eliminating Polynomial w(t) (lw) from CoreParametrization[0]
        if (!CoreParametrization[0].is_array() || CoreParametrization[0].size() != 2 ||
            !CoreParametrization[0][1].is_array()) {
            std::cerr << "[MSolveSolver] Error: Invalid structure for eliminating polynomial w(t) (lw_struct) in JSON."
                      << std::endl;
            return {};
        }
        const auto &lw_struct = CoreParametrization[0];
        const auto &lw_coeffs_json = lw_struct[1];
        std::vector<std::string> w_coeffs_str;
        for (const auto &coeff_item : lw_coeffs_json) {
            if (coeff_item.is_number()) {
                w_coeffs_str.push_back(std::to_string(coeff_item.get<double>()));
            } else if (coeff_item.is_string()) {
                w_coeffs_str.push_back(coeff_item.get<std::string>());
            } else {
                std::cerr << "[MSolveSolver] Error: Coefficient for w(t) is not a number or string." << std::endl;
                return {};
            }
        }
        std::cout << "  [MSolveSolver] Extracted w(t) coefficients (as strings): ";
        for (size_t i = 0; i < w_coeffs_str.size(); ++i) {
            std::cout << w_coeffs_str[i] << (i == w_coeffs_str.size() - 1 ? "" : ", ");
        }
        std::cout << std::endl;

        // Extract Derivative Polynomial w'(t) (lwp) from CoreParametrization[1]
        if (!CoreParametrization[1].is_array() || CoreParametrization[1].size() != 2 ||
            !CoreParametrization[1][1].is_array()) {
            std::cerr << "[MSolveSolver] Error: Invalid structure for derivative polynomial w'(t) (lwp_struct) in JSON."
                      << std::endl;
            return {};
        }
        const auto &lwp_struct = CoreParametrization[1];
        const auto &lwp_coeffs_json = lwp_struct[1];
        std::vector<std::string> wp_coeffs_str;
        for (const auto &coeff_item : lwp_coeffs_json) {
            if (coeff_item.is_number()) {
                wp_coeffs_str.push_back(std::to_string(coeff_item.get<double>()));
            } else if (coeff_item.is_string()) {
                wp_coeffs_str.push_back(coeff_item.get<std::string>());
            } else {
                std::cerr << "[MSolveSolver] Error: Coefficient for w'(t) is not a number or string." << std::endl;
                return {};
            }
        }
        std::cout << "  [MSolveSolver] Extracted w'(t) coefficients (as strings): ";
        for (size_t i = 0; i < wp_coeffs_str.size(); ++i) {
            std::cout << wp_coeffs_str[i] << (i == wp_coeffs_str.size() - 1 ? "" : ", ");
        }
        std::cout << std::endl;

        // Structure to hold parsed parametrizing polynomial v_i(t) = Numerator_i(t) / Denominator_i
        struct ParametrizingPolynomial {
            std::vector<std::string> numerator_coeffs_str;
            double denominator_val;
        };

        // Extract Parametrizing Polynomials v_i(t) list (params_list_struct) from CoreParametrization[2]
        std::vector<ParametrizingPolynomial> v_polynomials_list;
        const nlohmann::json *params_list_struct_ptr = nullptr;

        if (CoreParametrization.size() >= 3 && CoreParametrization[2].is_array()) {
            params_list_struct_ptr = &CoreParametrization[2];
            // Handle [null] case explicitly for empty list
            if (params_list_struct_ptr->size() == 1 && (*params_list_struct_ptr)[0].is_null()) {
                std::cout << "  [MSolveSolver] Parametrizing polynomials list (params_list_struct) is [null], treating "
                             "as empty."
                          << std::endl;
                // v_polynomials_list will remain empty
            } else {
                v_polynomials_list.reserve(params_list_struct_ptr->size());
                for (size_t i = 0; i < params_list_struct_ptr->size(); ++i) {
                    const auto &v_poly_data_json = (*params_list_struct_ptr)[i];
                    if (!v_poly_data_json.is_array() || v_poly_data_json.size() != 2 ||
                        !v_poly_data_json[0].is_array() || v_poly_data_json[0].size() != 2 ||
                        !v_poly_data_json[0][1].is_array() || !v_poly_data_json[1].is_number()) {
                        std::cerr << "[MSolveSolver] Error: Invalid structure for parametrizing polynomial v_" << i
                                  << "(t) in JSON." << std::endl;
                        return {};
                    }
                    const auto &numerator_coeffs_json = v_poly_data_json[0][1];
                    std::vector<std::string> current_v_num_coeffs_str;
                    for (const auto &coeff_item : numerator_coeffs_json) {
                        if (coeff_item.is_number()) {
                            current_v_num_coeffs_str.push_back(std::to_string(coeff_item.get<double>()));
                        } else if (coeff_item.is_string()) {
                            current_v_num_coeffs_str.push_back(coeff_item.get<std::string>());
                        } else {
                            std::cerr << "[MSolveSolver] Error: Numerator coefficient for v_" << i
                                      << "(t) is not a number or string." << std::endl;
                            return {};
                        }
                    }
                    double denominator = v_poly_data_json[1].get<double>();
                    v_polynomials_list.push_back({ current_v_num_coeffs_str, denominator });
                }
            }
        } else {
            std::cout << "  [MSolveSolver] Parametrizing polynomials list (params_list_struct) is missing or not an "
                         "array, treating as empty."
                      << std::endl;
            // v_polynomials_list will remain empty
        }

        std::cout << "  [MSolveSolver] Processed v_polynomials_list. Found " << v_polynomials_list.size() << " entries."
                  << std::endl;
        for (size_t i = 0; i < v_polynomials_list.size(); ++i) {
            const auto &v_poly_entry = v_polynomials_list[i];
            std::cout << "    v_poly_" << i << ": NumCoeffs=[ ";
            for (size_t j = 0; j < v_poly_entry.numerator_coeffs_str.size(); ++j) {
                std::cout << v_poly_entry.numerator_coeffs_str[j]
                          << (j == v_poly_entry.numerator_coeffs_str.size() - 1 ? "" : ", ");
            }
            std::cout << " ], Denominator=" << v_poly_entry.denominator_val << std::endl;
        }

        // Validate the size of v_polynomials_list
        // It should have N-1 entries if stored_poly_vars has N entries (where the Nth is 't')
        // Or 0 entries if stored_poly_vars has 1 entry (which is 't')
        size_t expected_v_list_size =
          stored_poly_vars.empty() ? 0 : (stored_poly_vars.size() > 0 ? stored_poly_vars.size() - 1 : 0);
        if (v_polynomials_list.size() != expected_v_list_size) {
            std::cerr << "[MSolveSolver] Error: Size of v_polynomials_list (" << v_polynomials_list.size()
                      << ") does not match expected size (" << expected_v_list_size
                      << ", which is stored_poly_vars.size()-1 = " << stored_poly_vars.size() << "-1)." << std::endl;
            // This might indicate an issue with JSON structure or msolve output interpretation.
            return {};
        }

        // TEMP: Remove this clear once solution reconstruction is complete
        // solutions.clear();

        // --- Phase MS-3: Univariate Solving and Solution Reconstruction ---
        std::cout << "  [MSolveSolver] Phase MS-3: Univariate Solving and Solution Reconstruction..." << std::endl;

        // 1. Call univariate_solver::find_roots for w(t)
        // Assuming a default precision for now. This might need to be configurable.
        int flint_precision = 128; // Example precision
        std::vector<std::complex<double>> w_roots;
        try {
            std::cout << "    Calling UnivariateSolverFlint::find_roots for w(t)..." << std::endl;
            w_roots = poly_ode::univariate_solver::find_roots(w_coeffs_str, flint_precision);
            std::cout << "    w(t) solve completed. Found " << w_roots.size() << " root(s)." << std::endl;
        } catch (const std::exception &e) {
            std::cerr << "[MSolveSolver] Error during univariate solve of w(t): " << e.what() << std::endl;
            return {}; // Return empty if w(t) solve fails
        }

        if (w_roots.empty()) {
            std::cout << "  [MSolveSolver] w(t) has no roots. No solutions to reconstruct." << std::endl;
            // solutions is already empty or cleared if we started with it, or will be returned empty.
            return solutions; // Should be an empty set
        }

        // 2. For each root theta_j of w(t):
        for (size_t root_idx = 0; root_idx < w_roots.size(); ++root_idx) {
            const auto &theta_j = w_roots[root_idx];
            std::cout << "    Processing root theta_" << root_idx << " = " << theta_j << std::endl;

            // Evaluate w'(theta_j)
            std::complex<double> w_prime_at_theta_j;
            try {
                w_prime_at_theta_j = MSolveSolver::evaluate_poly_at_complex(wp_coeffs_str, theta_j);
            } catch (const std::exception &e) {
                std::cerr << "[MSolveSolver] Error evaluating w'(theta_" << root_idx << "): " << e.what() << std::endl;
                continue; // Skip this root if w' evaluation fails
            }
            std::cout << "      w'(theta_" << root_idx << ") = " << w_prime_at_theta_j << std::endl;

            // Handle case where w'(theta_j) is zero (or very close to zero)
            double zero_tolerance = 1e-12; // Tolerance for checking against zero
            if (std::abs(w_prime_at_theta_j.real()) < zero_tolerance &&
                std::abs(w_prime_at_theta_j.imag()) < zero_tolerance) {
                std::cerr << "[MSolveSolver] Warning: w'(theta_" << root_idx
                          << ") is close to zero. Skipping this root to avoid division by zero." << std::endl;
                continue;
            }

            PolynomialSolutionMap current_solution_map;
            bool reconstruction_ok_for_this_root = true;

            // Determine which variable in stored_poly_vars is 't'
            // As per msolve-tutorial, the last variable in 'vars' (our stored_poly_vars) is 't'
            std::string t_var_name_msolve;
            if (!stored_poly_vars.empty()) {
                t_var_name_msolve = stored_poly_vars.back();
            } else {
                std::cerr << "[MSolveSolver] Error: stored_poly_vars is empty, cannot determine 't' for reconstruction."
                          << std::endl;
                reconstruction_ok_for_this_root = false; // Should not happen if w_roots is not empty
            }

            std::map<std::string, std::complex<double>> msolve_var_values; // Temp store for values of msolve's vars

            // Assign theta_j to msolve's 't' variable
            if (reconstruction_ok_for_this_root) {
                msolve_var_values[t_var_name_msolve] = theta_j;
                std::cout << "        Assigned msolve var '" << t_var_name_msolve << "' (t) = " << theta_j << std::endl;
            }

            // Reconstruct the values for the other original variables using v_i(t) = Num_i(t) / Den_i
            // The formula from msolve-tutorial p.14, for x_k(theta_j) = v_k(theta_j) / w'(theta_j)
            // where v_k(t) = Numerator_k(t) / Denominator_k
            // So, x_k(theta_j) = (Numerator_k(theta_j) / Denominator_k) / w'(theta_j)

            // num_orig_vars_for_params was for dim_flag == 1, for dim_flag == 0, we iterate over v_polynomials_list
            for (size_t var_idx = 0; var_idx < v_polynomials_list.size(); ++var_idx) { // Corrected loop limit
                const std::string &msolve_var_s =
                  stored_poly_vars[var_idx]; // Assumes v_polynomials align with first N-1 stored_poly_vars
                const auto &v_poly_entry = v_polynomials_list[var_idx];

                // Numerator_i(theta_j) / ( Denominator_i * w'(theta_j) )
                std::complex<double> v_num_at_theta_j;
                try {
                    v_num_at_theta_j =
                      MSolveSolver::evaluate_poly_at_complex(v_poly_entry.numerator_coeffs_str, theta_j);
                } catch (const std::exception &e) {
                    std::cerr << "[MSolveSolver] Error evaluating numerator of v_" << var_idx << "(theta_" << root_idx
                              << "): " << e.what() << std::endl;
                    reconstruction_ok_for_this_root = false;
                    break;
                }

                if (std::abs(v_poly_entry.denominator_val) < 1e-12) {
                    std::cerr << "[MSolveSolver] Error: Denominator for v_" << var_idx << " is zero." << std::endl;
                    reconstruction_ok_for_this_root = false;
                    break;
                }

                // Detailed calculation for debugging
                // HYPOTHESIS: Negate the v_k(t) term based on empirical results with msinp3.txt
                std::complex<double> num_div_den = -(v_num_at_theta_j / v_poly_entry.denominator_val); // ADDED NEGATION
                std::complex<double> final_calc_val = num_div_den / w_prime_at_theta_j;

                std::cout << "        DEBUG_RECON_DETAIL for '" << msolve_var_s << "':" << std::endl;
                std::cout << "          v_poly.num_coeffs_str: [";
                for (size_t k = 0; k < v_poly_entry.numerator_coeffs_str.size(); ++k) {
                    std::cout << "'" << v_poly_entry.numerator_coeffs_str[k] << "'"
                              << (k < v_poly_entry.numerator_coeffs_str.size() - 1 ? ", " : "");
                }
                std::cout << "]" << std::endl;
                std::cout << "          theta_j: " << theta_j << std::endl;
                std::cout << "          v_num_evaluated (v_num_at_theta_j): " << v_num_at_theta_j << std::endl;
                std::cout << "          v_poly.den_val (denominator_val): " << v_poly_entry.denominator_val
                          << std::endl;
                std::cout << "          num_div_den (-(v_num_at_theta_j / denominator_val)): " << num_div_den
                          << std::endl; // Note negation in comment
                std::cout << "          w_prime_at_theta_j: " << w_prime_at_theta_j << std::endl;
                std::cout << "          FINAL_CALCULATED_VALUE (num_div_den / w_prime_at_theta_j): " << final_calc_val
                          << std::endl;

                msolve_var_values[msolve_var_s] = final_calc_val;
                std::cout << "        Reconstructed msolve var '" << msolve_var_s << "' = " << final_calc_val
                          << std::endl;
            }

            // Map reconstructed msolve_var_values to the original system.unknowns
            if (reconstruction_ok_for_this_root) {
                for (const auto &unknown_var : system.unknowns) {
                    // Find the msolve name corresponding to this system.unknown (Variable object)
                    auto map_it = var_to_msolve_name_map.find(unknown_var);
                    if (map_it == var_to_msolve_name_map.end()) {
                        std::cerr << "[MSolveSolver] Error: Could not find system unknown '" << unknown_var.name
                                  << "' in var_to_msolve_name_map." << std::endl;
                        reconstruction_ok_for_this_root = false;
                        break;
                    }
                    const std::string &msolve_name_for_system_var = map_it->second;

                    // Now find this msolve_name in our reconstructed values
                    auto value_it = msolve_var_values.find(msolve_name_for_system_var);
                    if (value_it != msolve_var_values.end()) {
                        current_solution_map[unknown_var] = value_it->second;
                    } else {
                        // This should ideally not happen if msolve_var_values contains all msolve names
                        // including the one that became 't'.
                        std::cerr << "[MSolveSolver] Error: Could not find reconstructed value for msolve var '"
                                  << msolve_name_for_system_var << "' (derived from system unknown '"
                                  << unknown_var.name << "')." << std::endl;
                        reconstruction_ok_for_this_root = false;
                        break;
                    }
                }
            }

            if (reconstruction_ok_for_this_root && current_solution_map.size() == system.unknowns.size()) {
                solutions.push_back(current_solution_map);
            }
        }

        /* OLD JSON parsing logic for real solutions - TO BE REMOVED/REPLACED
        std::string status = parsed_json.value("status", "error");

        if (status == "success") {
            // ... old logic for iterating parsed_json["solutions"] ...
        } else {
            std::string msg = parsed_json.value("message", "No additional message.");
            std::cout << "  [MSolveSolver] Status from Python script: " << status << " - Message: " << msg << std::endl;
        }
        */

    } catch (const nlohmann::json::parse_error &e) {
        std::cerr << "[MSolveSolver] Error: Failed to parse JSON output from Python script: " << e.what() << std::endl;
        std::cerr << "[MSolveSolver] Raw JSON string was: " << json_output_str << std::endl;
        return {};
    }

    if (solutions.empty() && system.unknowns.size() == 1 && system.polynomials.size() == 1) {
        std::cout << "  [MSolveSolver] INFO: No solutions from msolve; attempting local univariate solver fallback."
                  << std::endl;

        const auto &poly = system.polynomials[0];
        const auto &var = system.unknowns[0];

        // Determine polynomial degree and build coefficient vector (ascending powers)
        int max_degree = 0;
        for (const auto &mono : poly.monomials) {
            int exp = 0;
            if (!mono.vars.empty()) { exp = mono.vars.begin()->second; }
            max_degree = std::max(max_degree, exp);
        }

        Eigen::VectorXd coeffs(max_degree + 1);
        coeffs.setZero();
        for (const auto &mono : poly.monomials) {
            int exp = 0;
            if (!mono.vars.empty()) { exp = mono.vars.begin()->second; }
            coeffs[exp] += mono.coeff;
        }

        // Ensure leading coefficient (highest power) is non-zero; otherwise, reduce degree.
        int actual_degree = max_degree;
        while (actual_degree > 0 && std::abs(coeffs[actual_degree]) < 1e-12) { --actual_degree; }
        coeffs.conservativeResize(actual_degree + 1);

        if (actual_degree == 0) {
            std::cerr << "  [MSolveSolver] Fallback aborted: polynomial is constant." << std::endl;
        } else {
            Eigen::PolynomialSolver<double, Eigen::Dynamic> psolver;
            psolver.compute(coeffs);
            auto eigen_roots = psolver.roots();

            std::ostringstream oss_key;
            oss_key << var;
            std::string var_key = oss_key.str();

            for (int i = 0; i < eigen_roots.size(); ++i) {
                std::complex<double> root = eigen_roots[i];
                PolynomialSolutionMap sol_map;
                sol_map[var_key] = root;
                solutions.push_back(sol_map);
            }
            std::cout << "  [MSolveSolver] Fallback solver produced " << eigen_roots.size() << " root(s)." << std::endl;
        }
    }

    std::cout << "  [MSolveSolver] Solve finished. Found " << solutions.size() << " solution(s) after JSON processing."
              << std::endl;

    return solutions;
}

// Helper function to convert string coefficient to double
// Supports integers, decimals, and fractions like "num/den"
double
MSolveSolver::string_coeff_to_double(const std::string &s_coeff) {
    if (s_coeff.find('/') != std::string::npos) {
        std::string num_str = s_coeff.substr(0, s_coeff.find('/'));
        std::string den_str = s_coeff.substr(s_coeff.find('/') + 1);
        try {
            double num = std::stod(num_str);
            double den = std::stod(den_str);
            if (den == 0) { throw std::runtime_error("Division by zero in fractional coefficient: " + s_coeff); }
            return num / den;
        } catch (const std::exception &e) {
            throw std::runtime_error("Invalid fractional coefficient format: " + s_coeff + " Error: " + e.what());
        }
    } else {
        try {
            double val_before_return = std::stod(s_coeff);
            std::cout << "    DEBUG_STOD: input_str='" << s_coeff << "' output_double=" << val_before_return
                      << std::endl;
            return val_before_return;
        } catch (const std::exception &e) {
            throw std::runtime_error("Invalid non-fractional coefficient format: " + s_coeff + ", error: " + e.what());
        }
    }
}

// Helper function to evaluate a polynomial (coeffs given as strings) at a complex point t
// Polynomial is c0 + c1*t + c2*t^2 + ...
std::complex<double>
MSolveSolver::evaluate_poly_at_complex(const std::vector<std::string> &coeffs_str, std::complex<double> t) {
    std::complex<double> result = { 0.0, 0.0 };
    std::complex<double> t_power = { 1.0, 0.0 }; // Starts at t^0

    for (const std::string &s_coeff : coeffs_str) {
        double coeff_val = MSolveSolver::string_coeff_to_double(s_coeff); // Call static method
        result += coeff_val * t_power;
        t_power *= t; // Next power of t
    }
    return result;
}

} // namespace poly_ode