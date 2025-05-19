#include "phc_solver.hpp"
#include "polynomial.hpp" // For Variable, Polynomial, Monomial
#include <nlohmann/json.hpp>

#include <array> // For popen read buffer
#include <chrono>
#include <complex>
#include <cstdio>  // For std::remove, popen, pclose
#include <cstdlib> // For std::system (though preferring popen)
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/types.h>
#include <unistd.h> // For getpid()
#include <vector>

namespace poly_ode {

// Constructor
PHCSolver::PHCSolver(std::string phc_executable_path, std::string python_script_path)
  : phc_executable_path_(std::move(phc_executable_path))
  , python_script_path_(std::move(python_script_path)) {}

// Get solver name
std::string
PHCSolver::name() const {
    return "PHCSolver";
}

// Main solve method
PolynomialSolutionSet
PHCSolver::solve(const AlgebraicSystem &system) {
    std::cout << "  [PHCSolver] Starting solve..." << std::endl;
    // --- Check if system is square --- //
    size_t n_unknowns = system.unknowns.size();
    size_t n_polys = system.polynomials.size();
    std::cout << "    [PHCSolver] System check: " << n_unknowns << " unknowns, " << n_polys << " polynomials."
              << std::endl;
    if (n_unknowns != n_polys) {
        std::cerr << "    [PHCSolver] WARNING: System is not square! PHC might behave unexpectedly or hang."
                  << std::endl;
        // Depending on requirements, could throw an error here instead.
    }
    // --- End Check --- //

    // --- DEBUG: Print Input System --- //
    std::cout << "    [PHCSolver] Input System:" << std::endl;
    std::cout << "      Unknowns (" << system.unknowns.size() << "): ";
    for (const auto &uk : system.unknowns) { std::cout << uk << " "; }
    std::cout << std::endl;
    std::cout << "      Polynomials (" << system.polynomials.size() << "):" << std::endl;
    for (size_t i = 0; i < system.polynomials.size(); ++i) {
        std::cout << "        P" << i << ": " << system.polynomials[i] << " = 0" << std::endl;
    }
    // --- END DEBUG ---

    // 1. Create temporary file names
    auto [input_file_path, output_file_base] = create_temp_files();
    std::cout << "    [PHCSolver] Temp input file: " << input_file_path << std::endl;
    std::cout << "    [PHCSolver] Temp output base: " << output_file_base << std::endl;
    std::string const output_file_phc = output_file_base + ".phc_out";
    std::string const output_file_dict = output_file_base + ".dict";
    std::string const output_file_json = output_file_base + ".json";

    // Ensure cleanup using RAII or try-finally block
    // Using a simple struct for RAII cleanup
    struct FileCleaner {
        std::vector<std::string> files_to_delete;
        ~FileCleaner() {
            for (const auto &f : files_to_delete) { std::remove(f.c_str()); }
        }
    };
    FileCleaner cleaner;
    cleaner.files_to_delete.push_back(input_file_path);
    cleaner.files_to_delete.push_back(output_file_phc);
    cleaner.files_to_delete.push_back(output_file_dict);
    cleaner.files_to_delete.push_back(output_file_json);

    // 2. Convert system to PHC format and write to input file
    std::map<Variable, size_t> variable_map; // Map poly_ode::Variable to index for PHC
    std::string phc_input_content = convert_to_phc_format(system, variable_map);
    // --- DEBUG: Print PHC Input --- //
    std::cout << "    [PHCSolver] Generated PHC Input Format:\n------\n" << phc_input_content << "------" << std::endl;
    // --- END DEBUG ---

    {
        std::ofstream input_file(input_file_path);
        if (!input_file) { throw std::runtime_error("Failed to open temporary input file: " + input_file_path); }
        input_file << phc_input_content;
    } // RAII closes file

    // 3. Run PHC blackbox solver: phc -b <input> <output>
    std::string cmd_solve = phc_executable_path_ + " -b \"" + input_file_path + "\" \"" + output_file_phc + "\"";
    std::cout << "    [PHCSolver] Running command: " << cmd_solve << std::endl;
    std::string solve_stdout, solve_stderr;
    int ret_solve = run_command(cmd_solve, solve_stdout, solve_stderr);
    // --- DEBUG: Print Command Output --- //
    std::cout << "      [PHCSolver] phc -b exit code: " << ret_solve << std::endl;
    if (!solve_stdout.empty())
        std::cout << "      [PHCSolver] phc -b stdout/stderr:\n------\n" << solve_stdout << "------" << std::endl;
    // if(!solve_stderr.empty()) std::cout << "      [PHCSolver] phc -b stderr:\n------\n" << solve_stderr << "------"
    // << std::endl;
    // --- END DEBUG ---
    if (ret_solve != 0) {
        throw std::runtime_error("PHC solve command ('" + cmd_solve + "') failed with code " +
                                 std::to_string(ret_solve) + ".\nStderr: " + solve_stderr +
                                 "\nStdout: " + solve_stdout);
    }

    // Check if PHC produced output (basic check)
    if (!std::filesystem::exists(output_file_phc) || std::filesystem::file_size(output_file_phc) == 0) {
        // Sometimes PHC exits 0 but creates an empty file if the system is trivial or inconsistent
        std::cerr << "Warning: PHC solve command finished but output file '" << output_file_phc
                  << "' is empty or missing. Assuming no solutions." << std::endl;
        // Check stderr for clues
        if (!solve_stderr.empty()) { std::cerr << "PHC Stderr was:\n" << solve_stderr << std::endl; }
        if (!solve_stdout.empty()) { std::cerr << "PHC Stdout was:\n" << solve_stdout << std::endl; }
        return {}; // Return empty solution set
    }

    // 4. Run PHC dictionary converter: phc -x <phc_output> <dict_output>
    std::string cmd_dict = phc_executable_path_ + " -x \"" + output_file_phc + "\" \"" + output_file_dict + "\"";
    std::cout << "    [PHCSolver] Running command: " << cmd_dict << std::endl;
    std::string dict_stdout, dict_stderr;
    int ret_dict = run_command(cmd_dict, dict_stdout, dict_stderr);
    // --- DEBUG: Print Command Output --- //
    std::cout << "      [PHCSolver] phc -x exit code: " << ret_dict << std::endl;
    if (!dict_stdout.empty())
        std::cout << "      [PHCSolver] phc -x stdout/stderr:\n------\n" << dict_stdout << "------" << std::endl;
    // if(!dict_stderr.empty()) std::cout << "      [PHCSolver] phc -x stderr:\n------\n" << dict_stderr << "------" <<
    // std::endl;
    // --- END DEBUG ---
    if (ret_dict != 0) {
        throw std::runtime_error("PHC dict command ('" + cmd_dict + "') failed with code " + std::to_string(ret_dict) +
                                 ".\nStderr: " + dict_stderr + "\nStdout: " + dict_stdout);
    }

    // 5. Read dict file content
    std::string dict_content;
    {
        std::ifstream dict_file(output_file_dict);
        if (!dict_file) { throw std::runtime_error("Failed to open temporary dict file: " + output_file_dict); }
        std::stringstream buffer;
        buffer << dict_file.rdbuf();
        dict_content = buffer.str();
    } // RAII closes file
    // --- DEBUG: Print Dict Content --- //
    std::cout << "    [PHCSolver] Dictionary File Content:\n------\n" << dict_content << "------" << std::endl;
    // --- END DEBUG ---
    if (dict_content.empty()) {
        throw std::runtime_error("PHC dict file '" + output_file_dict + "' was created but is empty.");
    }

    // ---- START MODIFIED DEBUG LOG FOR DICT CONTENT ----
    std::cout << "    [PHCSolver] Raw Dictionary File Content (first ~2000 chars):\n------\n";
    std::cout << dict_content.substr(0, 2000) << (dict_content.length() > 2000 ? "..." : "") << "\n------" << std::endl;
    // ---- END MODIFIED DEBUG LOG ----

    // 6. Run Python script to convert dict string to JSON using std::system
    std::string json_error_file = output_file_json + ".err";
    std::string cmd_json_system = "cat \"" + output_file_dict + "\" | python3 \"" + python_script_path_ + "\" > \"" +
                                  output_file_json + "\" 2> \"" + json_error_file + "\"";
    std::cout << "    [PHCSolver] Running command: " << cmd_json_system << std::endl;
    int ret_json_system = std::system(cmd_json_system.c_str());
    std::cout << "      [PHCSolver] Python script exit code: " << WEXITSTATUS(ret_json_system)
              << " (raw: " << ret_json_system << ")" << std::endl;

    // Check exit code
    if (ret_json_system != 0) {
        // Try to read the stderr file for clues
        std::string python_stderr_content;
        std::ifstream stderr_file(json_error_file);
        if (stderr_file) {
            std::stringstream err_buf;
            err_buf << stderr_file.rdbuf();
            python_stderr_content = err_buf.str();
            cleaner.files_to_delete.push_back(json_error_file); // Add stderr file to cleanup
        } else {
            python_stderr_content = "[Could not read stderr file: " + json_error_file + "]";
        }
        throw std::runtime_error("Python JSON conversion script ('" + cmd_json_system + "') failed with exit code " +
                                 std::to_string(WEXITSTATUS(ret_json_system)) +
                                 " (raw: " + std::to_string(ret_json_system) +
                                 ")."
                                 "\nScript Stderr:\n" +
                                 python_stderr_content);
    }

    // Check if JSON file was created and is not empty
    if (!std::filesystem::exists(output_file_json) || std::filesystem::file_size(output_file_json) == 0) {
        throw std::runtime_error("Python JSON conversion script executed successfully, but output file '" +
                                 output_file_json + "' is missing or empty.");
    }

    // Clean up the (likely empty) error file if it exists
    if (std::filesystem::exists(json_error_file)) {
        std::remove(json_error_file.c_str());
        // Optional: Remove from cleaner if already added, but harmless if left
    }

    // 7. Parse JSON output
    std::cout << "    [PHCSolver] Parsing JSON file: " << output_file_json << std::endl;
    PolynomialSolutionSet solutions = parse_phc_json_output(output_file_json, system.unknowns);
    std::cout << "    [PHCSolver] Parsed " << solutions.size() << " solutions." << std::endl;
    std::cout << "  [PHCSolver] Solve finished." << std::endl;
    return solutions;
}

// Convert system to PHC format string
std::string
PHCSolver::convert_to_phc_format(const AlgebraicSystem &system, std::map<Variable, size_t> &variable_map) const {
    std::stringstream ss;
    variable_map.clear();

    // Assign indices 0, 1, 2... to variables based on order in unknowns
    for (size_t i = 0; i < system.unknowns.size(); ++i) { variable_map[system.unknowns[i]] = i; }

    // Write number of equations and variables
    ss << system.polynomials.size() << " " << system.unknowns.size() << "\n";

    // Write each polynomial equation
    for (const auto &poly : system.polynomials) {
        if (poly.monomials.empty()) {
            ss << "0;\n"; // Equation is just 0 = 0
            continue;
        }

        bool first_mono = true;
        for (const auto &mono : poly.monomials) {
            double coeff = mono.coeff;
            if (std::abs(coeff) < 1e-15) { // Tolerance for zero coefficient
                continue;
            }

            // Sign handling
            if (!first_mono) {
                ss << (coeff > 0.0 ? " + " : " - ");
            } else if (coeff < 0.0) {
                ss << "-";
            }

            double abs_coeff = std::abs(coeff);
            bool has_vars = !mono.vars.empty();

            // Coefficient (only if not 1.0, unless it's a constant term)
            if (std::abs(abs_coeff - 1.0) > 1e-15 || !has_vars) {
                ss << abs_coeff;
                if (has_vars) ss << "*";
            }

            // Variables and powers
            bool first_var_in_mono = true;
            for (const auto &var_pair : mono.vars) {
                const Variable &var = var_pair.first;
                int exponent = var_pair.second;

                auto it = variable_map.find(var);
                if (it == variable_map.end()) {
                    // This should not happen if AlgebraicSystem is well-formed
                    throw std::runtime_error("Variable '" + var.name +
                                             "' found in polynomial but not in unknowns list.");
                }
                size_t var_index = it->second;

                if (!first_var_in_mono) { ss << "*"; }
                ss << "x" << var_index; // PHC uses x0, x1, ...
                if (exponent != 1) { ss << "^" << exponent; }
                first_var_in_mono = false;
            }
            first_mono = false;
        }
        ss << ";\n"; // End of equation marker
    }
    return ss.str();
}

// Parse JSON output
PolynomialSolutionSet
PHCSolver::parse_phc_json_output(const std::string &json_file_path,
                                 const std::vector<Variable> &original_unknowns) const {
    std::ifstream json_file(json_file_path);
    if (!json_file.is_open()) { throw std::runtime_error("Could not open JSON file for parsing: " + json_file_path); }

    nlohmann::json json_data;
    try {
        json_data = nlohmann::json::parse(json_file);
    } catch (const nlohmann::json::parse_error &e) {
        throw std::runtime_error("JSON parsing error in file '" + json_file_path + "': " + std::string(e.what()));
    }

    PolynomialSolutionSet solutions;

    if (!json_data.is_array()) {
        throw std::runtime_error("Expected JSON root to be an array of solutions in file: " + json_file_path);
    }

    // Process each solution object in the JSON array
    for (const auto &solution_obj : json_data) {
        if (!solution_obj.is_object()) continue; // Skip non-object entries

        PolynomialSolutionMap current_solution;
        bool solution_valid = true;

        for (auto it = solution_obj.begin(); it != solution_obj.end(); ++it) {
            const std::string &key = it.key();
            // Skip metadata like time, multiplicity, err, rco, res
            if (key == "time" || key == "multiplicity" || key == "err" || key == "rco" || key == "res") { continue; }

            // Expect variable keys like "x0", "x1", ...
            if (key.length() > 1 && key[0] == 'x') {
                try {
                    size_t var_index = std::stoul(key.substr(1));
                    if (var_index >= original_unknowns.size()) {
                        std::cerr << "Warning: Parsed variable index " << var_index
                                  << " out of bounds for original unknowns (size " << original_unknowns.size()
                                  << "). Skipping variable." << std::endl;
                        continue;
                    }

                    const Variable &original_var = original_unknowns[var_index];

                    // Expect value to be an object { "re": ..., "im": ... }
                    if (it.value().is_object() && it.value().contains("re") && it.value().contains("im")) {
                        double real_part = it.value()["re"].get<double>();
                        double imag_part = it.value()["im"].get<double>();
                        current_solution[original_var] = std::complex<double>(real_part, imag_part);
                    } else {
                        std::cerr << "Warning: Expected complex object for variable '" << key
                                  << "' but got different type. Skipping variable." << std::endl;
                        solution_valid = false;
                        break;
                    }
                } catch (const std::invalid_argument &) {
                    std::cerr << "Warning: Could not parse index from variable key '" << key << "'. Skipping variable."
                              << std::endl;
                } catch (const std::out_of_range &) {
                    std::cerr << "Warning: Index parsed from variable key '" << key
                              << "' is out of range for size_t. Skipping variable." << std::endl;
                }
            } else {
                std::cerr << "Warning: Unexpected key '" << key << "' in solution object. Skipping." << std::endl;
            }
        }

        // Add the solution if it was parsed correctly and is not empty
        if (solution_valid && !current_solution.empty()) {
            // TODO: Handle multiplicity if needed (PHC JSON format puts it in the object)
            int multiplicity = 1;
            if (solution_obj.contains("multiplicity")) {
                try {
                    multiplicity = solution_obj["multiplicity"].get<int>();
                } catch (...) { /* ignore error, default to 1 */
                }
            }
            for (int i = 0; i < multiplicity; ++i) { solutions.push_back(current_solution); }
        }
    }

    return solutions;
}

// Run external command helper function (using popen for stdout/stderr capture, no stdin)
int
PHCSolver::run_command(const std::string &command,
                       /* const std::string &std_input, // Removed */
                       std::string &std_output,
                       std::string &std_error) const {
    std_output.clear();
    std_error.clear(); // std_error is primarily for popen/pclose errors now

    std::string command_with_stderr = command + " 2>&1";  // Redirect stderr to stdout
    FILE *pipe = popen(command_with_stderr.c_str(), "r"); // Open for reading combined stdout/stderr
    if (!pipe) {
        std_error = "popen() failed for command: " + command;
        perror(std_error.c_str()); // Print system error
        return -1;                 // Indicate popen failure
    }

    std::array<char, 256> buffer;
    while (fgets(buffer.data(), buffer.size(), pipe) != nullptr) {
        std_output += buffer.data(); // Capture combined output
    }

    int exit_status = pclose(pipe);

    if (exit_status == -1) {
        // pclose failed
        std_error = "pclose() failed for command: " + command;
        perror(std_error.c_str());
        return -1; // Indicate pclose failure
    }

    // std_output contains the combined output. No reliable way to separate stderr here.

    return WEXITSTATUS(exit_status); // Return the command's actual exit code
}

// Create temporary files helper
std::pair<std::string, std::string>
PHCSolver::create_temp_files() const {
    std::filesystem::path const temp_dir = std::filesystem::temp_directory_path();
    auto now = std::chrono::system_clock::now();
    auto timestamp = std::chrono::duration_cast<std::chrono::nanoseconds>(now.time_since_epoch()).count();
    pid_t const pid = getpid();

    // Use thread-local random generator for uniqueness across threads if needed later
    static thread_local std::mt19937 generator(std::random_device{}() + pid);
    std::uniform_int_distribution<int> distribution(0, 99999);
    int random_num = distribution(generator);

    std::string unique_id = std::to_string(timestamp) + "_" + std::to_string(pid) + "_" + std::to_string(random_num);
    std::string input_filename = "polyode_phc_input_" + unique_id + ".txt";
    std::string output_base = "polyode_phc_output_" + unique_id;

    std::filesystem::path input_path = temp_dir / input_filename;
    std::filesystem::path output_path_base = temp_dir / output_base;

    // Return full path for input, base name for output (extensions added later)
    return { input_path.string(), output_path_base.string() };
}

} // namespace poly_ode