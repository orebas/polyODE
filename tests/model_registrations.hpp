#ifndef MODEL_REGISTRATIONS_HPP
#define MODEL_REGISTRATIONS_HPP

#include "model_test_framework.hpp"

namespace poly_ode {
namespace test_framework {

/**
 * @brief Register all available test models with the global registry
 * 
 * This function should be called once to populate the global model registry
 * with all available test models organized by category.
 */
void register_all_models();

/**
 * @brief Register models from specific categories
 */
void register_identifiability_models();
void register_simple_models();
void register_classical_models();
void register_biological_models();
void register_advanced_models();

/**
 * @brief Individual model registration functions
 */
namespace models {

// Identifiability test models
void register_trivial_unident();
void register_sum_test();
void register_global_unident_test();
void register_substr_test();

// Simple models  
void register_simple();
void register_onesp_cubed();
void register_threesp_cubed();
void register_simple_linear_combination();

// Classical models
void register_lotka_volterra();
void register_harmonic();
void register_vanderpol();
void register_brusselator();

// Biological models
void register_seir();
void register_treatment();
void register_hiv();

// Advanced models
void register_daisy_ex3();
void register_daisy_mamil3();
void register_fitzhugh_nagumo();
void register_crauste();

} // namespace models

} // namespace test_framework
} // namespace poly_ode

#endif // MODEL_REGISTRATIONS_HPP