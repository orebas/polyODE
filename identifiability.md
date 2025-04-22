Algorithm: Multipoint Local Identifiability Analysis via Sensitivity Matrix Rank

1. Introduction and Motivation

Parameter estimation is a critical step in constructing and validating Ordinary Differential Equation (ODE) models of biological systems. However, due to model structure or limitations in experimental data, not all model parameters (kinetic constants, initial conditions) may be uniquely determinable from the available measurements. This property is known as identifiability. Assessing parameter identifiability prior to parameter estimation is crucial to avoid futile computational effort, guide model refinement, inform experimental design, and correctly interpret estimation results. This algorithm assesses local identifiability by evaluating the rank of the system's sensitivity matrix at multiple points in the parameter space. Local non-identifiability at multiple points is a strong indicator of practical difficulties in achieving unique parameter estimates.

2. Algorithm Overview

The core idea is to numerically evaluate the local sensitivity of the model outputs (including their time derivatives) with respect to the parameters (including initial conditions designated for estimation). This sensitivity information is captured in a Jacobian matrix. If this matrix has full column rank, all tested parameters independently influence the outputs locally, suggesting local identifiability. If the matrix is rank-deficient, the nullspace reveals linear dependencies between parameter sensitivities, indicating local non-identifiability. The algorithm iteratively identifies and removes non-identifiable parameters by analyzing the nullspace, repeating the sensitivity calculation until a full-rank matrix (for the remaining parameters) is obtained. The analysis is performed at multiple, randomly selected points (differing in initial conditions) to increase robustness against coincidental rank deficiencies at a single point. Finally, the algorithm determines the minimum number of output time derivatives required to achieve this maximal rank.

3. Definitions and Setup

    ODE Model: The biological system is described by a set of ODEs: dx/dt = f(x(t), θ_params, t) where x(t) is the vector of state variables (e.g., species concentrations) at time t, and θ_params is the vector of time-invariant model parameters (e.g., kinetic constants).
    Initial Conditions: The state vector at t=0 is x(0) = x_0. Some or all elements of x_0 may be considered unknown parameters to be identified.
    Parameter Vector: Define the full vector of parameters to be assessed for identifiability as θ. This vector typically includes the model parameters θ_params and the subset of initial conditions x_0 that are considered unknown: θ = [θ_params; x_0].
    Measurement Equations: The experimentally measured quantities y_k are functions of the states and parameters: y_k(t) = g_k(x(t), θ_params, t) for k = 1, ..., N_obs.
    Output Vector and Derivatives: Define an extended output vector Y(t, θ) that includes the measured quantities and their time derivatives up to a specified order N: Y(t, θ) = [y_1, dy_1/dt, ..., d^N y_1/dt^N, y_2, ..., d^N y_{N_obs}/dt^N]^T.
    Sensitivity Matrix: The core object of analysis is the sensitivity matrix S, representing the partial derivatives of the extended output vector Y with respect to the parameter vector θ: S(t, θ) = ∂Y(t, θ) / ∂θ^T.

4. Detailed Algorithm Steps

    Input:
        ODE model functions f and g.
        The initial parameter vector θ containing all parameters and initial conditions subject to identifiability analysis.
        Number of random test points M (e.g., M=10-100). Points differ in their initial condition values.
        Maximum derivative order N for outputs (e.g., N = ceil(length(θ) / N_obs) + 2).
        Numerical tolerances for rank determination (ε_rank) and nullspace vector magnitude (ε_nullspace).

    Step 1: Initialization
        Maintain a list of parameters currently considered identifiable, θ_identifiable, initialized to θ.
        Maintain a dictionary of parameters found to be non-identifiable, θ_fixed, initialized as empty.
        Generate one set of random numerical values for the kinetic parameters θ_params.
        Generate M sets of random numerical values for the initial conditions x_0^(i), i=1..M.
        Define the M test points P_i = (θ_params, x_0^(i)). These represent different initial states for the same underlying kinetic model.

    Step 2: Pre-computation of Symbolic Derivatives (Optional but Recommended)
        To efficiently compute sensitivities, it is beneficial to pre-compute symbolic expressions for the time derivatives d^n y_k / dt^n up to order N.
        This involves repeatedly differentiating y_k = g_k(x, θ_params) with respect to time t, using the chain rule and substituting dx/dt = f(x, θ_params) whenever dx/dt terms appear. This yields expressions for d^n y_k / dt^n in terms of x, θ_params, and potentially lower-order time derivatives of x.
        Note: If symbolic pre-computation is infeasible, sensitivities involving higher derivatives can be computed using nested automatic differentiation or finite differences during Step 4, albeit potentially less efficiently. Any parameters already in θ_fixed should be substituted as constants during this stage.

    Step 2.1: Refined Strategy for Sensitivity Computation (Avoiding Symbolic Explosion)
        Direct symbolic substitution of dx/dt into higher derivatives of y_k can lead to extremely large expressions (expression swell), making computation infeasible for complex models or higher derivative orders.
        A more computationally tractable approach involves:
        a. Symbolic Differentiation: Compute symbolic derivatives d^n(x_i)/dt^n and d^n(y_k)/dt^n up to order N using `differentiate_wrt_t`, storing each derivative as a distinct `RationalFunction`. Do *not* perform symbolic forward substitutions (e.g., replacing dx/dt terms within d^2y/dt^2).
        b. Numerical Evaluation Function: Implement a function `compute_Y_numerical(θ_values)` that takes numerical values for the parameters/ICs (`θ_values`) and calculates the numerical value of the extended output vector Y. This function internally performs *numerical* forward substitution: it evaluates derivatives level by level (d^0x, d^0y, then d^1x, d^1y using d^0x, etc.) using the stored symbolic `RationalFunction`s and a map of already computed numerical derivative values.
        c. Automatic Differentiation: Wrap the `compute_Y_numerical` function using an Automatic Differentiation framework (like Ceres). When evaluated with `ceres::Jet` types as inputs for the parameters in `θ_identifiable`, this wrapper efficiently computes the numerical sensitivity matrix S_i = ∂Y / ∂θ_identifiable^T without constructing the full symbolic Jacobian.

    Step 3: Iterative Identifiability Assessment
        Start an iterative loop: while True.
        a. Construct Overall Sensitivity Matrix S:
            Initialize an empty list S_list.
            For each test point P_i = (θ_params, x_0^(i)) (i = 1 to M):
                Numerically evaluate the extended output vector Y(t, P_i) using the pre-computed derivative expressions (from Step 2) or direct numerical evaluation. Assign a specific time point t (often t=0 is convenient if derivatives are expressed in terms of x(0) and θ, but other choices are possible).
                Compute the local sensitivity matrix S_i = ∂Y / ∂θ_identifiable^T evaluated at P_i. The derivatives are taken only with respect to the parameters currently in θ_identifiable. Use Automatic Differentiation (AD) libraries or numerical finite differences. Ensure columns correspond consistently to parameters in θ_identifiable.
                Append S_i to S_list.
            Stack the matrices vertically to form the overall sensitivity matrix: S = vstack(S_list).
        b. Rank and Nullspace Analysis:
            Compute the Singular Value Decomposition (SVD) of S.
            Determine the numerical rank r of S by counting singular values greater than ε_rank * max(singular_values).
            If r equals the current number of parameters in θ_identifiable, then all remaining parameters are deemed locally identifiable at the tested points. Break the while loop.
            If r is less than the number of parameters in θ_identifiable, compute the basis vectors Z for the nullspace of S (corresponding to small singular values).
        c. Identify and Fix Non-Identifiable Parameter:
            Examine the nullspace vectors Z. Identify the columns j of S (corresponding to parameters θ_j in θ_identifiable) that have a significant contribution to the nullspace vectors (e.g., check if the norm of the j-th element across all nullspace vectors exceeds ε_nullspace).
            If no such parameters are clearly identified (e.g., nullspace involves complex combinations), report ambiguity and break the loop.
            Select one such parameter θ_fix (e.g., the one with the largest contribution).
            Remove θ_fix from θ_identifiable.
            Add θ_fix to θ_fixed, storing its numerical value (e.g., from P_1).
            Continue to the next iteration of the while loop (implicitly reducing the dimensionality for the next sensitivity calculation).

    Step 4: Determine Minimal Required Output Derivatives
        Let S_final be the sensitivity matrix computed in the last iteration of Step 3 (corresponding only to identifiable parameters θ_identifiable), potentially recomputed using all derivatives up to N. Let r_max be its rank.
        Initialize a dictionary n_orders mapping each observable index k to N.
        a. Reduce Overall Order: Iteratively decrease n_test from N-1 down to 0. In each step, construct a view S_view of S_final containing only rows corresponding to derivatives up to order n_test for all observables. If rank(S_view) < r_max, then n_max = n_test + 1 is the highest derivative order needed across all outputs. Break this reduction step. If n_test reaches 0 and rank is still r_max, set n_max = 0. Update n_orders to set all values to n_max.
        b. Refine Individual Orders: Start a loop: while True. Set improvement_found = false.
            Iterate through each observable k = 1 to N_obs.
            If n_orders[k] > 0:
                Temporarily decrease the order for observable k: n_temp = n_orders[k] - 1.
                Construct S_view using rows corresponding to derivative orders {n_orders[1], ..., n_temp (for k), ..., n_orders[N_obs]}.
                If rank(S_view) == r_max:
                    Permanently update n_orders[k] = n_temp.
                    Set improvement_found = true.
                    Break the inner loop (over k) and restart the while loop.
            If the inner loop completes without improvement_found = true, break the while loop.

    Step 5: Output
        Return:
            The list θ_identifiable of parameters deemed locally identifiable.
            The dictionary θ_fixed mapping locally non-identifiable parameters to the numerical value they were fixed at during analysis.
            The dictionary n_orders mapping each observable index k to the minimum required derivative order n_k.

5. Implementation Notes and Limitations

    The algorithm relies heavily on numerical libraries for linear algebra (SVD) and sensitivity calculation (Automatic Differentiation is highly recommended for accuracy and efficiency compared to finite differences).
    Symbolic math toolboxes can greatly facilitate the pre-computation of derivatives (Step 2).
    The choice of random points P_i and tolerances (ε_rank, ε_nullspace) can influence the results. Running the analysis multiple times with different random seeds is advisable.
    This method assesses local identifiability at the tested points. It does not guarantee structural identifiability (identifiability for almost all parameter values). However, consistent local non-identifiability across multiple points strongly suggests structural or practical non-identifiability relevant to parameter estimation.
    The algorithm identifies parameters one by one for fixing. The specific order can depend on numerical factors and which parameter exhibits the "strongest" contribution to the nullspace in each iteration. The final set of identifiable/non-identifiable parameters should ideally be independent of this order, but the specific values in θ_fixed might vary.

6. Required Data Structures and Library Integration

    To ensure consistency between identifiability analysis and parameter estimation, a common representation for the observed ODE system is crucial.

    a. `Observable` Struct:
        - Define a simple struct `Observable { std::string name; /* ... equality, hash ... */ };`.
        - This provides a unique, named identifier for each measured output, avoiding reliance on vector ordering.

    b. `ObservedOdeSystem` Struct:
        - Encapsulates the complete system definition:
            - ODEs `f`: `std::vector<RationalFunction<double>> equations`
            - State Variables `x`: `std::vector<Variable> state_variables`
            - Model Parameters `p`: `std::vector<Variable> parameters`
            - Observable Definitions `g`: `std::map<Observable, RationalFunction<double>> observable_definitions`
        - This struct serves as the primary input for both `IdentifiabilityAnalyzer` and `ParameterEstimationProblem`.

    c. `ExperimentalData` Struct (Revised):
        - Needs to associate measurements with named observables:
            - `std::vector<double> times`
            - `std::map<Observable, std::vector<double>> measurements`
        - Used primarily by `ParameterEstimationProblem`.

    d. `IdentifiabilityAnalyzer` Class:
        - Input: `ObservedOdeSystem`, list of `Variable`s to test (parameters + unknown ICs).
        - Internally computes and stores symbolic derivatives of `x` and `g_k`.
        - Implements the `compute_Y_numerical` function and its AD wrapper.
        - Performs iterative SVD/nullspace analysis using an external linear algebra library (e.g., Eigen).
        - Manages `θ_identifiable`, `θ_fixed`, `n_orders` results.

    e. `ParameterEstimationProblem` Class (Revised):
        - Input: `ObservedOdeSystem`, `ExperimentalData`, lists of parameters/ICs to estimate vs. fix.
        - Cost function accesses observable definitions from `ObservedOdeSystem` and aligns them with data from `ExperimentalData` using `Observable` keys.

    f. External Dependencies:
        - Requires a numerical linear algebra library (Eigen recommended) for SVD, rank, nullspace calculation.
        - Continues to rely on Ceres for Automatic Differentiation.
        - Continues to rely on Boost.Odeint (implicitly via ODE solver utilities).