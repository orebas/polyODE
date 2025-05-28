using HomotopyContinuation

# Define the variables
@var a b c x1 d1x1 d2x1 d3x1 d4x1 x2 d1x2 d2x2 d3x2 d4x2 x3 d1x3 d2x3 d3x3 d4x3 d5x3

# Define the polynomial system
f = [
    -540439/65536 + x3,
    -127869/65536 + d1x3,
    -16387/65536 + d2x3,
    -4196/65536 + d3x3,
    -747/65536 + d4x3,
    -159/65536 + d5x3,
    -c*x1 - c*x2 + d1x3,
    -c*d1x1 - c*d1x2 + d2x3,
    -c*d2x1 - c*d2x2 + d3x3,
    -c*d3x1 - c*d3x2 + d4x3,
    -c*d4x1 - c*d4x2 + d5x3,
    a*x1 + d1x1,
    -b*x2 + d1x2,
    a*d1x1 + d2x1,
    -b*d1x2 + d2x2,
    a*d2x1 + d3x1,
    -b*d2x2 + d3x2,
    a*d3x1 + d4x1,
    -b*d3x2 + d4x2
]

# Solve the system multiple times with different settings for better accuracy
println("Solving the polynomial system with improved settings...")

# First solve with default settings
result1 = solve(f, show_progress=true)

# Try to improve solutions using a different start system
println("Attempting to refine solutions...")
result2 = solve(f, show_progress=false, compile=false)

# Take the best residual solutions
println("Selecting best solutions based on residuals...")
vars = [a, b, c, x1, d1x1, d2x1, d3x1, d4x1, x2, d1x2, d2x2, d3x2, d4x2, x3, d1x3, d2x3, d3x3, d4x3, d5x3]

best_result = []
for i in 1:length(result1)
    sol1 = result1[i]
    sol2 = result2[i]
    
    res1 = norm([evaluate(eq, vars => solution(sol1)) for eq in f])
    res2 = norm([evaluate(eq, vars => solution(sol2)) for eq in f])
    
    if res1 <= res2
        push!(best_result, sol1)
    else
        push!(best_result, sol2)
    end
end

result = best_result

# Display the results
println("\nNumber of solutions found: ", length(result))
println("\n" * "="^50)

# Display each solution
for (i, sol) in enumerate(result)
    println("\nSolution $i:")
    println("-"^50)
    
    # Extract values with variable names
    vars = [a, b, c, x1, d1x1, d2x1, d3x1, d4x1, x2, d1x2, d2x2, d3x2, d4x2, x3, d1x3, d2x3, d3x3, d4x3, d5x3]
    var_names = ["a", "b", "c", "x1", "d1x1", "d2x1", "d3x1", "d4x1", "x2", "d1x2", "d2x2", "d3x2", "d4x2", "x3", "d1x3", "d2x3", "d3x3", "d4x3", "d5x3"]
    
    sol_values = solution(sol)
    
    for (j, name) in enumerate(var_names)
        val = sol_values[j]
        if imag(val) ≈ 0
            println("$name = $(real(val))")
        else
            println("$name = $val")
        end
    end
end

# Check if solutions are real
println("\n" * "="^50)
println("\nChecking for real solutions...")
real_solutions = [sol for sol in result if is_real(sol)]
println("Number of real solutions: ", length(real_solutions))

if length(real_solutions) > 0
    println("\nReal solutions:")
    for (i, sol) in enumerate(real_solutions)
        println("\nReal Solution $i:")
        println("-"^50)
        
        var_names = ["a", "b", "c", "x1", "d1x1", "d2x1", "d3x1", "d4x1", "x2", "d1x2", "d2x2", "d3x2", "d4x2", "x3", "d1x3", "d2x3", "d3x3", "d4x3", "d5x3"]
        
        sol_values = solution(sol)
        
        for (j, name) in enumerate(var_names)
            println("$name = $(real(sol_values[j]))")
        end
    end
end

# Verify solutions by substituting back
println("\n" * "="^50)
println("\nVerifying solutions...")
for (i, sol) in enumerate(result)
    vars = [a, b, c, x1, d1x1, d2x1, d3x1, d4x1, x2, d1x2, d2x2, d3x2, d4x2, x3, d1x3, d2x3, d3x3, d4x3, d5x3]
    
    sol_values = solution(sol)
    
    residual = norm([evaluate(eq, vars => sol_values) for eq in f])
    println("Solution $i residual norm: $residual")
end
