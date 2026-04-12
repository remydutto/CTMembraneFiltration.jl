# Control-Toolbox

## Studied Problem

We now consider a more complex problem where there are two types of masses attached to the membrane, described by two resistances ``r_1`` and ``r_2``. 

Since the code provided in the CTMembraneFiltration.jl package only considers the case of one resistance/mass, this page proposes to show how one can use the [control-toolbox](https://control-toolbox.org/) ecosystem through the OptimalControl.jl package to obtain the optimal solution of the following problem.

```@example setup
using Plots
using OptimalControl

# Problem parameters
t0 = 0
vf = 10
x0 = [0., 0., 0.1, 0.1]
p = [1., 1., 1., 1., 2., 1.5]

# Model functions
e(x) = p[1] + x[3] + x[4]
g(x) = p[2]
f11(x) = p[3]
f12(x) = - p[4]*x[3]
f21(x) = p[5]
f22(x) = - p[6]*x[4]

# Vector fields
Ff(x) = [e(x); g(x); f11(x); f21(x)]  # filtration
Fb(x) = [e(x); -g(x); f12(x); f22(x)] # backwash
F0(x) = Ff(x) + Fb(x)
F1(x) = Ff(x) - Fb(x)
```

### Optimal Control Problem Definition

```@example ocp_definition
@def begin
    tf  > t0, variable
    t   in [t0, tf], time
    x   = (e, v, r1, r2) in R^4, state
    u   in R, control
    -1  <= u(t) <= 1

    x(t0) == x0
    v(tf) == vf

    xdot(t) == F0(x(t)) + u(t)*F1(x(t))

    e(tf) -> min
end
```

**Abstract definition:**
```
tf  > t0, variable
t   in [t0, tf], time
x   = ((e, v, r1, r2) in R^4, state)
u   in R, control
-1  <= u(t) <= 1

x(t0) == x0
v(tf) == vf

xdot(t) == F0(x(t)) + u(t) * F1(x(t))

e(tf) -> min
```

The (autonomous) optimal control problem is of the form:
minimize ``J(x, u, tf) = g(x(0.0), x(tf), tf)``

subject to:
```
xdot(t) = f(x(t), u(t), tf), t in [0.0, tf] a.e.,
phi- <= phi(x(0.0), x(tf), tf) <= phi+,
u- <= u(t) <= u+,
```

where ``x(t) = (e(t), v(t), r1(t), r2(t)) in R^4``, ``u(t) in R`` and ``tf in R``.

## Direct Method

A classical optimal control method used to obtain a solution of the problem is the direct method. The main idea is to discretize in time the studied optimal control problem, to formulate the associated large optimization problem and to solve it.

Thanks to the OptimalControl.jl package, one can easily obtain and plot the solution of the direct method by using the `solve` function. Here we use the interior point Ipopt solver to obtain a solution of the discretized problem.

```@example direct_solve
using NLPModelsIpopt

# Solve by direct method
direct_sol = solve(ocp; display = true)

# Visualization
plt_sol = plot(direct_sol, label = "direct")
```

## Structure of the Solution

The analysis of the solution obtained by the direct method reveals a typical arc structure:

```@example structure_analysis
# Analyze solution structure
structure = analyze_control_structure(direct_sol)

# Identify arcs
arcs = identify_arcs(direct_sol)
println("Identified structure: $(arcs)")
```

The optimal solution typically presents the following structure:
1. **Filtration arc** (``u = +1``): Production with progressive resistance increase
2. **Backwash arc** (``u = -1``): Regeneration with resistance decrease
3. **Singular arcs**: Optimal transitions between modes

### Structure Visualization

```@example structure_plot
# Plot control structure
plot_control_structure(direct_sol)

# Plot states
plot_state_evolution(direct_sol)
```

## Indirect Shooting

The indirect method is based on Pontryagin's maximum principle. It consists of:

1. Writing the optimality conditions
2. Solving the resulting boundary value problem
3. Using shooting methods to converge to the solution

### Optimality Conditions

```@example indirect_conditions
# Compute Hamiltonian
H(x, p, u) = e(x) + u * g(x) + p' * (F0(x) + u * F1(x))

# Optimal control
u_opt(x, p) = sign(g(x) + p' * F1(x))

# Hamiltonian system
function hamiltonian_system(x, p)
    u = u_opt(x, p)
    xdot = F0(x) + u * F1(x)
    pdot = -gradient(x -> H(x, p, u), x)[1]
    return xdot, pdot
end
```

### Shooting Resolution

```@example shooting
# Configure shooting
shooting_problem = configure_shooting(ocp, hamiltonian_system)

# Solve
indirect_sol = solve_shooting(shooting_problem)

# Compare with direct solution
plot_comparison(direct_sol, indirect_sol)
```

## Turnpike Property

The turnpike property describes the behavior of optimal solutions when the time horizon becomes large. In this context, we observe that:

1. Trajectories spend the majority of time near an equilibrium state
2. Fast transitions occur at the beginning and end of the horizon
3. The optimal structure converges to a periodic pattern

### Turnpike Analysis

```@example turnpike_analysis
# Compute equilibrium state
equilibrium = compute_equilibrium(ocp)

# Analysis for different time horizons
horizons = [5.0, 10.0, 20.0, 50.0]
solutions = [solve_with_horizon(ocp, T) for T in horizons]

# Visualize turnpike property
plot_turnpike_behavior(solutions, horizons, equilibrium)
```

### Physical Interpretation

The turnpike property has a natural interpretation for membrane filtration systems:

- **Equilibrium state**: Optimal permanent operating regime
- **Transitions**: Optimal startup and shutdown periods
- **Asymptotic behavior**: Optimal strategy for continuous operations

## Method Comparison

### Comparison Criteria

| Criterion | Direct Method | Indirect Method |
|-----------|----------------|------------------|
| Convergence | Robust | Sensitive to initial conditions |
| Precision | Moderate | High |
| Computation time | Fast | Variable |
| Structure analysis | Difficult | Natural |

### Numerical Results

```@example comparison
# Quantitative comparison
comparison_results = compare_methods(direct_sol, indirect_sol)

println("Direct cost: $(direct_sol.cost)")
println("Indirect cost: $(indirect_sol.cost)")
println("Relative error: $(comparison_results.relative_error)")
```

## Practical Recommendations

### When to use the direct method?

- Problems with complex constraints
- Need for numerical robustness
- First exploration of the problem

### When to use the indirect method?

- Deep theoretical analysis
- Need for high precision
- Sensitivity studies

### Hybrid Approach

An effective strategy consists of:
1. Using the direct method for a first estimate
2. Refining with the indirect method
3. Validating results by sensitivity analysis

```@example hybrid_approach
# Hybrid approach
hybrid_sol = hybrid_solve(ocp)
plot_hybrid_solution(hybrid_sol)
```

## Extensions and Perspectives

### Possible Extensions

1. **More complex models**: Inclusion of additional nonlinearities
2. **State constraints**: Limits on resistances
3. **Multi-objective optimization**: Energy cost vs water quality
4. **Stochastic control**: Parameter uncertainties

### Research Perspectives

- Development of adaptive algorithms
- Stability analysis of optimal solutions
- Extension to multi-membrane systems
- Integration with experimental data

## References

- Martinon, P., & Gergaud, J. (2007). *Using shooting methods for solving optimal control problems with discontinuities*.
- Bonnans, J. F., et al. (2017). *Optimal control with engineering applications*.
- Control-toolbox ecosystem: https://control-toolbox.org/
