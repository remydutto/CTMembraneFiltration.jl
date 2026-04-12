# Optimal Synthesis Construction

## Introduction

The construction of optimal synthesis for production-regeneration systems is a crucial step that allows determining the complete structure of optimal solutions depending on initial and final conditions. This approach provides a global understanding of the system behavior and enables rapid calculation of optimal solutions without having to solve complex optimization problems each time.

## Synthesis Theory

### Fundamental Principle

Optimal synthesis is based on Pontryagin's maximum principle and switching arc analysis. For our production-regeneration system, we characterize all possible types of arcs:

- **Production arc** ``\sigma_+`` : ``u = +1``
- **Regeneration arc** ``\sigma_-`` : ``u = -1``  
- **Singular arc** ``\sigma_s`` : ``u = u_s`` with ``\dot x = 0``

### Hamiltonian Function

The system Hamiltonian is written as:

```math
H(x, p, u) = \frac{1 + u}{2} (l_p(x) + p \cdot f_p(x) + q \cdot g_p(x)) + \frac{1 - u}{2} (l_r(x) + p \cdot f_r(x) + q \cdot g_r(x))
```

where ``p`` and ``q`` are the adjoint variables associated with states ``x`` and ``y``.

### Optimality Conditions

The optimal control is given by:

```math
u^*(x, p) = \begin{cases}
+1 & \text{if } H(x, p, +1) < H(x, p, -1) \\
-1 & \text{if } H(x, p, +1) > H(x, p, -1) \\
u_s & \text{if } H(x, p, +1) = H(x, p, -1)
\end{cases}
```

## Practical Construction

### Step 1: Switching Point Analysis

Switching points between different arc types are determined by the equality of Hamiltonians:

```math
H(x, p, +1) = H(x, p, -1)
```

This condition defines a switching surface in the phase space ``(x, p)``.

### Step 2: Singular Arc Calculation

Singular arcs satisfy:

```math
\frac{d}{dt} \frac{\partial H}{\partial u} = 0
```

For our system, this leads to the condition ``\dot x = 0``, which means the fouling state remains constant during the singular arc.

### Step 3: Structure Assembly

The possible optimal structures are:

1. **Pure production** : ``\sigma_+``
2. **Regeneration followed by production** : ``\sigma_- \sigma_+``
3. **Singular arc followed by production** : ``\sigma_s \sigma_+``
4. **Regeneration, singular arc, production** : ``\sigma_- \sigma_s \sigma_+``
5. **Production, singular arc, production** : ``\sigma_+ \sigma_s \sigma_+``

## Construction Algorithm

```@example synthesis
using CTMembraneFiltration
using Plots

# System parameters
params = (
    lp = x -> x^2,           # production cost
    lr = x -> 0.5*x,         # regeneration cost  
    fp = x -> -0.1*x,        # production dynamics
    fr = x -> 0.2*(1-x),     # regeneration dynamics
    gp = x -> 1.0,           # production
    gr = x -> 0.0            # no production during regeneration
)

# Construct synthesis
synthesis = construct_optimal_synthesis(params)

# Visualize control regions
plot_control_regions(synthesis)
```

### Region Determination

```@example regions
# Compute optimal control regions
regions = compute_control_regions(synthesis)

# Display boundaries
plot_switching_curves(regions)
```

### Optimal Trajectory Calculation

```@example trajectories
# Initial and final points
x0 = 1.5
y0 = 0.0
T = 10.0

# Compute optimal trajectory
optimal_trajectory = compute_optimal_trajectory(x0, y0, T, synthesis)

# Visualization
plot(optimal_trajectory.time, optimal_trajectory.state, 
     label="Optimal trajectory", xlabel="Time", ylabel="State")
```

## Practical Applications

### Case 1: Low Initial Fouling

For ``x_0 < x_s`` (where ``x_s`` is the singular state), the optimal structure is typically ``\sigma_+`` (pure production).

### Case 2: Moderate Fouling

For ``x_s < x_0 < x_c`` (where ``x_c`` is the switching point), the optimal structure is ``\sigma_- \sigma_+``.

### Case 3: High Fouling

For ``x_0 > x_c``, the optimal structure may include a singular arc: ``\sigma_- \sigma_s \sigma_+``.

## Numerical Validation

```@example validation
# Validate synthesis against direct resolution
test_points = [(x0=0.5, T=5.0), (x0=1.0, T=8.0), (x0=1.5, T=10.0)]

for (i, point) in enumerate(test_points)
    # Solution by synthesis
    sol_synthesis = solve_by_synthesis(point.x0, 0.0, point.T, synthesis)
    
    # Direct solution
    sol_direct = solve_direct(point.x0, 0.0, point.T)
    
    # Comparison
    error = compute_error(sol_synthesis, sol_direct)
    println("Point $i: error = $error")
end
```

## Advantages of Synthesis Approach

1. **Computational speed**: Once synthesis is constructed, solutions are obtained instantly
2. **Global understanding**: Provides an overview of all possible solutions
3. **Robustness**: Less sensitive to numerical errors than direct methods
4. **Parametric analysis**: Easy to study as a function of system parameters

## Limitations and Extensions

### Current Limitations

- Regularity assumptions on dynamics
- Simple quadratic costs
- Absence of additional state constraints

### Possible Extensions

- Inclusion of state constraints ``x(t) \leq x_{max}``
- More general costs (non-convex)
- Systems with multiple operating modes
- Parameter uncertainties

## References

- Dutto, R., Cots, O., & Martinon, P. (2026). *Optimal control of production-regeneration systems with application to membrane filtration*. HAL-05493075.
- Bryson, A. E., & Ho, Y. C. (1975). *Applied optimal control*. Taylor & Francis.
- Pontryagin, L. S., et al. (1962). *The mathematical theory of optimal processes*. Interscience Publishers.
