# Optimal Synthesis

This section demonstrates the setup process for optimal synthesis analysis. The following code briefly loads the required packages for mathematical computations and visualization.

```@example main
using Plots
using OptimalControl
using ForwardDiff
using Roots
```

# Model definition 

This section covers the mathematical modeling phase of optimal synthesis. The following code briefly defines the system parameters and dynamics functions for the membrane filtration model.

```@example main
μ = 1; R₀ = 1; η = 1; Jf = 1; Jb = 1; C = 1; b = 1; ω = 1; A = 1;
f(Rc) = [ 
    (μ * (R₀ + Rc) / (2 * η))*(Jf^2 + Jb^2),
    (Jf * b * C - Jb * ω * Rc)/2,
    A*(Jf - Jb)/2
]
g(Rc) = [ 
    (μ * (R₀ + Rc) / (2 * η))*(Jf^2 - Jb^2),
    (Jf * b * C + Jb * ω * Rc)/2,
    A*(Jf + Jb)/2
]
```

# Singular locus

This section addresses the singular control analysis for optimal synthesis. The following code briefly computes the singular locus and generates a visualization of the switching function behavior.

```@example main
Δ(Rc) = f(Rc)[2]*g(Rc)[3] - f(Rc)[3]*g(Rc)[2]
α(Rc) = (f(Rc)[1]*g(Rc)[3] - f(Rc)[3]*g(Rc)[1]) / Δ(Rc)
β(Rc) = (f(Rc)[2]*g(Rc)[1] - f(Rc)[1]*g(Rc)[2]) / Δ(Rc)

dβ(Rc) = ForwardDiff.derivative(β, Rc)
ψ(Rc) = -dβ(Rc) /Δ(Rc)

Rc_sing = Roots.find_zero(ψ, 2.5)
println("Singular state : ", Rc_sing)

plt = plot(xlim = (0, 10), ylim = (-1, 1))
plot!(plt, ψ, 0, 10, label = "ψ")
plot!(plt, [0, 10], [0, 0], label = nothing, color = :black, ls = :dash)
scatter!(plt, [Rc_sing], [0], label = "Singular state", color = :red)
```
