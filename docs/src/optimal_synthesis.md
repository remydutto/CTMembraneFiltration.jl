# Optimal Synthesis

This section shows how to construct optimal synthesis associated to Problem (OCP). Roughtly speaking, an optimal synthesis corresponds to a partition of the state space into regions where different control strategies are optimal. To do this, we need to compute the singular locus and the switching locus. 

!!! remark
    In this page, some results are given quickly and without detailed explanations. If you want to deeply understand the mathematical background, please refer to the original paper [Dutto et al., 2026](https://hal.science/hal-05493075).

First of all, let's import the required packages for this page.

```@example main
# Import required packages for optimal control analysis
using Plots              # For plotting and visualization
using OptimalControl     # For optimal control tools
using ForwardDiff        # For automatic differentiation
using Roots              # For root finding algorithms
using OrdinaryDiffEq     # For solving differential equations
```

# Model definition

Let remark that the (OCP) problem is affine with respect to the control. The state and the cost dynamic can be written as

```math
\dot x = f(Rc) + u \cdot g(Rc)
```
where ``x = [e, R_c, v]`` is the augmented state, and where functions ``f = [f_0, f_1, f_2]`` and ``g = [g_0, g_1, g_2]`` are given in the defintion of (OCP) and defined in the follwing code. 

```@example main
# Model parameters for the membrane filtration system
μ = 1; R₀ = 1; η = 1; Jf = 1; Jb = 1; C = 1; b = 1; ω = 1; A = 1;
vf = 10;

# Define the vector field f(Rc) for the membrane filtration dynamics
# This function represents the system dynamics with control u = +1
f(Rc) = [
    (μ * (R₀ + Rc) / (2 * η))*(Jf^2 + Jb^2),
    (Jf * b * C - Jb * ω * Rc)/2,
    A*(Jf - Jb)/2
]

# Define the vector field g(Rc) for the membrane filtration dynamics
# This function represents the system dynamics with control u = -1
g(Rc) = [
    (μ * (R₀ + Rc) / (2 * η))*(Jf^2 - Jb^2),
    (Jf * b * C + Jb * ω * Rc)/2,
    A*(Jf + Jb)/2
]

nothing; # hide
```

# Singular locus

The goal now is to compute the singular locus. To do this, we need first to analyse the system to determine candidates for singular states. 

Since the system is affine with respect to the control, the hamitonian is also affine with respect to the control, and can be written as

```math
H(R_c, p, u) = H_0(R_c, p) + u \cdot H_1(R_c, p),
```
where 
```math
H_0(R_c, p) = p \cdot f(R_c) \quad \text{and} \quad H_1(R_c, p) = p \cdot g(R_c)
```
where ``p = [-1, p_1, p_2]`` is the augmented costate.

Thanks to the maximization condition of the Pontryagin maximum principle, the optimal control is given for almost all ``t \in [t_0, t_f]`` by 

```math 
u(t) \left\{ \begin{array}{ll}
 = +1 & \text{if } H_1(R_c(t), p(t)) > 0 \\[0.5em]
 = -1 & \text{if } H_1(R_c(t), p(t)) < 0 \\[0.5em]
\in [-1, 1] & \text{if } H_1(R_c(t), p(t)) = 0
\end{array} \right.
```

A singular arc corresponds to an interval of time `` I \subset [t_0, t_f]`` such that ``H_1(R_c(t), p(t)) = 0`` for all ``t \in I``. Since the final time ``t_f`` is free, one has ``H(R_c(t_f), p(t_f)) = 0`` for all ``t \in [t_0, t_f]``. By using 

```math
H_0(R_c, p) = H_1(R_c, p) = 0
```

one can deduce that ``p_1(t) = \alpha(R_c(t))`` and ``p_2(t) = \beta(R_c(t))`` for all ``t \in I``, where functions ``\alpha`` and ``\beta`` are defined on the code below.

Since ``H_1(R_c(t), p(t)) = 0`` on the singular arc, one has 
```math 
\frac{\mathrm d H_1}{\partial t}(R_c(t), p(t)) = \{H_0, H_1\}(R_c(t), p(t)) = 0 \quad \text{for all } t \in I
``` 

The following code finds the singular state ``R_c^\star`` such that 
```math
\psi(R_c^\star) = \{H_0, H_1\}(R_c^\star, p_s(R_c^\star)) = 0 \quad \text{with } p_s(R_c) = [-1, \alpha(R_c), \beta(R_c)].
``` 

```@example main
# Compute the main functions
Δ(Rc) = f(Rc)[2]*g(Rc)[3] - f(Rc)[3]*g(Rc)[2]
α(Rc) = (f(Rc)[1]*g(Rc)[3] - f(Rc)[3]*g(Rc)[1]) / Δ(Rc)
β(Rc) = (f(Rc)[2]*g(Rc)[1] - f(Rc)[1]*g(Rc)[2]) / Δ(Rc)

# Compute the derivative of β with respect to Rc using automatic differentiation
dβ(Rc) = ForwardDiff.derivative(β, Rc)

# Compute the function ψ(Rc) = -dβ/Δ, used to find the singular state
ψ(Rc) = -dβ(Rc) /Δ(Rc)

# Find the singular state by solving ψ(Rc) = 0
Rc_star = Roots.find_zero(ψ, 2.5)

# Create a visualization plot of the ψ function
plt = plot(xlim = (0, 10), ylim = (-1, 1))
plot!(plt, ψ, 0, 10, label = "ψ")
plot!(plt, [0, 10], [0, 0], label = nothing, color = :black, ls = :dash)
scatter!(plt, [Rc_star], [0], label = "Singular state", color = :red)
```

Since the state ``R_c`` is constant on the singular arc, the singular control cancels ``\dot R_c`` which leads to 

```math 
u(t) = u_s(R_c^\star) \quad \text{where} \quad u_s(R_c) = -\frac{f_1(R_c)}{g_1(R_c)}
```

We can now easily compute the Hamiltonians ``H_+``, ``H_-`` and ``H_s`` and the flow ``\phi_+``, ``\phi_-`` and ``\phi_s`` respectively associated to ``u = +1``, ``u = -1`` and ``u = u_s(R_c)``.

The singular locus corresponds to the set of points in the state space where there can exsits a singular arc. To find this set, it rest only to find the value ``v^\star`` on the ``v`` space where it is optimal to switch from ``u = u_s(R_c^\star)`` to ``u = +1``. To do this, and by using the transversality condition of the maximum principle, such points ``v^\star`` correponds to the root of the function ``S`` defined below. 

```@example main
# Define the Hamiltonian function H(x,p,u) = p'*(f(x[2]) + u*g(x[2]))
H(x,p,u) = p'*(f(x[2]) + u*g(x[2]))

# Define the singular control function uₛ(Rc) = -f₂(Rc)/g₂(Rc)
uₛ(Rc) = -f(Rc)[2]/g(Rc)[2]

# Define the Hamiltonian associated to u = -1, u = +1 and u = uₛ(x)
H₋(x,p) = H(x,p,-1)
H₊(x,p) = H(x,p,+1)
Hₛ(x,p) = H(x,p,uₛ(x[2]))

# Create flow objects for each Hamiltonian
ϕ₋ = Flow(OptimalControl.Hamiltonian(H₋))
ϕ₊ = Flow(OptimalControl.Hamiltonian(H₊))
ϕₛ = Flow(OptimalControl.Hamiltonian(Hₛ))

# Function to find the end of the singular arc
function S(v0)
    x0 = [0, Rc_star, v0]                                       # Initial state
    p0 = [-1, α(Rc_star), β(Rc_star)]                           # Initial costate
    cb = ContinuousCallback((z,_,_) -> z[3] - vf, terminate!)   # Stop when v = vf
    xf, pf = ϕ₊(0, x0, p0, 100, callback = cb)                  # Integrate forward
    return pf[2]                                                # Return final costate
end

# Compute the end of the singular arc by finding the zero of the function S
v_star = Roots.find_zero(S, 5)

# Plot the function S
plt = plot(xlim = (0, 9.5))
plot!(plt, S, label = "S")
plot!(plt, [0,10], [0,0], label = nothing, color = :black, ls = :dash)
scatter!(plt, [v_star], [0], label = "Singular state", color = :red)
```

# Switching locus

The goal now is to construct the switching locus ``\mathcal S``, which is the set of points ``(R_c, v)`` in the state space where the optimal control switches from ``u = -1`` to ``u = +1``. If ``(R_c, v) \in \mathcal S``, then one has ``S(R_c, v) = 0``, where the function ``S \colon \mathbb R^2 \to \mathbb R`` is defined below. Moreover, one can use differential continuation method to find the switching locus. 

To do this, let assume that there exists a function ``\gamma`` such that points ``(Rc, v)`` satisfy ``\gamma(Rc) = v``. One has thus ``S(Rc, γ(Rc)) = 0``, and then 

```math
\frac{\partial S}{\partial R_c}(R_c, \gamma(R_c)) + \frac{\partial S}{\partial v}(R_c, \gamma(R_c)) \gamma'(R_c) = 0, 
```
which gives

```math
\gamma'(R_c) = F(R_c, \gamma(R_c)) = -\left(\frac{\partial S}{\partial v}(R_c, \gamma(R_c))\right)^{-1} \frac{\partial S}{\partial R_c}(R_c, \gamma(R_c)).
```
The differential continuation method consists to solve the follwing cauchy problem 

```math
\gamma'(R_c) = F(R_c, \gamma(R_c)) \quad \gamma(R_c^\star) = v^\star, \quad R_c \geq R_c^\star.
```
This is done in the following code.

```@example main
# Define the plot limits for the synthesis
xlim = [0, vf]; ylim = [2, 5];

# Define the function S(Rc,v) such that if (Rc,v) is on the switching locus, then S(Rc,v) = 0
function S(Rc,v)
    x0 = [0., Rc, v]                                            # Initial state
    p0 = [-1, α(Rc), β(Rc)]                                     # Initial costate
    cb = ContinuousCallback((z,_,_) -> z[3] - vf, terminate!)   # Stop when v = vf
    xf, pf = ϕ₊(0, x0, p0, 100, callback = cb)                  # Integrate forward
    return pf[2]                                                # Return final costate
end

# Compute partial derivatives of S using automatic differentiation
dSdRc(Rc, v) = ForwardDiff.derivative(Rc -> S(Rc,v), Rc)        # Derivative w.r.t Rc
dSdv(Rc, v) = ForwardDiff.derivative(v -> S(Rc,v), v)           # Derivative w.r.t v

# Define the differential equation associated to the switching locus: dv/dRc = -dSdRc/dSdv
lhs(v,_,Rc) = -dSdRc(Rc, v)/dSdv(Rc, v)

# Solve the ODE to compute the switching locus from Rc_star to ylim[2]
prob = ODEProblem(lhs, v_star, (Rc_star, ylim[2]), saveat=0.01)
SL = solve(prob, Tsit5())

nothing; # hide
```

The following code simply creates the optimal feedback synthesis, based on the construction of the singular locus and the switching locus.

```@example main
# Create the synthesis plot showing the optimal control regions
synthesis = plot(xlim = xlim, ylim = ylim, xlabel = "v", ylabel = "Rc", legend = :topleft)

# Plot the region where control u = -1 is optimal (red shape)
plot!(synthesis,
    [xlim[1]; v_star;  SL.u; xlim[1]],
    [Rc_star; Rc_star; SL.t; ylim[2]],
    st=:shape, linecolor = nothing, label = "u = -1", fillcolor = :red, fillalpha = 0.35)

# Plot the region where control u = +1 is optimal (blue shape)
plot!(synthesis,
    [xlim[1]; v_star;  SL.u; xlim[2]; xlim[2]; xlim[1]],
    [Rc_star; Rc_star; SL.t; ylim[2]; ylim[1]; ylim[1]],
    st=:shape, linecolor = nothing, label = "u = +1", fillcolor = :blue, fillalpha = 0.35)

# Plot the singular locus (green line)
plot!(synthesis, [xlim[1], v_star], [Rc_star, Rc_star], label = "Singular locus", c=:green, lw=3)

# Plot the switching locus (gold line)
plot!(synthesis, SL.u, SL.t, label = "Switching locus", c=:gold, lw=3)
```
