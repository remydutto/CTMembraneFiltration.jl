# Introduction

Production–regeneration systems describe processes that must alternate between a productive phase and a recovery phase. During production, the system generates an output (for instance filtered water), but this phase also progressively degrades an internal state (such as fouling or resistance). The regeneration phase does the opposite: it restores the system’s internal condition, but does not directly produce useful output.

This package focuses on the analysis and optimization of such systems, particularly in the context of membrane filtration processes where fouling accumulates during filtration and must be periodically removed through regeneration (e.g., backwashing). Such system is schematized through the following figure:

```@raw html
    <img src="./assets/backwash.jpg" width="800px">
```

# Description of membrane filtration system

The goal now is to briefly describe the modelling of a membrane filtration process, in order to introduce two examples of problems which can be handled with [CTMembraneFiltration.jl](https://remydutto.github.io/CTMembraneFiltration.jl) package.

The permeate flux ``J~(\mathrm m. \mathrm s^{-1})`` corresponds to the permeate volume ``v~(\mathrm m^3)`` flow per membrane area ``A~(\mathrm m^2)`` and is related to the total resistance ``R~(\mathrm m^{-1})`` of the membrane according to the Darcy's law 
```math
    J = \frac{\dot v}{A} = \frac{\Delta p}{\mu R}
```
where ``\Delta p~(\mathrm{Pa})`` corresponds to the trans-membrane pressure and ``\mu~(\mathrm{Pa}.\mathrm s)`` is the permeate viscosity. The total resistance of the membrane corresponds to the sum between the the intrinsic resistance ``R_0~(\mathrm m^{-1})`` and the residual cake resistance ``R_c~(\mathrm m^{-1})`` thanks to the "resistances-in-series" concept
```math
    R = R_0 + R_c. 
```
The pump consumes an energy ``e~(\mathrm W.\mathrm h)`` flow, whose variation is proportional to the hydraulic power ``P_h~(\mathrm W)``, the later being the product of the permeate flow ``J`` by the trans-membrane pressure ``\Delta p``
```math 
    \dot e = \frac{J\Delta p}{\eta},
```
where ``\eta`` is the efficiency parameter of the pump.

The dynamics of the resistance ``R_c`` rely on complex phenomena, but an efficient simple formulation is given by the following: 
- In filtration phases (when ``u = +1``), this resistance variation is proportional to the flux ``J`` and the total suspended solids concentration ``C~(\mathrm{Kg}.\mathrm m^{-3})`` 
```math
    \dot R_c = J \beta C \quad \text{when} \quad u = +1,
```
where ``\beta~(\mathrm m.\mathrm{Kg}^{-1})`` is a parameter that describes the resistance of the cake layer. 
- In backwash phases (when ``u = -1``), this resistance variation is proportional to the flux ``J`` and the resistance cake layer ``R_c``
```math
    \dot R_c = J \omega R_c  \quad \text{when} \quad u = -1,
```
where ``\omega~(\mathrm m^{-1})`` models the detachment resistance of fouling. 

# Problem statement

Let us assume that the flux $J$ is constant during filtration (``J = J_f > 0``) and backwash (``J = -J_b < 0``) phases. The goal is to minimize the total power used to produce a targeted permeate volume ``v_f`` at a free final time ``t_f``. Using equations given previously, the dynamics of the cost and the state variables are given in filtration and backwash modes by
```math
\mathrm{Filtration} : \left\{ \begin{array}{rl}
\dot e & = \frac{J_f^2 \mu}{\eta} (R_0 + R_c(t)),
\\[0.5em]
\dot R_c & = J_f \beta C, 
\\[0.5em]
\dot v & = J_f A, 
\end{array} \right. 
\quad \text{and} \quad 
\mathrm{Backwash} : \left\{ \begin{array}{rl}
\dot e & = \frac{J_b^2 \mu}{\eta} (R_0 + R_c), 
\\[0.5em]
\dot R_c & = -J_b \omega R_c,
\\[0.5em]
\dot v & = - J_b A.
\end{array} \right.  
```

Let ``u \in [-1, 1]`` be the control variable, where `` u = +1`` denotes filtration mode and ``u = -1`` denotes backwash mode, we consider the following optimal control problem

```math
\text{(OCP)} \quad
\left\{ 
\begin {array}{ll}
\displaystyle \min_{x,y,t_f} \int_{t_0}^{t_f} \frac{\mu(R_0 + R_c)}{2\eta}(J_f^2 + J_b^2 + u(t)(J_f^2 - J_b^2)) \, \mathrm dt, 
& t \in [t_0, t_f] \ \mathrm{a.e.}, \\[1em]
\displaystyle \mathrm{s.t.} \ \dot x(t) = \frac{1 + u(t)}{2} f_p(x(t)) + \frac{1-u(t)}{2} f_r(x(t)), & t \in [t_0, t_f] \ \mathrm{a.e.}, \\[1em] 
\displaystyle \phantom{\mathrm{s.t.} \ } \dot y(t) = \frac{1 + u(t)}{2} g_p(x(t)) + \frac{1-u(t)}{2} g_r(x(t)), \, & t \in [t_0, t_f] \ \mathrm{a.e.}, \\[1em]
\phantom{\mathrm{s.t.} \ } u(t) \in [-1, 1], & t \in [t_0, t_f], \\[1em]
\phantom{\mathrm{s.t.} \ } x(t_0) = x_0, \quad y(t_0) = y_0, \quad y(t_f) = T,
\end{array}
\right.
```
where ``t_0 \in \mathbb R``, ``x_0 > 0 ``, ``y_0 > 0`` and ``T > 0`` are provided.  

# Main theoretical results

Based on the theoretical developments presented in [Dutto et al., 2026](https://hal.science/hal-05493075), all possible optimal solution structures are characterized by the following result:

!!! tip "Theorem"
    Under standard regularity assumptions, and denoting 
    - ``\sigma_-`` a regeneration arc associated to ``u = -1``,
    - ``\sigma_+`` a production arc associated to ``u = +1``,
    - ``\sigma_s`` a singular arc associated to ``u = u_s``,
    the structure of an optimal solution can only be one of the following:
    ```math
    \sigma_+, \
    \sigma_-\sigma_+, \ 
    \sigma_s\sigma_+, \
    \sigma_-\sigma_s\sigma_+ ~
    \text{or} ~
    \sigma_+\sigma_s\sigma_+.
    ```

Here, ``u_s`` denotes a singular control such that ``\dot x = 0``. See [Dutto et al., 2026](https://hal.science/hal-05493075) for a detailed presentation and proof.

## Reproducibility

```@setup main
using Pkg
using InteractiveUtils
using Markdown

# Download links for the benchmark environment
function _downloads_toml(DIR)
    link_manifest = joinpath("assets", DIR, "Manifest.toml")
    link_project = joinpath("assets", DIR, "Project.toml")
    return Markdown.parse("""
    You can download the exact environment used to build this documentation:
    - 📦 [Project.toml]($link_project) - Package dependencies
    - 📋 [Manifest.toml]($link_manifest) - Complete dependency tree with versions
    """)
end
```

```@example main
_downloads_toml(".") # hide
```

```@raw html
<details style="margin-bottom: 0.5em; margin-top: 1em;"><summary>ℹ️ Version info</summary>
```

```@example main
versioninfo() # hide
```

```@raw html
</details>
```

```@raw html
<details style="margin-bottom: 0.5em;"><summary>📦 Package status</summary>
```

```@example main
Pkg.status() # hide
```

```@raw html
</details>
```

```@raw html
<details style="margin-bottom: 0.5em;"><summary>📚 Complete manifest</summary>
```

```@example main
Pkg.status(; mode = PKGMODE_MANIFEST) # hide
```

```@raw html
</details>
```

# Partners and Fundings

```@raw html
<div style="display: flex; align-items: center; justify-content: center; gap: 20px;">
    <img src="./assets/Logo-INRAE.jpg" width="200px">
    <img src="./assets/Logo-WOC.png" width="200px">
    <img src="./assets/Logo-Defi-Region.jpg" width="200px">
</div>
```