## Introduction

Production–regeneration systems describe processes that must alternate between a productive phase and a recovery phase. During production, the system generates an output (for instance filtered water), but this phase also progressively degrades an internal state (such as fouling or resistance). The regeneration phase does the opposite: it restores the system’s internal condition, but does not directly produce useful output.

Our goal is to showh how the [Control-Toolbox](https://control-toolbox.org/) ecosystem, in particulary the [OptimalControl.jl](https://github.com/Control-toolbox/OptimalControl.jl) package, can be used in the context of membrane filtration processes where fouling accumulates during filtration and must be periodically removed through regeneration (e.g., backwashing). Such system is schematized through the following figure.

```@raw html
    <img src="./assets/backwash.jpg" width="800px">
```

A package dedicated to production-regeneration systems, specialized on such application is available and denoted [Filtration.jl](https://remydutto.github.io/doc-Filtration.jl/dev/).

## Description of membrane filtration system

The goal now is to briefly describe the modelling of a membrane filtration process. The permeate flux ``J~(\mathrm m. \mathrm s^{-1})`` corresponds to the permeate volume ``v~(\mathrm m^3)`` flow per membrane area ``A~(\mathrm m^2)`` and is related to the total resistance ``R~(\mathrm m^{-1})`` of the membrane according to the Darcy's law 
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

## Problem statement

Let us assume that the flux $J$ is constant during filtration (``J = J_f > 0``) and backwash (``J = -J_b < 0``) phases. The goal is to minimize the total power used to produce a targeted permeate volume ``v_f`` at a free final time ``t_f``. Using equations given previously, the dynamics of the cost and the state variables are given in filtration and backwash modes by
```math
\mathrm{Filtration} : \left\{ \begin{array}{rl}
\dot e & = \frac{J_f^2 \mu}{\eta} (R_0 + R_c),
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
\displaystyle \min_{x,y,t_f} \int_{t_0}^{t_f} \frac{\mu(R_0 + R_c(t))}{2\eta}\big(J_f^2 + J_b^2 + u(t)(J_f^2 - J_b^2)\big) \, \mathrm dt, &  \\[1em]
\displaystyle \mathrm{s.t.} \ \dot R_c(t) = J_f\beta C - J_b \omega R_c(t) + u(t)\big(J_f \beta C + J_b \omega R_c(t)\big), & t \in [t_0, t_f] \ \mathrm{a.e.}, \\[1em] 
\displaystyle \phantom{\mathrm{s.t.} \ } \dot v(t) = A\big((J_f - J_b) + u(t)(J_f + J_b)\big), \, & t \in [t_0, t_f] \ \mathrm{a.e.}, \\[1em]
\phantom{\mathrm{s.t.} \ } u(t) \in [-1, 1], & t \in [t_0, t_f], \\[1em]
\phantom{\mathrm{s.t.} \ } R_c(t_0) = R_{c0}, \quad v(t_0) = v_0, \quad v(t_f) = v_f,
\end{array}
\right.
```
where ``t_0 \in \mathbb R``, ``R_{c0} > 0``, ``v_0 > 0`` and ``v_f > 0`` are provided.  

## Main theoretical results

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

## References

- [Filtration.jl](https://remydutto.github.io/doc-Filtration.jl/dev/)
- [OptimalControl.jl](https://juliadynamics.github.io/OptimalControl.jl/dev/)
- Rémy Dutto, Jérôme Harmand, Alain Rapaport (2026). [Optimal control synthesis for a class of production-regeneration systems -Application to membrane filtration](https://hal.science/hal-05493075).
- F. Aichouche, N. Kalboussi, A. Rapaport, J. Harmand (2020). [Modeling and optimal control for production-regeneration systems - preliminary results -](https://ieeexplore.ieee.org/document/9143741), _2020 European Control Conference (ECC)_
- B. Benyahia, A. Charfi, N. Benamar, M. Heran, A. Grasmick, B. Cherki, J. Harmand (2013). [A simple model of anaerobic membrane bioreactor for control design: coupling the “AM2b” model with a simple membrane fouling dynamics](https://www.researchgate.net/publication/272506325_A_simple_model_of_anaerobic_membrane_bioreactor_for_control_design_coupling_the_AM2b_model_with_a_simple_membrane_fouling_dynamics), _World Congress on Anerobic Digestion: Recovering (bio) Ressources for the World_
- F. Ellouze, Y. Kammoun, N. Kalboussi, A. Rapaport, J. Harmand, S. Nasr, N. Ben Amar (2023) [Optimal control of backwash scheduling for pumping energy saving: Application to the treatment of urban wastewater](https://www.sciencedirect.com/science/article/abs/pii/S221471442300898X), _Journal of Water Process Engineering_
- N. Kalboussi, A. Rapaport, T. Bayen, N. Ben Amar, F. Ellouze, J. Harmand (2017) [Optimal control of a membrane filtration system](https://doi.org/10.1016/j.ifacol.2017.08.1554), _IFAC-PapersOnLine_
- N. Kalboussi, J. Harmand, A. Rapaport, T. Bayen, F. Ellouze, N. Ben Amar (2018) [Optimal control of physical backwash strategy - towards the enhancement of membrane filtration process performance](https://www.sciencedirect.com/science/article/abs/pii/S0376738817319166), _Journal of Membrane Science_

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