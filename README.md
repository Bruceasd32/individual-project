# SCIG Wind-Turbine Hamiltonian-Passivity Controller and Operating-Range Analysis

## 1. Introduction

This repository contains the MATLAB code used to produce the simulation and
analysis results in the Year-3 individual project *"Methodology and
Operating-Range Analysis of a Hamiltonian-Passivity-Based Bounded Controller for
a Grid-Connected SCIG Wind Turbine"*.

The code does two things. First, it implements the full 14-state averaged model
of a squirrel-cage induction generator (SCIG) wind turbine connected to the grid
through a back-to-back AC/DC/AC converter, together with the parameter-free
Hamiltonian-passivity-based bounded controller reconstructed in Part I of the
report. Second, it implements the original contribution of Part II: the
*feasibility-envelope* analysis (Concept A) that computes the admissible
steady-state operating range of the four controlled quantities
`(ids, wr, id, Vdc)` and verifies it against the full nonlinear model.

## 2. Contextual overview (data flow)

```
            2-MW SCIG parameters (Table I) + hardware/standard limits
                                  |
        +-------------------------+--------------------------+
        |                         |                          |
        v                         v                          v
 scig_wt_nonlinear_        feasibility_envelope.m     feasibility_envelope_3d.m
   control_2.m             (Concept A, 2-D:           (Concept A, 3-D:
 (14-state closed-loop      sweep v, Vdc_ref;          sweep v, ids_ref, Vdc_ref)
  model + controller;       fixes ids_ref = 286 A)            |
  reference simulation)            |                          |
        |                          v                          v
        |                  fig_feas_envelope_*       fig_feas_3d_slices_*
        |                  fig_feas_margins_*        fig_feas_3d_volume_*
        |                  fig_feas_verify_*         feasibility_envelope_3d_results.mat
        |                  feasibility_envelope_results.mat
        v
 bounds_analysis.m  (bounds B1-B4 of the four controlled states)
        |
        v
 fig_bounds_*, fig_sim_errors_*, fig_theory_validation_*, fig_wind_profile_*
 bounds_analysis_results.txt / .mat
```

Every figure file is generated in two language variants: `*_en.png` (used in the
English report) and `*_zh.png` (used in the Chinese version). All four scripts
re-use the same physical parameter set (`setupParams`, defined locally inside
each script), so the scripts are self-contained and can be run independently in
any order.

## 3. Installation

- **MATLAB R2020a or later** (the code uses `exportgraphics`, introduced in
  R2020a). Developed and tested on a recent release under Windows.
- **Control System Toolbox** — required by `bounds_analysis.m`, which uses
  `ss`, `step`, and `lyap`. The three feasibility scripts and
  `scig_wt_nonlinear_control_2.m` use only base-MATLAB functions
  (`ode23tb`, `eig`, `exportgraphics`).
- No other dependencies. Clone or download the folder and add it to the MATLAB
  path (or simply `cd` into it).

```matlab
cd code
```

## 4. How to run

Each script is a standalone entry point. Run any of the following from the MATLAB
command window; figures and result files are written to the current folder.

```matlab
scig_wt_nonlinear_control_2   % closed-loop reference simulation of the 14-state model
feasibility_envelope          % 2-D feasibility envelope + boundary-corner verification
feasibility_envelope_3d       % 3-D feasibility envelope (v, ids_ref, Vdc_ref)
bounds_analysis               % bounds B1-B4 of the four controlled states
```

A typical full reproduction of the report figures is:

```matlab
feasibility_envelope;          % Part II, 2-D envelope figures
feasibility_envelope_3d;       % Part II, 3-D envelope figures
bounds_analysis;               % bounds-comparison figures
```

Run time is dominated by the parameter sweeps and the boundary ODE checks; each
script completes in roughly a minute on a desktop machine.

## 5. Technical details

- **Model.** 14 states: stator `dq` currents, rotor `dq` flux, rotor speed, grid
  `dq` currents, DC-link voltage, and six controller states. The averaged
  converter model uses loss-less duty-ratio relations. Parameters follow Table I
  of the source paper (2-MW SCIG).
- **Controller.** Parameter-free Hamiltonian-passivity-based bounded controller.
  The duty ratios are confined to the unit sphere
  (`m_ds^2 + m_qs^2 <= 1`, `m_dg^2 + m_qg^2 <= 1`), i.e. the linear-modulation
  region, by construction.
- **Feasibility envelope (Concept A).** For each point of a reference sweep the
  closed-form equilibrium is computed, every physical/hardware and
  controller-side (sphere) constraint is evaluated, and the cell is marked
  feasible / infeasible with the binding limit recorded. The feasible set is then
  projected onto the four controlled states to give the *state operating box*. A
  boundary corner is fed into the full nonlinear model and integrated to confirm
  convergence inside the envelope and divergence outside it.
- **Bounds (B1-B4).** B1 empirical bound from the full nonlinear ODE simulation;
  B2 peak step response of the linearised closed loop; B3 Lyapunov ultimate
  bound (`lyap`); B4 quasi-static reference-jump plus dynamic-peak estimate.
- **Integrator.** Stiff solver `ode23tb` with tightened tolerances near the
  envelope boundary.

## 6. Known issues and future improvements

- The feasibility envelope is computed on the **averaged** converter model and
  treats the duty-ratio relations as loss-less; switching ripple is not modelled.
- A **linear-iron** magnetic model is used; saturation is imposed only as a
  reference-level current ceiling, not as a dynamic effect.
- Only one **equilibrium family** is characterised (steady-state feasibility),
  not the full transient envelope.
- Grid voltage `Um` is held fixed; sweeping it (to study low-voltage ride-through)
  is left as future work.
- Figure auto-zoom occasionally needs manual axis adjustment for near-flat
  DC-link traces.
