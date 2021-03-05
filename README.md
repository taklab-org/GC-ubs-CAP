# Codes of "A geometric characterization of unstable blow-up solutions with computer-assisted proof"

This repository contains the codes associated with the paper: "A geometric characterization of unstable blow-up solutions with computer-assisted proof" by J-P Lessard, K Matsue and A Takayasu. (arXiv:)

**Abstract** In this paper, blow-up solutions of autonomous ordinary differential equations (ODEs) which are unstable under perturbations of initial conditions are studied. Combining dynamical systems machinery (e.g. phase space compactifications, time-scale desingularizations of vector fields) with tools from computer-assisted proofs (e.g. rigorous integrators, parameterization method for invariant manifolds) these {\em unstable} blow-up solutions are obtained as trajectories on stable manifolds of hyperbolic (saddle) equilibria at infinity. In this process, important features are obtained: smooth dependence of blow-up times on initial conditions near blow-up, level set distribution of blow-up times, singular behavior of blow-up times on unstable blow-up solutions, organizing the phase space via separatrices (stable manifolds). In particular, we show that unstable blow-up solutions and bounded solutions connected by those blow-up solutions can separate initial conditions into two regions where solution trajectories are either globally bounded or blow-up, no matter how the large initial conditions are.

These codes require MATLAB with INTLAB - INTerval LABoratory (MATLAB toolbox for interval arithmetic) version 11 and the kv library (a C++ Library for rigorous numerics) version 0.4.48.

---

## Demonstration 1

In `/Demo1` folder, the codes of Section 5 is contained.

- `scrip_prove_manifold.m` provides rigorous enclosure of stable manifold and plots Figure 3.
- `script_compute_t_max.m` provides the rigorous bound of the blow-up time and plots Figure 4.
- `ex1_integrate_stable_manifold.cc` extends the stable manifold by integrating (5.7) in backward time and computes rigorous blow-up time. This code executes by compiling the following:
```
c++ -I../.. -O3 -DNDEBUG -DKV_FASTROUND ex1_integrate_stable_manifold.cc
```
- `script_plot3_blowup_times.m` plots Figure 5 and 6. It also provides the data of Table 1.


## Demonstration 2

In `/Demo2` folder, the codes of Section 6 is contained.

- `script_plot_manifolds.m` plots Figure 7.
- `script_get_distribution_of_blowup_time.m` plots Figure 8.
- `ex2_integrate_backward_from_**.cc` proves connecting orbits by integrating (6.3) in backward time, which connects ** (p0, p1, p2) to pb. For example, the connecting orbit from p0 to pb is proved by executing the following:
```
c++ -I../.. -O3 -DNDEBUG -DKV_FASTROUND ex2_integrate_backward_from_p0.cc
```
- `script_2d_stable_global_ex2.m` plots Figure 9.

## Demonstration 3

In `/Demo3` folder, the codes of Section 7 is contained.

- `scrip_prove_manifold.m` provides rigorous enclosure of 1D stable manifolds (attached to p_0 and p_{infty,s} equilibria) and plots Figure 11.
- `script_plot_fig12.m` plots Figure 12.
- `ex3_integrate_forward_with_blowup_time.cc` provides data of blow-up times, which gives the data of Figure 13.
- `script_plot_blowup_times.m` plots Figure 13.
- `script_plot_separatorix.m` plots Figure 14 and 15.
- Related to Figure 14, the following codes provides 6 computer-assisted proofs of connecting orbits, global solutions, and blow-up solutions
  - `ex3_integrate_backward_from_p0.cc`: connecting orbit from p_0 to p_b^+
  - `ex3_integrate_backward_from_pinfs.cc`: connecting orbit from p_{infty,s} to p_b^+
  - `ex3_integrate_forward_global_solution_near_p0.cc`: global solution converges to p_b^-
  - `ex3_integrate_forward_global_solution_near_pinfs.cc`: global solution converges to p_b^-
  - `ex3_integrate_forward_with_blowup_time_near_p0.cc`: blow-up solution near p_0
  - `ex3_integrate_forward_with_blowup_time_near_pinfs.cc`: blow-up solution near p_{infty,s}



Copyright (C) 2021 J-P Lessard, K Matsue and A Takayasu.
