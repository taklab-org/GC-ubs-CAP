# Codes of "Saddle-Type Blow-Up Solutions with Computer-Assisted Proofs: Validation and Extraction of Global Nature"

This repository contains the codes associated with the paper: "Saddle-Type Blow-Up Solutions with Computer-Assisted Proofs: Validation and Extraction of Global Nature" by J-P Lessard, K Matsue and A Takayasu. ([arXiv:2103.12390 [math.DS]](https://arxiv.org/abs/2103.12390))

**Abstract** In this paper, blow-up solutions of autonomous ordinary differential equations (ODEs) which are unstable under perturbations of initial points, referred to as saddle-type blow-up solutions, are studied. Combining dynamical systems machinery (e.g., compactifications, time-scale desingularizations of vector fields) with tools from computer-assisted proofs (e.g., rigorous integrators, the parameterization method for invariant manifolds), these blow-up solutions are obtained as trajectories on local stable manifolds of hyperbolic saddle equilibria at infinity. With the help of computer-assisted proofs, global trajectories on stable manifolds, inducing blow-up solutions, provide a global picture organized by global-in-time solutions and blow-up solutions simultaneously. Using the proposed methodology, intrinsic features of saddle-type blow-ups are observed: locally smooth dependence of blow-up times on initial points, level set distribution of blow-up times, and decomposition of the phase space playing a role as separatrixes among solutions, where the magnitude of initial points near those blow-ups does not matter for asymptotic behavior. Finally, singular behavior of blow-up times on initial points belonging to different family of blow-up solutions is addressed.

These codes require MATLAB with INTLAB - INTerval LABoratory (MATLAB toolbox for interval arithmetic) version 11 and the kv library (a C++ Library for rigorous numerics) version 0.4.48.

---

## Example 1

In `/Example1` folder, the codes of Section 5 is contained.

- `scrip_prove_manifold.m` provides rigorous enclosure of stable manifold and plots Figure 3.
- `script_compute_t_max.m` provides the rigorous bound of the blow-up time and plots Figure 4.
- `ex1_integrate_stable_manifold.cc` extends the stable manifold by integrating (5.7) in backward time and computes rigorous blow-up time. This code executes by compiling the following:
```
c++ -I../.. -O3 -DNDEBUG -DKV_FASTROUND ex1_integrate_stable_manifold.cc
```
- `script_plot3_blowup_times.m` plots Figure 5 and 6. It also provides the data of Table 1.


## Example 2

In `/Example2` folder, the codes of Section 6 is contained.

- `script_plot_manifolds.m` plots Figure 7.
- `script_get_distribution_of_blowup_time.m` plots Figure 8.
- `ex2_integrate_backward_from_**.cc` proves connecting orbits by integrating (6.3) in backward time, which connects ** (p0, p1, p2) to pb. For example, the connecting orbit from p0 to pb is proved by executing the following:
```
c++ -I../.. -O3 -DNDEBUG -DKV_FASTROUND ex2_integrate_backward_from_p0.cc
```
- `script_2d_stable_global_ex2.m` plots Figure 9.

## Example 3

In `/Example3` folder, the codes of Section 7 is contained.

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
