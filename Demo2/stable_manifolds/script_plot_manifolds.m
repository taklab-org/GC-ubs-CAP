addpath('2d_stable_manifold')
close all
clear 
clc

% 2d_stable_manifold_p1
load 2d_manifold_data_p1.mat

plot_2d_manifold(a,scaling,fp,nu)

hold on

% 2d_stable_manifold_p2
load 2d_manifold_data_p2.mat

plot_2d_manifold(a,scaling,fp,nu)

rmpath('2d_stable_manifold')

addpath('1d_stable_manifold')

% 1d_stable_manifold_p0
load 1d_manifold_data_p0.mat

plot_manifold(a,tx,nu)

rmpath('2d_stable_manifold')

view([-103,30])
