close all
clear 
clc

% scaling = [.05;.05]; N = 50; fp = 1;
scaling = [.07;.15]; N = 60; fp = 2;
 
a = zeros(3*((N+1)*(N+2)/2-3),1);
a = newton(a,scaling,fp);

tail_terms(a,scaling,fp)

plot_2d_manifold(a,scaling,fp,1)

a = intval(a);
[I,success] = rad_poly_ex3_2d(a,scaling,fp);

r0 = I(1);
disp(r0)


if fp == 1
  p1 = point_on_2d_manifold(a,scaling,fp,r0,-0.2,1);
  p2 = point_on_2d_manifold(a,scaling,fp,r0,-1,1);
  p3 = point_on_2d_manifold(a,scaling,fp,r0,-1,2/3);
  p4 = point_on_2d_manifold(a,scaling,fp,r0,-1,1/3);
  p5 = point_on_2d_manifold(a,scaling,fp,r0,-1,0);
  p6 = point_on_2d_manifold(a,scaling,fp,r0,-1,-1/3);
  p7 = point_on_2d_manifold(a,scaling,fp,r0,-1,-2/3);
  p8 = point_on_2d_manifold(a,scaling,fp,r0,-1,-1);
  p9 = point_on_2d_manifold(a,scaling,fp,r0,-0.2,-1);
  pp = [p1 p2 p3 p4 p5 p6 p7 p8 p9];
  save p1.mat pp
end

if fp == 2
  p1 = point_on_2d_manifold(a,scaling,fp,r0,0.2,1);
  p2 = point_on_2d_manifold(a,scaling,fp,r0,1,1);
  p3 = point_on_2d_manifold(a,scaling,fp,r0,1,2/3);
  p4 = point_on_2d_manifold(a,scaling,fp,r0,1,1/3);
  p5 = point_on_2d_manifold(a,scaling,fp,r0,1,0);
  p6 = point_on_2d_manifold(a,scaling,fp,r0,1,-1/3);
  p7 = point_on_2d_manifold(a,scaling,fp,r0,1,-2/3);
  p8 = point_on_2d_manifold(a,scaling,fp,r0,1,-1);
  p9 = point_on_2d_manifold(a,scaling,fp,r0,0.2,-1);
  pp = [p1 p2 p3 p4 p5 p6 p7 p8 p9];
  save p2.mat pp
end

% theta = [1/2;1/2];
% 
% t_max = rigorous_computation_tmax_ex3_2d(a,scaling,fp,theta,r0);
