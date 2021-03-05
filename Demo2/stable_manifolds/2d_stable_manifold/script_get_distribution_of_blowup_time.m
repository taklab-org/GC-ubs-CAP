close all
clear 
clc

% scaling = [.05;.05]; N = 50; fp = 1;
scaling = [.07;.15]; N = 60; fp = 2;
 
a = zeros(3*((N+1)*(N+2)/2-3),1);
a = newton(a,scaling,fp);
% 
tail_terms(a,scaling,fp)
% 
plot_distribution_of_blowuptime(a,scaling,fp,1),hold on

nu = 1;
a = mid(a);
N = round((-3+sqrt(25+8/3*length(a)))/2,0); % a = (a0,a1,a2) with ai = (ai_{n})_{|n|=2}^N \in R^{(N+1)*(N+2)/2-3}

i0 = (1:(N+1)*(N+2)/2-3); i1 = i0+(N+1)*(N+2)/2-3; i2 = i1+(N+1)*(N+2)/2-3;
a0 = a(i0); a1 = a(i1); a2 = a(i2);

if fp == 1
    tx = [0.718092800211620;0.695947361719431;0];
    %lambda = [-1.031314539431529;-0.114370869183374];
    v1 = scaling(1)*[0.023249227542643;0.999241674152604;-0.03123379668521];
    v2 = scaling(2)*[0.664951647788332;-0.686110785140320;0.295112345756185];
end

if fp == 2
    tx = [0.998562893876901;-0.053592415248714;0];
    %lambda = [-1.994255706055617;-0.187090101867703];
    v1 = scaling(1)*[-0.994185185162860;0.107676008418592;-0.001301850117418];
    v2 = scaling(2)*[0.052670683572220;0.981388690288807;-0.184667370331782];
end

a0 = reshape_data_2d_manifold_vec2matrix(a0,tx(1),v1(1),v2(1));
a1 = reshape_data_2d_manifold_vec2matrix(a1,tx(2),v1(2),v2(2));
a2 = reshape_data_2d_manifold_vec2matrix(a2,tx(3),v1(3),v2(3));

% r = linspace(0,nu,100);
% t = linspace(0,2*pi,100) ;
% [R,T] = meshgrid(r,t);
% Theta1 = R.*cos(T); Theta2 = R.*sin(T);

% Draw slow manifold
theta1 = 0*linspace(-nu,nu,100); % =0 corresponds to Slow manifold
theta2 = linspace(-nu,nu,100); % Fast

[M1,M2] = size(theta1);

P0 = zeros(M1,M2);
P1 = zeros(M1,M2);
P2 = zeros(M1,M2);

for alpha1 = 0:N
    for alpha2 = 0:N
        P0 = P0 + a0(alpha1+1,alpha2+1)*(theta1.^alpha1).*(theta2.^alpha2);
        P1 = P1 + a1(alpha1+1,alpha2+1)*(theta1.^alpha1).*(theta2.^alpha2);
        P2 = P2 + a2(alpha1+1,alpha2+1)*(theta1.^alpha1).*(theta2.^alpha2);
    end
end

plot3(P0,P1,P2,'Linewidth',3)

% Draw fast manifold
theta1 = linspace(-nu,nu,100); % Slow
theta2 = 0*linspace(-nu,nu,100); % =0 corresponds to Fast manifold

[M1,M2] = size(theta1);

P0 = zeros(M1,M2);
P1 = zeros(M1,M2);
P2 = zeros(M1,M2);
% P3 = zeros(M1,M2);

for alpha1 = 0:N
    for alpha2 = 0:N
        P0 = P0 + a0(alpha1+1,alpha2+1)*(theta1.^alpha1).*(theta2.^alpha2);
        P1 = P1 + a1(alpha1+1,alpha2+1)*(theta1.^alpha1).*(theta2.^alpha2);
        P2 = P2 + a2(alpha1+1,alpha2+1)*(theta1.^alpha1).*(theta2.^alpha2);
    end
end

plot3(P0,P1,P2,'Linewidth',3)

colorbar
