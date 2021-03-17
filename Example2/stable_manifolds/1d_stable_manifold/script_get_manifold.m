% close all
% clear 
% clc

scaling = .12; 

N = 160; 

a0 = zeros(N-1,1);  a1 = zeros(N-1,1); a2 = zeros(N-1,1);

a = [a0;a1;a2];
a = newton(a,scaling);

tx = [0.933378938375587;0.358892403648868;0]; % Fixed point
v = scaling*[0.743908505928788;0.668260964215875;0.005236268906990];
N = (length(a)+3)/3; % a = (a0,a1,a2) with ai = (ai_n)_{n=2}^N \in R^{N-1}
a0 = a(1:N-1); a1 = a(N:2*N-2); a2 = a(2*N-1:3*N-3);
a0 = [tx(1);v(1);a0]; a1 = [tx(2);v(2);a1]; a2 = [tx(3);v(3);a2];
disp([a0 a1 a2])
plot_manifold([a0;a1;a2],tx,1);

nu = 1;
a = intval(a);
[I,success] = rad_poly_ex3_1d(a,scaling,nu);

rmin = I(1);
a = intval([a0 a1 a2]);
P_at_1 = sum(a,1); error_at_1 = rmin + rad(P_at_1);
P_at_1 = midrad(P_at_1.mid,error_at_1);

P_at_minus_1 = sum(a.*(repmat((-1).^(0:N)',1,3)),1);
error_at_minus_1 = rmin + rad(P_at_minus_1);
P_at_minus_1 = midrad(P_at_minus_1.mid,error_at_minus_1);

% save end_points_1d_stable_manifold P_at_minus_1