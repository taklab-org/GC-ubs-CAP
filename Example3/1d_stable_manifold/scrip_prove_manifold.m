close all
clear
figure 
hold on

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Proof of the stable manifold of 0')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

disp('      ')

load p_0

nu = 1;
[I] = rad_poly(a,par,nu,0);

disp(['||P(theta) - P^(N)(theta)|| < = r_min = ',num2str(I(1))])

% return

%%%%%%%%%%%%%
%%%%%%%%%%%%%

disp('      ')

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Proof of the stable manifold of p_{infty,s}')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

disp('      ')

load p_infty_s

nu = 1;
[I] = rad_poly(a,par,nu,0);
disp(['||P(theta) - P^(N)(theta)|| < = r_min = ',num2str(I(1))])

tx = [0.7328506362011802;0.5370700549804747]; % p_b (the source) 
plot(tx(1),tx(2),'*','Linewidth',5,'color',[0 0 0])

% return

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Computing t_max %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
figure
r0 = I(1);
load initial_data_ex2
scaling = par(1);
lambda = tlambda;
v = scaling*tv;
N = (length(a)+2)/2; % a = (a1,a2) with ai = (ai_n)_{n=2}^N \in R^{N-1}
a1 = a(1:N-1); a2 = a(N:2*N-2);
a1 = [tx(1);v(1);a1]; a2 = [tx(2);v(2);a2];
a = [a1;a2];

t_max = rigorous_computation_tmax_ex2(a,lambda,r0,1);

% disp(['||P(theta) - P^(N)(theta)|| < = r_min = ',num2str(I(1))])


%
figure
plot_blowup_times(t_max(end-1)) % Plot Figure 13(a)