function Df = Df_stable_manifold(a,par)

N = (length(a)+2)/2; % a = (a1,a2) with ai = (ai_n)_{n=2}^N \in R^{N-1}

a1 = a(1:N-1); a2 = a(N:2*N-2);

scaling = par(1); fp = par(2);


if fp == 0
    
    %%%%%%%%%%%%%%%
    %%%% (0,0) %%%%
    %%%%%%%%%%%%%%%
    
    lambda = -1/4; % Stable eigenvalue
    tx = [0;0]; % Fixed point
    v = scaling*[1;1]; % Stable eigenvector

end

if fp == 1
   
    %%%%%%%%%%%%%%%%%%%%%
    %%%% p_{infty,s} %%%%
    %%%%%%%%%%%%%%%%%%%%%
    
    % Stable eigenvalue (+ rigorous error bound for the proof)
    lambda = -0.187256681090725;
    
    % Fixed point (+ rigorous error bound for the proof)
    tx = [0.8861081289780320;0.6192579489210105];
    
    % Stable eigenvector (+ rigorous error bound for the proof)
    v = scaling*[-0.581870201791746;-0.813281666009282];
    
end

if fp == 2
   
    %%%%%%%%%%%%%
    %%%% p_b %%%%
    %%%%%%%%%%%%%
    
    % Unstable eigenvalue (+ rigorous error bound for the proof)
    lambda = 0.211448482247263;
    
    % Fixed point (+ rigorous error bound for the proof)
    tx = [0.732850636201180;0.537070054980475];
    
    % Stable eigenvector (+ rigorous error bound for the proof)
    v = scaling*[-0.981046226514401;-0.193773840963773];
    
end

a1 = [tx(1);v(1);a1]; a2 = [tx(2);v(2);a2];

% one = [1;zeros(N,1)];
% a1_3 = cauchy([a1 a1 a1]);
% a1_4 = cauchy([a1 a1 a1 a1]);
% a2_2 = cauchy([a2 a2]);
% 
% % p4 = x1^4 + x2^2;
% p4 = a1_4 + a2_2;
% 
% % q = x1^3/3 - (1 - p4)^2*x1;
% q = (1/3)*a1_3 - cauchy([one-p4 one-p4 a1]);
% 
% % r = x1^2 - x2;
% r = cauchy([a1 a1]) - a2;
% 
% % F = (1/4)*(1 + 3*p4);
% F = (1/4)*(one + 3*p4);
% 
% % G = x1^3*r + (x2/2)*q;
% G = cauchy([a1_3 r]) + (1/2)*cauchy([a2 q]);
% 
% % F_x1 = 3*x1^3;
% F_x1 = 3*a1_3;
% 
% % F_x2 = (3/2)*x2;
% F_x2 = (3/2)*a2;
%  
% % r_x1 = 2*x1; 
% r_x1 = 2*a1;
% 
% % r_x2 = -1;
% r_x2 = -one;
% 
% % q_x1 = x1^2 - (1-p4)*(7*x1^4-x2^2+1);
% q_x1 = cauchy([a1 a1]) - cauchy([one-p4 7*a1_4-a2_2+one]);
% 
% % q_x2 = 4*x1*x2*(1-p4);
% q_x2 = 4*cauchy([a1 a2 one-p4]);
% 
% % G_x1 = 3*x1^2*r + x1^3*r_x1 + (x2/2)*q_x1;
% G_x1 = 3*cauchy([a1 a1 r]) + cauchy([a1_3 r_x1]) + (1/2)*cauchy([a2 q_x1]);
% 
% % G_x2 = x1^3*r_x2 + (1/2)*q + (x2/2)*q_x2;
% G_x2 = cauchy([a1_3 r_x2]) + (1/2)*q + (1/2)*cauchy([a2 q_x2]);
% 
% c11 = cauchy([r_x1 F]) + cauchy([r F_x1]) - G - cauchy([a1 G_x1]);
% c12 = cauchy([r_x2 F]) + cauchy([r F_x2]) - cauchy([a1 G_x2]);
% c21 = cauchy([q_x1 F]) + cauchy([q F_x1]) - 2*cauchy([a2 G_x1]);
% c22 = cauchy([q_x2 F]) + cauchy([q F_x2]) - 2*G - 2*cauchy([a2 G_x2]);

%%%%%%%%%%%%%%%%%
%%% Version 2 %%%
%%%%%%%%%%%%%%%%%

one = [1;zeros(N,1)];

a1_2 = cauchy([a1 a1]);
a1_4 = cauchy([a1 a1 a1 a1]);
a1_5 = cauchy([a1 a1 a1 a1 a1]);
a1_6 = cauchy([a1 a1 a1 a1 a1 a1]);
a1_8 = cauchy([a1 a1 a1 a1 a1 a1 a1 a1]);
a1_10 = cauchy([a1 a1 a1 a1 a1 a1 a1 a1 a1 a1]);
a1_12 = cauchy([a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1]);
a1a2 = cauchy([a1 a2]);
a1a2_5 = cauchy([a1 a2 a2 a2 a2 a2]);
a1_2a2 = cauchy([a1 a1 a2]);
a1_2a2_4 = cauchy([a1 a1 a2 a2 a2 a2]);
a1_2a2_2 = cauchy([a1 a1 a2 a2]);
a1_4a2_2 = cauchy([a1 a1 a1 a1 a2 a2]);
a1_4a2 = cauchy([a1 a1 a1 a1 a2]);
a1_4a2_4 = cauchy([a1 a1 a1 a1 a2 a2 a2 a2]);
a1_6a2_2 = cauchy([a1 a1 a1 a1 a1 a1 a2 a2]);
a1a2_2 = cauchy([a1 a2 a2]);
a1a2_3 = cauchy([a1 a2 a2 a2]);
a1_3a2 = cauchy([a1 a1 a1 a2]);
a1_5a2 = cauchy([a1 a1 a1 a1 a1 a2]);
a1_9a2 = cauchy([a1 a1 a1 a1 a1 a1 a1 a1 a1 a2]);
a1_5a2_3 = cauchy([a1 a1 a1 a1 a1 a2 a2 a2]);
a1_8a2_2 = cauchy([a1 a1 a1 a1 a1 a1 a1 a1 a2 a2]);
a2_2 = cauchy([a2 a2]);
a2_4 = cauchy([a2 a2 a2 a2]);
a2_6 = cauchy([a2 a2 a2 a2 a2 a2]);

c11 = (1/2)*a1 - (3/2)*a1_5 + a1a2 + (1/3)*a1_3a2 - 6*a1_5a2 + 5*a1_9a2 + (3/2)*a1a2_2 - 2*a1a2_3 + 6*a1_5a2_3 + a1a2_5;
c12 = -(1/4)*one + (1/2)*a1_2 + (1/12)*a1_4 - a1_6 + (1/2)*a1_10 + (3/2)*a1_2a2 - (9/4)*a2_2 - 3*a1_2a2_2 + 3*a1_6a2_2 + (5/2)*a1_2a2_4;
c21 = -(1/4)*one + (1/4)*a1_2 - (5/4)*a1_4 + (7/4)*a1_6 + (45/4)*a1_8 - (39/4)*a1_12 - 10*a1_4a2 + (3/4)*a2_2 + (23/4)*a1_2a2_2 + (5/2)*a1_4a2_2 - (45/4)*a1_8a2_2 - (3/4)*a2_4 - (5/4)*a1_4a2_4 + (1/4)*a2_6;
c22 = -2*a1_5 + (3/2)*a1a2 + (23/6)*a1_3a2 + a1_5a2 - (5/2)*a1_9a2 - 3*a1a2_3 - a1_5a2_3 + (3/2)*a1a2_5;

dphi_11 = intval(zeros(N-1));
dphi_12 = intval(zeros(N-1));
dphi_21 = intval(zeros(N-1));
dphi_22 = intval(zeros(N-1));

for n = 2:N
    j = (2:n);
    k = n-j+1;
    dphi_11(n-1,j-1) = c11(k);
    dphi_12(n-1,j-1) = c12(k);
    dphi_21(n-1,j-1) = c21(k);
    dphi_22(n-1,j-1) = c22(k);
end

n = (2:N)';

Df_11 = lambda*diag(n) - dphi_11 ;
Df_12 = - dphi_12 ;
Df_21 = - dphi_21 ;
Df_22 = lambda*diag(n) - dphi_22 ;

Df = [[Df_11 Df_12];[Df_21 Df_22]];

end

