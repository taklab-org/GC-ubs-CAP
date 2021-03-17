function f = f_stable_manifold(a,par)

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

n = (2:N)';
% 
% one = [1;zeros(N,1)];
% a1_3 = cauchy([a1 a1 a1]);
% 
% % p4 = x1^4 + x2^2;
% p4 = cauchy([a1 a1 a1 a1]) + cauchy([a2 a2]);
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
% phi1 = cauchy([r F]) - cauchy([a1 G]) ;
% phi2 = cauchy([q F]) - 2*cauchy([a2 G]);

a1_2 = cauchy([a1 a1]);
a1_3 = cauchy([a1 a1 a1]);
a1_5 = cauchy([a1 a1 a1 a1 a1]);
a1_6 = cauchy([a1 a1 a1 a1 a1 a1]);
a1_7 = cauchy([a1 a1 a1 a1 a1 a1 a1]);
a1_9 = cauchy([a1 a1 a1 a1 a1 a1 a1 a1 a1]);
a1_13 = cauchy([a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1]);
a1a2_2 = cauchy([a1 a2 a2]);
a1a2_4 = cauchy([a1 a2 a2 a2 a2]);
a1a2_6 = cauchy([a1 a2 a2 a2 a2 a2 a2]);
a1_2a2 = cauchy([a1 a1 a2]);
a1_2a2_2 = cauchy([a1 a1 a2 a2]);
a1_2a2_3 = cauchy([a1 a1 a2 a2 a2]);
a1_2a2_5 = cauchy([a1 a1 a2 a2 a2 a2 a2]);
a1_3a2_2 = cauchy([a1 a1 a1 a2 a2]);
a1_5a2_2 = cauchy([a1 a1 a1 a1 a1 a2 a2]);
a1_4a2 = cauchy([a1 a1 a1 a1 a2]);
a1_5a2 = cauchy([a1 a1 a1 a1 a1 a2]);
a1_5a2_4 = cauchy([a1 a1 a1 a1 a1 a2 a2 a2 a2]);
a1_6a2 = cauchy([a1 a1 a1 a1 a1 a1 a2]);
a1_6a2_3 = cauchy([a1 a1 a1 a1 a1 a1 a2 a2 a2]);
a1_9a2_2 = cauchy([a1 a1 a1 a1 a1 a1 a1 a1 a1 a2 a2]);
a1_10a2 = cauchy([a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a2]);
a2_3 = cauchy([a2 a2 a2]);

phi1 = (1/4)*a1_2 - (1/4)*a1_6 - (1/4)*a2 + (1/2)*a1_2a2 + (1/12)*a1_4a2 - a1_6a2 + (1/2)*a1_10a2 + (3/4)*a1_2a2_2 - (3/4)*a2_3 - a1_2a2_3 + a1_6a2_3 + (1/2)*a1_2a2_5;
phi2 = -(1/4)*a1 + (1/12)*a1_3 - (1/4)*a1_5 + (1/4)*a1_7 + (5/4)*a1_9 - (3/4)*a1_13 -  2*a1_5a2 + (3/4)*a1a2_2 + (23/12)*a1_3a2_2 + (1/2)*a1_5a2_2 - (5/4)*a1_9a2_2 - (3/4)*a1a2_4 - (1/4)*a1_5a2_4 + (1/4)*a1a2_6;

f = [lambda*n.*a1(3:end) - mid(phi1(3:end)) ; lambda*n.*a2(3:end) - mid(phi2(3:end))];

end

