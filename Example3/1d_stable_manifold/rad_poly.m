function [I,success] = rad_poly(a,par,nu,plus)

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
    
    load initial_data_ex2
    
    % Stable eigenvalue (includes the rigorous error bound)
    lambda = tlambda; %-0.187256681090725;
    
    % Fixed point (includes the rigorous error bound)
    % loads tx = [0.8861081289780320;0.6192579489210105] + error;
    
    % Stable eigenvector (+ rigorous error bound for the proof)
    v = scaling*tv;
    
end

I = [-1 -1];
success = 0 ;

[f1,f2] = f_stable_manifold_ext(a,par);

f1F = f1(1:N-1); f1_tail = f1(N:end); % size of f1_tail = 3*N
f2F = f2(1:N-1); f2_tail = f2(N:end); % size of f2_tail = 4*N

disp('Evaluating the jacobian with interval arithmetic ...')
Df = Df_stable_manifold(a,par);
A = inv(mid(Df));

%%%%%%%%%%
%%% Y0 %%%
%%%%%%%%%%

y0_F = A*[f1F;f2F];

Y0_1 = norma([y0_F(1:N-1);1./(lambda*(N+1:11*N)').*f1_tail],nu);
Y0_2 = norma([y0_F(N:2*N-2);1./(lambda*(N+1:13*N)').*f2_tail],nu);

Y0 = max([sup(Y0_1) sup(Y0_2)]);

disp(['Y0 = ',num2str(Y0)])

%%%%%%%%%%
%%% Z0 %%%
%%%%%%%%%%

iA = intval(A);
B = intval(eye(2*(N-1)) - iA*Df);

i1 = (1:N-1); i2 = (N:2*N-2);

B11 = B(i1,i1); B12 = B(i1,i2); B21 = B(i2,i1); B22 = B(i2,i2);

norm_B11 = mnorma(B11,nu);
norm_B12 = mnorma(B12,nu);
norm_B21 = mnorma(B21,nu);
norm_B22 = mnorma(B22,nu);

Z0_1 = norm_B11 + norm_B12;
Z0_2 = norm_B21 + norm_B22;

Z0 = max([sup(Z0_1) sup(Z0_2)]);

disp(['Z0 = ',num2str(Z0)])

%%%%%%%%%%
%%% Z1 %%%
%%%%%%%%%%

a1 = [tx(1);v(1);a1]; a2 = [tx(2);v(2);a2];

a1_2 = cauchy_ext([a1 a1],12*N+1);
a1_4 = cauchy_ext([a1 a1 a1 a1],12*N+1);
a1_5 = cauchy_ext([a1 a1 a1 a1 a1],12*N+1);
a1_6 = cauchy_ext([a1 a1 a1 a1 a1 a1],12*N+1);
a1_8 = cauchy_ext([a1 a1 a1 a1 a1 a1 a1 a1],12*N+1);
a1_10 = cauchy_ext([a1 a1 a1 a1 a1 a1 a1 a1 a1 a1],12*N+1);
a1_12 = cauchy_ext([a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1],12*N+1);
a1a2 = cauchy_ext([a1 a2],12*N+1);
a1a2_5 = cauchy_ext([a1 a2 a2 a2 a2 a2],12*N+1);
a1_2a2 = cauchy_ext([a1 a1 a2],12*N+1);
a1_2a2_4 = cauchy_ext([a1 a1 a2 a2 a2 a2],12*N+1);
a1_2a2_2 = cauchy_ext([a1 a1 a2 a2],12*N+1);
a1_4a2_2 = cauchy_ext([a1 a1 a1 a1 a2 a2],12*N+1);
a1_4a2 = cauchy_ext([a1 a1 a1 a1 a2],12*N+1);
a1_4a2_4 = cauchy_ext([a1 a1 a1 a1 a2 a2 a2 a2],12*N+1);
a1_6a2_2 = cauchy_ext([a1 a1 a1 a1 a1 a1 a2 a2],12*N+1);
a1a2_2 = cauchy_ext([a1 a2 a2],12*N+1);
a1a2_3 = cauchy_ext([a1 a2 a2 a2],12*N+1);
a1_3a2 = cauchy_ext([a1 a1 a1 a2],12*N+1);
a1_5a2 = cauchy_ext([a1 a1 a1 a1 a1 a2],12*N+1);
a1_9a2 = cauchy_ext([a1 a1 a1 a1 a1 a1 a1 a1 a1 a2],12*N+1);
a1_5a2_3 = cauchy_ext([a1 a1 a1 a1 a1 a2 a2 a2],12*N+1);
a1_8a2_2 = cauchy_ext([a1 a1 a1 a1 a1 a1 a1 a1 a2 a2],12*N+1);
a2_2 = cauchy_ext([a2 a2],12*N+1);
a2_4 = cauchy_ext([a2 a2 a2 a2],12*N+1);
a2_6 = cauchy_ext([a2 a2 a2 a2 a2 a2],12*N+1);

one = [1;zeros(12*N,1)];

c11 = (1/2)*[a1;zeros(11*N,1)] - (3/2)*a1_5 + a1a2 + (1/3)*a1_3a2 - 6*a1_5a2 + 5*a1_9a2 + (3/2)*a1a2_2 - 2*a1a2_3 + 6*a1_5a2_3 + a1a2_5;
c12 = -(1/4)*one + (1/2)*a1_2 + (1/12)*a1_4 - a1_6 + (1/2)*a1_10 + (3/2)*a1_2a2 - (9/4)*a2_2 - 3*a1_2a2_2 + 3*a1_6a2_2 + (5/2)*a1_2a2_4;
c21 = -(1/4)*one + (1/4)*a1_2 - (5/4)*a1_4 + (7/4)*a1_6 + (45/4)*a1_8 - (39/4)*a1_12 - 10*a1_4a2 + (3/4)*a2_2 + (23/4)*a1_2a2_2 + (5/2)*a1_4a2_2 - (45/4)*a1_8a2_2 - (3/4)*a2_4 - (5/4)*a1_4a2_4 + (1/4)*a2_6;
c22 = -2*a1_5 + (3/2)*a1a2 + (23/6)*a1_3a2 + a1_5a2 - (5/2)*a1_9a2 - 3*a1a2_3 - a1_5a2_3 + (3/2)*a1a2_5;

Z1_1 = (1/(abs(lambda)*(N+1)))*(norma(c11,nu)+norma(c12,nu));
Z1_2 = (1/(abs(lambda)*(N+1)))*(norma(c21,nu)+norma(c22,nu));

Z1 = max([sup(Z1_1) sup(Z1_2)]);

disp(['Z1 = ',num2str(Z1)])

%%%%%%%%%%
%%% Z2 %%%
%%%%%%%%%%

A11 = A(i1,i1); A12 = A(i1,i2); A21 = A(i2,i1); A22 = A(i2,i2);

r_star = 1e-6;

h = infsup(-r_star,r_star)*nu.^-(0:N)';

b1 = a1 + h;
b2 = a2 + h;

b1_2 = cauchy([b1 b1]);
b1_3 = cauchy([b1 b1 b1]);
b1_4 = cauchy([b1 b1 b1 b1]);
b1_5 = cauchy([b1 b1 b1 b1 b1]);
b1_7 = cauchy([b1 b1 b1 b1 b1 b1 b1]);
b1_9 = cauchy([b1 b1 b1 b1 b1 b1 b1 b1 b1]);
b1_11 = cauchy([b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1]);
b1_2b2 = cauchy([b1 b1 b2]);
b1_3b2 = cauchy([b1 b1 b1 b2]);
b1_3b2_2 = cauchy([b1 b1 b1 b2 b2]);
b1_2b2_3 = cauchy([b1 b1 b2 b2 b2]);
b1_6b2 = cauchy([b1 b1 b1 b1 b1 b1 b2]);
b1b2_2 = cauchy([b1 b2 b2]);
b1_4b2 = cauchy([b1 b1 b1 b1 b2]);
b1_3b2_4 = cauchy([b1 b1 b1 b2 b2 b2 b2]);
b1_5b2_2 = cauchy([b1 b1 b1 b1 b1 b2 b2]);
b1_7b2_2 = cauchy([b1 b1 b1 b1 b1 b1 b1 b2 b2]);
b1_4b2_3 = cauchy([b1 b1 b1 b1 b2 b2 b2]);
b1_8b2 = cauchy([b1 b1 b1 b1 b1 b1 b1 b1 b2]);
b1b2 = cauchy([b1 b2]);
b1b2_4 = cauchy([b1 b2 b2 b2 b2]);
b2_2 = cauchy([b2 b2]);
b2_3 = cauchy([b2 b2 b2]);
b2_5 = cauchy([b2 b2 b2 b2 b2]);

nb1 = norma(a1,nu) + r_star;
nb2 = norma(a2,nu) + r_star;

f1_x1x1 = 1/2*h - (15*b1_4)/2 + b2 + b1_2b2 - 30*b1_4b2 + 45*b1_8b2 + (3*b2_2)/2 - 2*b2_3 + 30*b1_4b2_3 + b2_5;
f1_x1x2 = b1 + b1_3/3 - 6*b1_5 + 5*b1_9 + 3*b1b2 - 6*b1b2_2 + 18*b1_5b2_2 + 5*b1b2_4;
f1_x2x2 = (3*b1_2)/2 - (9*b2)/2 - 6*b1_2b2 + 6*b1_6b2 + 10*b1_2b2_3;
f2_x1x1 = b1/2 - 5*b1_3 + (21*b1_5)/2 + 90*b1_7 - 117*b1_11 - 40*b1_3b2 + (23*b1b2_2)/2 + 10*b1_3b2_2 - 90*b1_7b2_2 - 5*b1_3b2_4;
f2_x1x2 = -10*b1_4 + (3*b2)/2 + (23*b1_2b2)/2 + 5*b1_4b2 - (45*b1_8b2)/2 - 3*b2_3 - 5*b1_4b2_3 + (3*b2_5)/2;
f2_x2x2 = (3*b1)/2 + (23*b1_3)/6 + b1_5 - (5*b1_9)/2 - 9*b1b2_2 - 3*b1_5b2_2 + (15*b1b2_4)/2;

nf1_x1x1 = 1/2*r_star + (15*nb1^4)/2 + nb2 + nb1^2*nb2 + 30*nb1^4*nb2 + 45*nb1^8*nb2 + (3*nb2^2)/2 + 2*nb2^3 + 30*nb1^4*nb2^3 + nb2^5;
nf1_x1x2 = nb1 + nb1^3/3 + 6*nb1^5 + 5*nb1^9 + 3*nb1*nb2 + 6*nb1*nb2^2 + 18*nb1^5*nb2^2 + 5*nb1*nb2^4;
nf1_x2x2 = (3*nb1^2)/2 + (9*nb2)/2 + 6*nb1^2*nb2 + 6*nb1^6*nb2 + 10*nb1^2*nb2^3;
nf2_x1x1 = nb1/2 + 5*nb1^3 + (21*nb1^5)/2 + 90*nb1^7 + 117*nb1^11 + 40*nb1^3*nb2 + (23*nb1*nb2^2)/2 + 10*nb1^3*nb2^2 + 90*nb1^7*nb2^2 + 5*nb1^3*nb2^4;
nf2_x1x2 = 10*nb1^4 + (3*nb2)/2 + (23*nb1^2*nb2)/2 + 5*nb1^4*nb2 + (45*nb1^8*nb2)/2 + 3*nb2^3 + 5*nb1^4*nb2^3 + (3*nb2^5)/2;
nf2_x2x2 = (3*nb1)/2 + (23*nb1^3)/6 + nb1^5 + (5*nb1^9)/2 + 9*nb1*nb2^2 + 3*nb1^5*nb2^2 + (15*nb1*nb2^4)/2;

Z2_1_F = norma(A11*f1_x1x1(3:N+1) + A12*f2_x1x1(3:N+1),nu) ...
     + 2*norma(A11*f1_x1x2(3:N+1) + A12*f2_x1x2(3:N+1),nu) ...
       + norma(A11*f1_x2x2(3:N+1) + A12*f2_x2x2(3:N+1),nu);

Z2_1_tail = (1/(lambda*(N+1)))*(nf1_x1x1 + 2*nf1_x1x2 + nf1_x2x2);

Z2_1 = Z2_1_F + Z2_1_tail;

Z2_2_F = norma(A21*f1_x1x1(3:N+1) + A22*f2_x1x1(3:N+1),nu) ...
     + 2*norma(A21*f1_x1x2(3:N+1) + A22*f2_x1x2(3:N+1),nu) ...
       + norma(A21*f1_x2x2(3:N+1) + A22*f2_x2x2(3:N+1),nu);

Z2_2_tail = (1/(lambda*(N+1)))*(nf2_x1x1 + 2*nf2_x1x2 + nf2_x2x2);

Z2_2 = Z2_2_F + Z2_2_tail;

Z2 = max([sup(Z2_1) sup(Z2_2)]);

disp(['Z2 = ',num2str(Z2)])

%%%%%%%%%%%%%%%%%%%%%%%%
%%% Radii polynomial %%%
%%%%%%%%%%%%%%%%%%%%%%%%

Y0 = intval(Y0); Z0 = intval(Z0); Z1 = intval(Z1); Z2 = intval(Z2);

if inf(1-Z0-Z1)>0 
  if inf((1-Z0-Z1)^2-4*Y0*Z2) > 0  
    rmin=sup(((1-Z0-Z1) - sqrt((1-Z0-Z1)^2-4*Y0*Z2))/(2*Z2));
    rmax=inf(((1-Z0-Z1) + sqrt((1-Z0-Z1)^2-4*Y0*Z2))/(2*Z2));
    if rmin<rmax 
      success=1;
      I=[rmin rmax];
      disp('success')
      plot_manifold(mid([a1;a2]),mid(tx),nu,plus)

    else
      disp('failure: rmin > rmax')
    end
  else
    disp('failure: discriminant is negative')  
  end
else
    disp('failure: linear term is positive')
end


end

