function [convolution] = cauchy_ext(a,n)

%%% n : size of the output we want

% Inputs : (p vectors of Taylor coefficients - a1, ..., ap) 
%  ai = (ai_0,ai_1,...,ai_{M-1}) \in R^M, for i=1,...,p

% This code computes a rigorous enclosure of a Cauchy product 
% (a1*...*ap)_k = sum_{k1+...+kp = k, ki >= 0} a1_k1 ... ap_kp for k = 0,...,p*(M-1)

[M,p] = size(a); % p : power of the nonlinearity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% First enclosure: Rigorous estimates based on the FFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We make sure that the inputs are powers of 2 %

exponent = ceil( log(2*p*M-1)/log(2) );
M1 = (2^exponent-2*p*M)/2; % Hence 2*p*M+2*M1 is a power of 2 %

a_ext = intval(zeros(2*M-1,p));
ta = intval(zeros(2*p*M+2*M1,p));
tu = intval(zeros(2*p*M+2*M1,p));
tb = ones(2*p*M+2*M1,1);

for i=1:p
    a_ext(:,i) = [zeros(M-1,1);a(:,i)];
    ta(:,i) = [intval(zeros((p-1)*M+M1+1,1));intval(a_ext(:,i));intval(zeros((p-1)*M+M1,1))];
    tu(:,i) = verifyfft(ifftshift(ta(:,i)),-1);
    tb = (2*p*M+2*M1)*tb.*tu(:,i);
end

F = fftshift(verifyfft(tb,1));

convolution1 = real(F(M1+p*M+1:M1+p*(2*M-1)+1))/(2*p*M+2*M1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Second enclosure: Rigorous estimates based on Banach Algebras
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho = zeros(p,1);

for i=1:p
    rho(i) = exp_decay_a_least_square(mid(a(:,i)));
end

rho_min = intval(min(rho));
omega = [intval(1) rho_min.^abs(1:M-1)]';

ia = intval(a);
norm_product = 1;

for i=1:p
    norm_product = norm_product*sum(abs(ia(:,i)).*omega);
end

k_p = (1:p*(M-1));
omega = [intval(1) rho_min.^k_p]';

convolution2 = sup(norm_product./omega);

convolution2 = infsup(-convolution2,convolution2);

convolution = intersect(convolution1,convolution2);

k = n - (p*(M-1)+1);

convolution = [convolution;zeros(k,1)];

end