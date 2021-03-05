function [s] = int_general_power_FFT(a)

% Inputs : (p vectors of Fourier coefficients - a1, ..., ap) 
%  ai = (ai_0,ai_1,...,ai_{M-1}) \in C^M, for i=1,...,p
% Assumption: the Fourier coeffcients represent a real function, that is we
% have the relation ai_{-k} = conj(ai_{k}). 

% This code computes a rigorous enclosure of a discrete convolution 
% (a1*...*ap)_k = sum_{k1+...+kp=k, ki in Z} a1_k1 ... ap_kp for k=0,...,p*(M-1)

[M,p] = size(a); % p : power of the nonlinearity

% We make sure that the inputs are powers of 2 %

exponent = ceil( log(2*p*M-1)/log(2) );
M1 = (2^exponent-2*p*M)/2; % Hence 2*p*M+2*M1 is a power of 2 %

a_ext = zeros(2*M-1,p);
ta = intval(zeros(2*p*M+2*M1,p));
tu = intval(zeros(2*p*M+2*M1,p));
tb = ones(2*p*M+2*M1,1);

for i=1:p
    a_ext(:,i) = [flip(conj(a(2:M,i)));a(:,i)];
    ta(:,i) = [intval(zeros((p-1)*M+M1+1,1));intval(a_ext(:,i));intval(zeros((p-1)*M+M1,1))];
    tu(:,i) = verifyfft(ifftshift(ta(:,i)),-1);
    tb = (2*p*M+2*M1)*tb.*tu(:,i);
end

F = fftshift(verifyfft(tb,1));

s = real(F(M1+p*M+1:M1+p*(2*M-1)+1))/(2*p*M+2*M1);

end