function [reduced,extended]=taylorconvolution(varargin)
% aout=taylorconvolution(a1,a2,a3,...)
% Uses fft to compute the convolution of a1, a2, a3, ...
% which each represents a Taylor polynomial.
% It is assumed all inputs have the same size.
% The number of inputs and dimensions is arbitrary.
%
% The input can be (tensors of) floats or intvals
% The first output has the same size as each of the inputs
% The second output is the full extended convolution (all nonzero terms)
% 
% Possible improvements are to use Banach algebra estimates
% and/or reduce sizes by using nested convolutions
% (working with inputs of different sizes)

a=varargin;
K=length(a);
N=size(a{1})-1;
M=length(N); 

% define some indices
for i=1:M
    inda{i}=1:(N(i)+1);
end
for i=1:M
    indq{i}=1:(K*N(i)+1);
end

zeroindices=K*N+1;
if exist('intval','file') && isintval(a{1}(1))
    zeroindices = 2.^nextpow2(zeroindices);
    extendedzeros=intval(zeros(zeroindices));
else
    extendedzeros=zeros(zeroindices);
end

for j=1:K
    % padding 
    aextended=extendedzeros;
    aextended(inda{1:M})=a{j};
    % fft
    Fa{j}=altfftn(aextended);
end

% product
Fq=Fa{1};
for j=2:K
    Fq=Fq.*Fa{j};
end

%ifft
extended=altifftn(Fq);
if exist('intval','file') && isintval(a{1}(1))
    %remove extraneous zeros
    extended=extended(indq{1:M});
end

% the part matching input size
reduced=extended(inda{1:M});

end

function y = altfftn(x)
% computes the FFT of a tensor of any dimensions
% for variables of double or intval type
%
% for intval type the size of the input (in all dimensions) 
% have to be powers of 2

if exist('intval','file') && isintval(x(1))
  n = size(x);
  m = length(n);
  y = x;   
  for dimension = 1:m
      if n(dimension)>1   % no fft if first dimension is singleton! 
          y = reshape(verifyfft(y,1),size(y));
      end
      y = permute(y,[2:m 1]);
  end
else
  y=fftn(x); 
end

end


function y = altifftn(x)
% computes the IFFT of a tensor of any dimensions
% for variables of double or intval type
%
% for intval type the size of the input (in all dimensions) 
% have to be powers of 2

if exist('intval','file') && isintval(x(1))
  n = size(x);
  m = length(n);
  y = x;   
  for dimension = 1:m
      if n(dimension)>1   % no ifft if first dimension is singleton! 
          y = reshape(verifyfft(y,-1),size(y));
      end
      y = permute(y,[2:m 1]);
  end
else
  y=ifftn(x); 
end

end


