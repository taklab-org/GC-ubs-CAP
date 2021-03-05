function g = g_vector_field(x)

if isintval(x) == 1
    g = intval(zeros(3,1));   
    % Fix the parameters
    a = intval('0.3'); c = intval('0.7'); delta = intval(9); w = intval('0.02');
else
    g = zeros(3,1);   
    % Fix the parameters
    a = 0.3; c = 0.7; delta = 9; w = 0.02;
end

x0 = x(1); x1 = x(2); x2 = x(3);

tf0 = x0^3 - (1-x0^2-x1^2-x2^2)*x0;
tf1 = x0^2*x1 + x0^2*x2;
tf2 = x0^2*x2 + (1/delta)*(c*x0^2*x2-x1*(x1-a*x0)*(x0-x1)+w*x0^3);

sum = x0*tf0+x1*tf1+x2*tf2;

%sum = 2*x0^4 - x0^2 + 2*x0^2*x1^2 + (c/delta+2)*x0^2*x2^2 + (1+a/delta)*x0^2*x1*x2 ...
%       + (1/delta)*x1^3*x2 - ((a+1)/delta)*x0*x1^2*x2 + (w/delta)*x0^3*x2;

g(1) = tf0 - x0*sum;
% g0 = 3*x0^3 - x0 + x0*x1^2 + x0*x2^2 - 2*x0^5 - 2*x0^3*x1^2 ... 
%      - (c/delta+2)*x0^3*x2^2 - (1+a/delta)*x0^3*x1*x2 - (1/delta)*x0*x1^3*x2 ...
%      + ((a+1)/delta)*x0^2*x1^2*x2 - (w/delta)*x0^4*x2;
  
g(2) = tf1 - x1*sum;
% g1 = 2*x0^2*x1 + x0^2*x2 - 2*x0^4*x1 - 2*x0^2*x1^3 - (c/delta+2)*x0^2*x1*x2^2 ...
%      - (1+a/delta)*x0^2*x1^2*x2 - (1/delta)*x1^4*x2 + ((a+1)/delta)*x0*x1^3*x2 ...
%      - (w/delta)*x0^3*x1*x2;
 
% g2 = (2+c/delta)*x0^2*x2 + (1/delta)*x1^3 - ((a+1)/delta)*x0*x1^2 + (a/delta)*x0^2*x1 ...
%      + (w/delta)*x0^3 - 2*x0^4*x2 - 2*x0^2*x1^2*x2 - (c/delta+2)*x0^2*x2^3 ...
%      - (1+a/delta)*x0^2*x1*x2^2 - (1/delta)*x1^3*x2^2 + ((a+1)/delta)*x0*x1^2*x2^2 ...
%      - (w/delta)*x0^3*x2^2;
 
g(3) = tf2 - x2*sum;

end